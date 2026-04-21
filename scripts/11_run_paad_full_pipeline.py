#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import math
import shutil
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import GroupKFold, KFold

try:
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
except Exception:
    pass

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import (  # noqa: E402
    NAME_NORMALIZER,
    admet_assessment,
    build_numeric_base,
    build_smiles_svd,
    build_strong_context,
    load_config,
    make_model,
    pipeline_paths,
    read_json,
    regression_metrics,
    safe_pearson,
    safe_spearman,
    split_genes,
    write_json,
    write_table,
)


def _load_script(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load helper script: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


SRC = _load_script("thyroid_source_build_helpers", REPO_ROOT / "scripts" / "02_build_model_ready_from_thyroid_raw.py")
DL = _load_script("thyroid_dl_helpers", REPO_ROOT / "scripts" / "07_train_additional_dl_models.py")

TCGA_CODE = "PAAD"
INPUT_NAME = "numeric_strong_context_smiles_pan_lincs"


def log(message: str) -> None:
    print(message, flush=True)


def _fmt(value: Any) -> str:
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.4f}"
    except Exception:
        return str(value)


def _md_table(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    if df is None or df.empty:
        return "_No rows._"
    sub = df[[c for c in cols if c in df.columns]].head(n).copy()
    if sub.empty:
        return "_No matching columns._"
    lines = ["| " + " | ".join(sub.columns) + " |", "| " + " | ".join(["---"] * len(sub.columns)) + " |"]
    for row in sub.itertuples(index=False):
        lines.append("| " + " | ".join(_fmt(v) for v in row) + " |")
    return "\n".join(lines)


def priority_genes(cfg: dict[str, Any]) -> set[str]:
    genes = {
        g.upper()
        for g in cfg.get("paad_biology_terms", [])
        if isinstance(g, str) and g.replace("-", "").replace(".", "").isalnum() and len(g) <= 16
    }
    genes.update(
        {
            "KRAS",
            "TP53",
            "CDKN2A",
            "SMAD4",
            "BRCA1",
            "BRCA2",
            "PALB2",
            "ATM",
            "ATR",
            "CHEK1",
            "CHEK2",
            "PARP1",
            "PARP2",
            "EGFR",
            "ERBB2",
            "ERBB3",
            "PIK3CA",
            "AKT1",
            "MTOR",
            "MAP2K1",
            "MAP2K2",
            "MAPK1",
            "MAPK3",
            "MYC",
            "GATA6",
            "MSLN",
            "TOP1",
            "TOP2A",
            "RRM1",
            "RRM2",
            "TUBB",
            "AURKA",
            "AURKB",
            "WEE1",
        }
    )
    return genes


def ensure_external_sources(cfg: dict[str, Any]) -> dict[str, Path]:
    paths = pipeline_paths(cfg)
    stage = paths.root / "data" / "paad_source_staging"
    tcga_dir = stage / "tcga_paad"
    ct_dir = stage / "clinical_trials"
    tcga_dir.mkdir(parents=True, exist_ok=True)
    ct_dir.mkdir(parents=True, exist_ok=True)

    sources = cfg.get("external_sources", {})
    mapping = {
        "tcga_expression": (sources.get("tcga_expression"), tcga_dir / "TCGA-PAAD.star_tpm.tsv.gz"),
        "tcga_clinical": (sources.get("tcga_clinical"), tcga_dir / "TCGA-PAAD.clinical.tsv.gz"),
        "tcga_survival": (sources.get("tcga_survival"), tcga_dir / "TCGA-PAAD.survival.tsv.gz"),
        "clinical_trials": (
            sources.get("clinical_trials"),
            ct_dir / "clinicaltrials_paad_cancer_drug_20260421.json",
        ),
    }
    out: dict[str, Path] = {}
    for key, (src_raw, dst) in mapping.items():
        src = Path(str(src_raw)) if src_raw else dst
        if not dst.exists():
            if not src.exists():
                raise FileNotFoundError(f"Missing PAAD external source for {key}: {src}")
            shutil.copy2(src, dst)
        out[key] = dst
    return out


def build_gdsc_response_and_base(staging: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gdsc_dir = staging / "gdsc"
    gdsc = pd.read_parquet(gdsc_dir / "gdsc2_annotation_normalized_20260406.parquet")
    paad = gdsc.loc[gdsc["TCGA_DESC"].astype(str).str.upper().eq(TCGA_CODE)].copy()
    paad["sample_id"] = paad["SANGER_MODEL_ID"].astype(str)
    paad["canonical_drug_id"] = paad["DRUG_ID"].astype(str)
    paad["cell_line_name"] = paad["CELL_LINE_NAME"].astype(str)
    paad["drug_name"] = paad["DRUG_NAME"].astype(str)
    paad["disease_label"] = "pancreatic cancer"
    response = paad[
        [
            "sample_id",
            "cell_line_name",
            "disease_label",
            "canonical_drug_id",
            "drug_name",
            "LN_IC50",
            "AUC",
            "Z_SCORE",
        ]
    ].copy()
    response = response.dropna(subset=["LN_IC50"]).drop_duplicates(["sample_id", "canonical_drug_id"])

    cell = pd.read_parquet(gdsc_dir / "gdsc2_cellline_annotation_table_20260406.parquet")
    cell = cell.loc[cell["TCGA_DESC"].astype(str).str.upper().eq(TCGA_CODE)].copy()
    drug_ann = pd.read_parquet(gdsc_dir / "gdsc2_drug_annotation_table_20260406.parquet").copy()
    drug_ann["canonical_drug_id"] = drug_ann["DRUG_ID"].astype(str)
    return response, cell, drug_ann


def build_sample_features(staging: Path, cell: pd.DataFrame, cfg: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    depmap_dir = staging / "depmap"
    model = pd.read_parquet(depmap_dir / "depmap_model_basic_clean_20260406.parquet")
    keep_model_cols = [
        "ModelID",
        "CellLineName",
        "SangerModelID",
        "COSMICID",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "OncotreeCode",
        "Age",
        "Sex",
        "PrimaryOrMetastasis",
        "SampleCollectionSite",
    ]
    model = model[[c for c in keep_model_cols if c in model.columns]].copy()
    cell = cell.merge(model, left_on="SANGER_MODEL_ID", right_on="SangerModelID", how="left")

    cosmic = cell["COSMIC_ID"] if "COSMIC_ID" in cell.columns else cell.get("COSMICID", np.nan)
    subtype = cell.get("OncotreeSubtype", pd.Series([""] * len(cell))).fillna("").astype(str)
    primary = cell.get("OncotreePrimaryDisease", pd.Series([""] * len(cell))).fillna("").astype(str)
    site = cell.get("SampleCollectionSite", pd.Series([""] * len(cell))).fillna("").astype(str)

    sample = pd.DataFrame(
        {
            "sample_id": cell["SANGER_MODEL_ID"].astype(str),
            "cell_line_name": cell["CELL_LINE_NAME"].astype(str),
            "disease_label": subtype.where(subtype.ne(""), "pancreatic cancer"),
            "paad_subtype": subtype,
            "sample_depmap_model_id": cell["ModelID"].astype(str).replace("nan", ""),
            "sample_cosmic_id": pd.to_numeric(cosmic, errors="coerce"),
            "sample_has_depmap_model": cell["ModelID"].notna().astype(int),
            "sample__depmap_age": pd.to_numeric(cell.get("Age"), errors="coerce") if "Age" in cell.columns else np.nan,
            "sample__is_male": cell.get("Sex", pd.Series([""] * len(cell))).astype(str).str.lower().eq("male").astype(int),
            "sample__is_female": cell.get("Sex", pd.Series([""] * len(cell))).astype(str).str.lower().eq("female").astype(int),
        }
    )
    sample["sample__is_pancreatic_adenocarcinoma"] = (
        subtype.str.contains("pancreatic|ductal|adenocarcinoma", case=False, regex=True)
        | primary.str.contains("pancreatic", case=False, regex=True)
    ).astype(int)
    sample["sample__is_exocrine"] = subtype.str.contains("exocrine|ductal|acinar", case=False, regex=True).astype(int)
    sample["sample__is_neuroendocrine"] = subtype.str.contains("neuroendocrine", case=False, regex=True).astype(int)
    sample["sample__is_metastatic"] = (
        cell.get("PrimaryOrMetastasis", pd.Series([""] * len(cell))).astype(str).str.contains("metastatic", case=False, na=False)
        | site.str.contains("metasta", case=False, regex=True)
    ).astype(int)

    priority = priority_genes(cfg)
    crispr = pd.read_parquet(depmap_dir / "depmap_crispr_gene_dependency_basic_clean_20260406.parquet")
    crispr_limit = int(cfg.get("sample_feature_enhancement", {}).get("crispr_feature_limit", 1024))
    selected = SRC.top_variance_columns(crispr, "ModelID", crispr_limit, priority)
    crispr_small = crispr[["ModelID"] + selected].copy()
    crispr_small = crispr_small.rename(
        columns={col: f"sample__depmap_dependency__{SRC.clean_gene_symbol(col)}" for col in selected}
    )

    sample_model = cell[["SANGER_MODEL_ID", "ModelID"]].rename(columns={"SANGER_MODEL_ID": "sample_id"}).copy()
    sample_model["sample_id"] = sample_model["sample_id"].astype(str)
    sample = sample.merge(sample_model, on="sample_id", how="left")
    sample = sample.merge(crispr_small, on="ModelID", how="left")
    sample["sample_has_crispr"] = sample["ModelID"].isin(set(crispr["ModelID"])).astype(int)

    model_ids = set(sample["ModelID"].dropna().astype(str))
    enh = cfg.get("sample_feature_enhancement", {})
    expr_small, expr_qc = SRC.read_depmap_wide_omics(
        depmap_dir / "OmicsExpressionProteinCodingGenesTPMLogp1_24Q2.csv",
        model_ids,
        "expr",
        int(enh.get("expression_feature_limit", 256)),
        priority,
    )
    cnv_small, cnv_qc = SRC.read_depmap_wide_omics(
        depmap_dir / "OmicsCNGene_24Q2.csv",
        model_ids,
        "cnv",
        int(enh.get("cnv_feature_limit", 256)),
        priority,
    )
    mut_small, mut_qc = SRC.read_depmap_mutation_features(
        depmap_dir / "OmicsSomaticMutations_24Q2.csv",
        model_ids,
        priority,
        int(enh.get("mutation_top_gene_limit", 96)),
    )
    sample = sample.merge(expr_small, on="ModelID", how="left")
    sample = sample.merge(cnv_small, on="ModelID", how="left")
    sample = sample.merge(mut_small, on="ModelID", how="left")
    expr_cols = [c for c in sample.columns if c.startswith("sample__depmap_expr__")]
    cnv_cols = [c for c in sample.columns if c.startswith("sample__depmap_cnv__")]
    mut_cols = [c for c in sample.columns if c.startswith("sample__depmap_mut__")]
    sample["sample_has_depmap_expression"] = sample[expr_cols].notna().any(axis=1).astype(int) if expr_cols else 0
    sample["sample_has_depmap_cnv"] = sample[cnv_cols].notna().any(axis=1).astype(int) if cnv_cols else 0
    if "sample_has_depmap_mutation" not in sample.columns:
        sample["sample_has_depmap_mutation"] = 0
    sample["sample_has_non_crispr_omics_fallback"] = (
        sample[["sample_has_depmap_expression", "sample_has_depmap_cnv", "sample_has_depmap_mutation"]].sum(axis=1).gt(0).astype(int)
    )
    sample = sample.drop(columns=["ModelID"], errors="ignore")

    qc = {
        "paad_cell_lines": int(cell["SANGER_MODEL_ID"].nunique()),
        "depmap_model_matches": int(cell["ModelID"].notna().sum()),
        "crispr_feature_rows": int(sample["sample_has_crispr"].sum()),
        "selected_crispr_features": int(len(selected)),
        **expr_qc,
        **cnv_qc,
        **mut_qc,
        "selected_non_crispr_omics_features": int(len(expr_cols) + len(cnv_cols) + len(mut_cols)),
    }
    return sample, qc


def build_drug_features(staging: Path, drug_ann: pd.DataFrame, cfg: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    refs = SRC.load_drug_reference_maps(staging)
    known = {SRC.norm_name(x) for x in cfg.get("known_paad_drugs", [])}
    ct_path = Path(cfg.get("external_sources", {}).get("clinical_trials", ""))
    if not ct_path.exists():
        ct_path = REPO_ROOT / "data" / "paad_source_staging" / "clinical_trials" / "clinicaltrials_paad_cancer_drug_20260421.json"
    clinical_trials_text = ct_path.read_text(encoding="utf-8", errors="ignore").lower() if ct_path.exists() else ""

    bits_n = int(cfg.get("drug_feature_enhancement", {}).get("morgan_bits", 256))
    records = []
    for row in drug_ann.itertuples(index=False):
        drug_id = str(row.canonical_drug_id)
        drug_name = str(row.DRUG_NAME)
        norm = SRC.norm_name(drug_name)
        smiles = ""
        source = "missing"
        drugbank_id = ""
        chembl_id = ""
        indication = ""
        mechanism = ""
        if "catalog_by_id" in refs and drug_id in refs["catalog_by_id"].index:
            cat_row = refs["catalog_by_id"].loc[drug_id]
            if int(cat_row.get("has_smiles", 0)) == 1 and pd.notna(cat_row.get("canonical_smiles")):
                smiles = str(cat_row["canonical_smiles"])
                source = "gdsc_catalog"
        if norm in refs["drugbank_by_norm"].index:
            db_row = refs["drugbank_by_norm"].loc[norm]
            drugbank_id = str(db_row.get("drugbank_id", ""))
            indication = "" if pd.isna(db_row.get("indication")) else str(db_row.get("indication"))
            mechanism = "" if pd.isna(db_row.get("mechanism_of_action")) else str(db_row.get("mechanism_of_action"))
            if not smiles and pd.notna(db_row.get("smiles")):
                smiles = str(db_row["smiles"])
                source = "drugbank_name"
        if norm in refs["chembl_by_norm"].index:
            ch_row = refs["chembl_by_norm"].loc[norm]
            chembl_id = str(ch_row.get("chembl_id", ""))
            if not smiles and pd.notna(ch_row.get("canonical_smiles")):
                smiles = str(ch_row["canonical_smiles"])
                source = "chembl_pref_name"
        can, ok = SRC.canonicalize_smiles(smiles)
        db_targets = refs["drugbank_targets"].get(drugbank_id, "")
        target_genes = SRC.parse_targets(db_targets, row.PUTATIVE_TARGET_NORMALIZED, row.PUTATIVE_TARGET)
        mentioned = bool(drug_name.lower() and drug_name.lower() in clinical_trials_text)
        indication_low = f"{indication} {mechanism}".lower()
        if norm in known:
            classification = "paad_standard_or_known"
        elif "pancrea" in indication_low or mentioned:
            classification = "paad_trial_or_indication"
        else:
            classification = "screened_candidate"
        bridge_bits = [ok, bool(target_genes), norm in refs["drugbank_by_norm"].index, norm in refs["chembl_by_norm"].index]
        bridge_strength = "multi_source" if sum(bridge_bits) >= 3 else "single_source" if sum(bridge_bits) >= 1 else "weak"
        rec = {
            "canonical_drug_id": drug_id,
            "drug_name": drug_name,
            "canonical_smiles": can,
            "smiles_source": source,
            "smiles_parse_ok": int(ok),
            "drugbank_id": drugbank_id,
            "chembl_id": chembl_id,
            "target_genes": target_genes,
            "PATHWAY_NAME_NORMALIZED": getattr(row, "PATHWAY_NAME_NORMALIZED", ""),
            "classification": classification,
            "drug_bridge_strength": bridge_strength,
            "stage3_resolution_status": "resolved" if ok and target_genes else "target_partial" if target_genes else "unresolved",
            "TCGA_DESC": TCGA_CODE,
            "drug__has_smiles": int(ok),
            "drug__has_target": int(bool(target_genes)),
            "drug__is_known_paad_control": int(norm in known),
            "drug__is_clinical_trial_mentioned": int(mentioned),
        }
        rec.update(SRC.smiles_descriptors(can))
        for i, bit in enumerate(SRC.morgan_bits(can, bits_n)):
            rec[f"drug_morgan_{i:04d}"] = bit
        records.append(rec)

    drugs = pd.DataFrame(records)
    lincs_path = staging / "lincs" / "lincs_drug_signature_pancancer_20260421.parquet"
    if not lincs_path.exists():
        lincs_path = staging / "lincs" / "lincs_drug_signature_normalized.parquet"
    lincs_qc = {"lincs_rows": 0, "selected_lincs_features": 0, "drug_lincs_matches": 0}
    if lincs_path.exists():
        lincs = pd.read_parquet(lincs_path)
        lincs["canonical_drug_id"] = lincs["canonical_drug_id"].astype(str)
        lincs_limit = int(cfg.get("drug_feature_enhancement", {}).get("lincs_feature_limit", 512))
        selected = SRC.top_variance_columns(lincs, "canonical_drug_id", lincs_limit, set())
        rename = {}
        for c in selected:
            if str(c).startswith("drug__lincs__"):
                rename[c] = c
            else:
                rename[c] = f"drug__lincs__{SRC.clean_gene_symbol(str(c).replace('crispr__', ''))}"
        lincs_small = lincs[["canonical_drug_id"] + selected].rename(columns=rename)
        lincs_small["drug__has_lincs_signature"] = 1
        drugs = drugs.merge(lincs_small, on="canonical_drug_id", how="left")
        drugs["drug__has_lincs_signature"] = drugs["drug__has_lincs_signature"].fillna(0).astype(int)
        lincs_feature_cols = [c for c in drugs.columns if c.startswith("drug__lincs__")]
        lincs_qc = {
            "lincs_source": str(lincs_path),
            "lincs_rows": int(len(lincs)),
            "selected_lincs_features": int(len(lincs_feature_cols)),
            "drug_lincs_matches": int(drugs["drug__has_lincs_signature"].sum()),
        }
    else:
        drugs["drug__has_lincs_signature"] = 0

    annotations = drugs[
        [
            "canonical_drug_id",
            "drug_name",
            "target_genes",
            "PATHWAY_NAME_NORMALIZED",
            "classification",
            "drug_bridge_strength",
            "stage3_resolution_status",
            "TCGA_DESC",
            "smiles_source",
            "drugbank_id",
            "chembl_id",
        ]
    ].copy()
    qc = {
        "drug_count": int(len(drugs)),
        "smiles_present": int(drugs["canonical_smiles"].astype(str).str.len().gt(0).sum()),
        "smiles_parse_ok": int(drugs["smiles_parse_ok"].sum()),
        "target_gene_present": int(drugs["target_genes"].astype(str).str.len().gt(0).sum()),
        "classification_counts": drugs["classification"].value_counts(dropna=False).to_dict(),
        "smiles_source_counts": drugs["smiles_source"].value_counts(dropna=False).to_dict(),
        "morgan_bits": int(bits_n),
        **lincs_qc,
    }
    return drugs, annotations, qc


def build_external_validation_files(
    cfg: dict[str, Any],
    external_sources: dict[str, Path],
    target_genes: set[str],
) -> dict[str, Any]:
    paths = pipeline_paths(cfg)
    output_dir = paths.external_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    gencode = paths.root / "data" / "source_staging" / "gencode" / "gencode.v36.annotation.gtf.gz"
    expression = external_sources["tcga_expression"]
    survival = external_sources["tcga_survival"]
    clinical = external_sources["tcga_clinical"]

    gene_map: dict[str, str] = {}
    with gzip.open(gencode, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            gid = SRC.re.search(r'gene_id "([^"]+)"', parts[8])
            gname = SRC.re.search(r'gene_name "([^"]+)"', parts[8])
            if gid and gname:
                gene_map[gid.group(1).split(".")[0]] = gname.group(1).upper()

    rows = []
    target_upper = {g.upper() for g in target_genes if g}
    with gzip.open(expression, "rt", encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            symbol = gene_map.get(fields[0].split(".")[0], "")
            if symbol and symbol in target_upper:
                rows.append([symbol] + fields[1:])
    expr_df = pd.DataFrame(rows, columns=["Hugo_Symbol"] + header[1:])
    if not expr_df.empty:
        value_cols = [c for c in expr_df.columns if c != "Hugo_Symbol"]
        expr_df[value_cols] = expr_df[value_cols].apply(pd.to_numeric, errors="coerce")
        expr_df = expr_df.groupby("Hugo_Symbol", as_index=False)[value_cols].mean()
    write_table(expr_df, output_dir / "paad_expression.csv")

    surv = pd.read_csv(survival, sep="\t", compression="gzip")
    clin = pd.read_csv(clinical, sep="\t", compression="gzip")
    clinical_out = surv.rename(columns={"sample": "PATIENT_ID", "OS.time": "OS_DAYS", "OS": "OS_STATUS"}).copy()
    clinical_out["OS_MONTHS"] = pd.to_numeric(clinical_out["OS_DAYS"], errors="coerce") / 30.4375
    if "sample" in clin.columns:
        clin_small_cols = [
            "sample",
            "gender.demographic",
            "age_at_index.demographic",
            "ajcc_pathologic_stage.diagnoses",
            "primary_diagnosis.diagnoses",
            "tumor_grade.diagnoses",
        ]
        clin_small = clin[[c for c in clin_small_cols if c in clin.columns]].copy()
        clin_small = clin_small.rename(columns={"sample": "PATIENT_ID"})
        clinical_out = clinical_out.merge(clin_small, on="PATIENT_ID", how="left")
    write_table(clinical_out, output_dir / "paad_clinical.csv")
    return {
        "target_gene_universe": int(len(target_upper)),
        "external_expression_genes_written": int(expr_df["Hugo_Symbol"].nunique()) if not expr_df.empty else 0,
        "external_expression_samples": int(max(0, expr_df.shape[1] - 1)),
        "clinical_rows": int(len(clinical_out)),
    }


def build_model_ready(cfg: dict[str, Any], staging: Path) -> dict[str, Any]:
    paths = pipeline_paths(cfg)
    external_sources = ensure_external_sources(cfg)
    response, cell, drug_ann = build_gdsc_response_and_base(staging)
    sample_features, sample_qc = build_sample_features(staging, cell, cfg)
    drug_features, drug_annotations, drug_qc = build_drug_features(staging, drug_ann, cfg)
    response, drug_features, drug_annotations, filter_qc = SRC.apply_primary_candidate_filters(
        response,
        drug_features,
        drug_annotations,
        cfg,
    )
    target_genes = set()
    for value in drug_features["target_genes"].dropna():
        target_genes.update([g for g in str(value).split(";") if g])
    target_genes.update(priority_genes(cfg))
    external_qc = build_external_validation_files(cfg, external_sources, target_genes)
    admet_qc = SRC.build_admet_csvs(staging, paths.admet_source_dir)

    write_table(response, paths.raw_dir / "paad_response_pairs.csv")
    write_table(sample_features, paths.raw_dir / "sample_features.csv")
    write_table(drug_features, paths.raw_dir / "drug_features.csv")
    write_table(drug_annotations, paths.raw_dir / "drug_annotations.csv")

    rows = response.merge(sample_features.drop_duplicates("sample_id"), on="sample_id", how="left", suffixes=("", "_sample"))
    rows = rows.merge(drug_features.drop_duplicates("canonical_drug_id"), on="canonical_drug_id", how="left", suffixes=("", "_drug"))
    duplicate_ann_cols = [c for c in drug_annotations.columns if c in rows.columns and c != "canonical_drug_id"]
    rows = rows.merge(drug_annotations.drop(columns=duplicate_ann_cols).drop_duplicates("canonical_drug_id"), on="canonical_drug_id", how="left")
    rows["TCGA_DESC"] = TCGA_CODE
    rows["disease_label_normalized"] = TCGA_CODE
    rows["paad_subtype"] = rows.get("paad_subtype", "").fillna("").astype(str)
    write_table(rows, paths.processed_dir / "paad_response_pairs.csv", also_parquet=True)
    write_table(rows, paths.processed_dir / "row_metadata.csv", also_parquet=True)

    qc = {
        "step": "paad_source_to_model_ready",
        "source_staging_dir": str(staging),
        "response_rows": int(len(response)),
        "response_cell_lines": int(response["sample_id"].nunique()),
        "response_drugs": int(response["canonical_drug_id"].nunique()),
        "row_metadata_shape": [int(rows.shape[0]), int(rows.shape[1])],
        "sample_features_shape": [int(sample_features.shape[0]), int(sample_features.shape[1])],
        "drug_features_shape": [int(drug_features.shape[0]), int(drug_features.shape[1])],
        "drug_annotations_shape": [int(drug_annotations.shape[0]), int(drug_annotations.shape[1])],
        "sample_qc": sample_qc,
        "drug_qc": drug_qc,
        "primary_filter_qc": filter_qc,
        "external_qc": external_qc,
        "admet_qc": admet_qc,
    }
    write_json(paths.reports_dir / "qc_paad_source_to_model_ready_20260421.json", qc)
    return qc


def write_feature_matrices(cfg: dict[str, Any]) -> dict[str, Any]:
    paths = pipeline_paths(cfg)
    rows = pd.read_csv(paths.processed_dir / "row_metadata.csv")
    X_numeric, y, numeric_names = build_numeric_base(cfg, rows)
    X_smiles, smiles_names = build_smiles_svd(cfg, rows)
    X_context, context_names = build_strong_context(cfg, rows)

    lincs_mask = np.array(["lincs" in str(name).lower() for name in numeric_names], dtype=bool)
    X_numeric_no_lincs = X_numeric[:, ~lincs_mask]
    numeric_no_lincs_names = [name for name, is_lincs in zip(numeric_names, lincs_mask) if not is_lincs]

    matrices = {
        "numeric_pan_lincs": (X_numeric, numeric_names),
        "numeric_no_lincs": (X_numeric_no_lincs, numeric_no_lincs_names),
        "numeric_smiles_pan_lincs": (np.concatenate([X_numeric, X_smiles], axis=1), numeric_names + smiles_names),
        "numeric_smiles_no_lincs": (
            np.concatenate([X_numeric_no_lincs, X_smiles], axis=1),
            numeric_no_lincs_names + smiles_names,
        ),
        "numeric_strong_context_smiles_pan_lincs": (
            np.concatenate([X_numeric, X_context, X_smiles], axis=1),
            numeric_names + context_names + smiles_names,
        ),
        "numeric_strong_context_smiles_no_lincs": (
            np.concatenate([X_numeric_no_lincs, X_context, X_smiles], axis=1),
            numeric_no_lincs_names + context_names + smiles_names,
        ),
    }
    summary = {
        "row_count": int(len(y)),
        "label_count": int(len(y)),
        "lincs_numeric_features": int(lincs_mask.sum()),
        "matrices": {},
    }
    for name, (matrix, features) in matrices.items():
        matrix = matrix.astype(np.float32)
        np.save(paths.processed_dir / f"X_{name}.npy", matrix)
        write_json(paths.processed_dir / f"{name}_feature_names.json", features)
        summary["matrices"][name] = {"shape": [int(matrix.shape[0]), int(matrix.shape[1])], "feature_count": int(len(features))}
    write_json(paths.reports_dir / "qc_paad_feature_matrices_20260421.json", summary)
    return summary


def _splits(cfg: dict[str, Any], X: np.ndarray, y: np.ndarray, rows: pd.DataFrame, cv: str) -> list[tuple[np.ndarray, np.ndarray]]:
    seed = int(cfg["training"].get("random_state", cfg["project"].get("seed", 42)))
    n_splits = int(cfg["training"].get("n_splits", 4))
    if cv == "random4":
        return list(KFold(n_splits=min(n_splits, len(y)), shuffle=True, random_state=seed).split(X))
    if cv == "groupcv4_drug":
        groups = rows["canonical_drug_id"].astype(str).to_numpy()
        return list(GroupKFold(n_splits=min(n_splits, len(np.unique(groups)))).split(X, y, groups))
    raise ValueError(f"Unknown cv: {cv}")


def train_ml_cv(cfg: dict[str, Any], input_name: str, cv: str, force: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = pipeline_paths(cfg)
    X = np.load(paths.processed_dir / f"X_{input_name}.npy").astype(np.float32)
    y = np.load(paths.processed_dir / "y_train.npy").astype(np.float32)
    rows = pd.read_csv(paths.processed_dir / "row_metadata.csv")
    out_root = paths.results_dir / "ml" / cv / input_name
    summary_path = out_root / "ml_metrics_summary.csv"
    fold_path = out_root / "ml_fold_metrics.csv"
    if summary_path.exists() and fold_path.exists() and not force:
        return pd.read_csv(summary_path), pd.read_csv(fold_path)

    out_root.mkdir(parents=True, exist_ok=True)
    oof_dir = out_root / "oof"
    fold_dir = out_root / "folds"
    oof_dir.mkdir(parents=True, exist_ok=True)
    fold_dir.mkdir(parents=True, exist_ok=True)
    seed = int(cfg["training"].get("random_state", cfg["project"].get("seed", 42)))
    splits = _splits(cfg, X, y, rows, cv)
    models = cfg["training"].get("models", [])
    summary_rows = []
    fold_rows = []
    for model_name in models:
        if make_model(model_name, seed, len(y)) is None:
            log(f"[ml:{cv}] skip {model_name}: dependency unavailable")
            continue
        log(f"[ml:{cv}] {model_name} start input={input_name}")
        model_start = time.time()
        oof = np.zeros(len(y), dtype=np.float32)
        model_folds = []
        for fold, (train_idx, valid_idx) in enumerate(splits, start=1):
            fold_start = time.time()
            model = make_model(model_name, seed + fold, len(y))
            model.fit(X[train_idx], y[train_idx])
            train_pred = np.asarray(model.predict(X[train_idx]), dtype=float)
            valid_pred = np.asarray(model.predict(X[valid_idx]), dtype=float)
            oof[valid_idx] = valid_pred.astype(np.float32)
            tr = regression_metrics(y[train_idx], train_pred)
            va = regression_metrics(y[valid_idx], valid_pred)
            row = {
                "cv": cv,
                "input_set": input_name,
                "model": model_name,
                "fold": fold,
                "train_rows": int(len(train_idx)),
                "valid_rows": int(len(valid_idx)),
                "valid_drug_groups": int(rows.iloc[valid_idx]["canonical_drug_id"].nunique()),
                "valid_cell_lines": int(rows.iloc[valid_idx]["sample_id"].nunique()),
                "train_spearman": tr["spearman"],
                "valid_spearman": va["spearman"],
                "valid_rmse": va["rmse"],
                "valid_r2": va["r2"],
                "elapsed_sec": float(time.time() - fold_start),
            }
            fold_rows.append(row)
            model_folds.append({"fold": fold, "train": tr, "valid": va, "elapsed_sec": row["elapsed_sec"]})
            log(f"[ml:{cv}] {model_name} fold={fold} spearman={va['spearman']:.4f} rmse={va['rmse']:.4f}")
        metrics = regression_metrics(y, oof)
        gap = float(np.nanmean([f["train"]["spearman"] for f in model_folds]) - metrics["spearman"])
        elapsed = float(time.time() - model_start)
        np.save(oof_dir / f"{model_name}.npy", oof)
        write_json(
            fold_dir / f"{model_name}_{cv}.json",
            {
                "cv": cv,
                "input_set": input_name,
                "model": model_name,
                "fold_metrics": model_folds,
                "oof_metrics": metrics,
                "train_oof_spearman_gap": gap,
                "prediction_variance": float(np.nanvar(oof)),
                "elapsed_sec": elapsed,
            },
        )
        summary_rows.append(
            {
                "cv": cv,
                "input_set": input_name,
                "model": model_name,
                **metrics,
                "train_oof_spearman_gap": gap,
                "prediction_variance": float(np.nanvar(oof)),
                "elapsed_sec": elapsed,
            }
        )
        log(f"[ml:{cv}] {model_name} done spearman={metrics['spearman']:.4f} rmse={metrics['rmse']:.4f}")
    summary = pd.DataFrame(summary_rows).sort_values("spearman", ascending=False, na_position="last")
    folds = pd.DataFrame(fold_rows)
    write_table(summary, summary_path)
    write_table(folds, fold_path)
    return summary, folds


def train_dl_cv(cfg: dict[str, Any], input_name: str, cv: str, force: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = pipeline_paths(cfg)
    X = np.load(paths.processed_dir / f"X_{input_name}.npy").astype(np.float32)
    y = np.load(paths.processed_dir / "y_train.npy").astype(np.float32)
    rows = pd.read_csv(paths.processed_dir / "row_metadata.csv")
    feature_names = read_json(paths.processed_dir / f"{input_name}_feature_names.json")
    blocks = DL.feature_blocks(feature_names)
    device = DL._torch_device("auto")
    out_root = paths.results_dir / "dl" / cv / input_name
    summary_path = out_root / "dl_metrics_summary.csv"
    fold_path = out_root / "dl_fold_metrics.csv"
    if summary_path.exists() and fold_path.exists() and not force:
        return pd.read_csv(summary_path), pd.read_csv(fold_path)

    oof_dir = out_root / "oof"
    fold_dir = out_root / "folds"
    oof_dir.mkdir(parents=True, exist_ok=True)
    fold_dir.mkdir(parents=True, exist_ok=True)
    splits = _splits(cfg, X, y, rows, cv)
    seed = int(cfg["training"].get("random_state", cfg["project"].get("seed", 42)))
    max_epochs = cfg["training"].get("dl_max_epochs")
    summary_rows = []
    fold_rows = []
    for model_name in cfg["training"].get("dl_models", []):
        params = dict(DL.MODEL_DEFAULTS[model_name])
        if max_epochs is not None:
            params["epochs"] = int(max_epochs)
        log(f"[dl:{cv}] {model_name} start input={input_name} device={device} params={params}")
        model_start = time.time()
        oof = np.zeros(len(y), dtype=np.float32)
        model_folds = []
        for fold, (train_idx, valid_idx) in enumerate(splits, start=1):
            pred, fold_result = DL.train_one_fold(
                model_name,
                X,
                y,
                train_idx,
                valid_idx,
                blocks,
                device,
                seed + fold * 1000 + len(model_name),
                params,
            )
            fold_result.fold = fold
            oof[valid_idx] = pred
            row = {
                "cv": cv,
                "input_set": input_name,
                "model": model_name,
                "fold": fold,
                "best_epoch": fold_result.best_epoch,
                "best_valid_rmse": fold_result.best_valid_rmse,
                "train_spearman": fold_result.train_metrics["spearman"],
                "valid_spearman": fold_result.valid_metrics["spearman"],
                "valid_rmse": fold_result.valid_metrics["rmse"],
                "valid_r2": fold_result.valid_metrics["r2"],
                "valid_drug_groups": int(rows.iloc[valid_idx]["canonical_drug_id"].nunique()),
                "valid_cell_lines": int(rows.iloc[valid_idx]["sample_id"].nunique()),
                "elapsed_sec": fold_result.elapsed_sec,
            }
            fold_rows.append(row)
            model_folds.append(
                {
                    "fold": fold,
                    "best_epoch": fold_result.best_epoch,
                    "train": fold_result.train_metrics,
                    "valid": fold_result.valid_metrics,
                    "elapsed_sec": fold_result.elapsed_sec,
                }
            )
            log(
                f"[dl:{cv}] {model_name} fold={fold} epoch={fold_result.best_epoch} "
                f"spearman={fold_result.valid_metrics['spearman']:.4f} rmse={fold_result.valid_metrics['rmse']:.4f}"
            )
        metrics = regression_metrics(y, oof)
        gap = float(np.nanmean([f["train"]["spearman"] for f in model_folds]) - metrics["spearman"])
        elapsed = float(time.time() - model_start)
        np.save(oof_dir / f"{model_name}.npy", oof)
        write_json(
            fold_dir / f"{model_name}_{cv}.json",
            {
                "cv": cv,
                "input_set": input_name,
                "model": model_name,
                "params": params,
                "fold_metrics": model_folds,
                "oof_metrics": metrics,
                "train_oof_spearman_gap": gap,
                "prediction_variance": float(np.nanvar(oof)),
                "elapsed_sec": elapsed,
                "device": str(device),
                "feature_blocks": {k: len(v) for k, v in blocks.items()},
            },
        )
        summary_rows.append(
            {
                "cv": cv,
                "input_set": input_name,
                "model": model_name,
                **metrics,
                "train_oof_spearman_gap": gap,
                "prediction_variance": float(np.nanvar(oof)),
                "elapsed_sec": elapsed,
            }
        )
        log(f"[dl:{cv}] {model_name} done spearman={metrics['spearman']:.4f} rmse={metrics['rmse']:.4f}")
    summary = pd.DataFrame(summary_rows).sort_values("spearman", ascending=False, na_position="last")
    folds = pd.DataFrame(fold_rows)
    write_table(summary, summary_path)
    write_table(folds, fold_path)
    return summary, folds


def build_cv_ensemble(cfg: dict[str, Any], input_name: str, cv: str) -> dict[str, Any]:
    paths = pipeline_paths(cfg)
    y = np.load(paths.processed_dir / "y_train.npy").astype(np.float32)
    rows = pd.read_csv(paths.processed_dir / "row_metadata.csv")
    pred_dirs = [
        ("ML", paths.results_dir / "ml" / cv / input_name / "oof"),
        ("DL", paths.results_dir / "dl" / cv / input_name / "oof"),
    ]
    preds: dict[str, np.ndarray] = {}
    for family, pred_dir in pred_dirs:
        for path in sorted(pred_dir.glob("*.npy")):
            preds[f"{family}_{path.stem}"] = np.load(path).astype(np.float32)
    if not preds:
        raise FileNotFoundError(f"No OOF predictions found for {cv}/{input_name}")

    individual = pd.DataFrame(
        [
            {"cv": cv, "member": name, **regression_metrics(y, pred), "prediction_variance": float(np.nanvar(pred))}
            for name, pred in preds.items()
        ]
    ).sort_values("spearman", ascending=False, na_position="last")
    raw_weights = {row.member: max(0.0, float(row.spearman)) if np.isfinite(row.spearman) else 0.0 for row in individual.itertuples()}
    total = sum(raw_weights.values())
    weighted = {k: (v / total if total > 0 else 1.0 / len(preds)) for k, v in raw_weights.items()}
    equal = {k: 1.0 / len(preds) for k in preds}

    def combine(weights: dict[str, float]) -> np.ndarray:
        out = np.zeros(len(y), dtype=np.float32)
        for name, weight in weights.items():
            out += float(weight) * preds[name]
        return out

    top4_members = individual["member"].head(min(4, len(individual))).tolist()
    top4_raw = {name: raw_weights[name] for name in top4_members}
    top4_total = sum(top4_raw.values())
    top4 = {k: (v / top4_total if top4_total > 0 else 1.0 / len(top4_raw)) for k, v in top4_raw.items()}

    ensemble_preds = {
        "spearman_weighted": combine(weighted),
        "equal_weight": combine(equal),
        "top4_spearman_weighted": combine(top4),
    }
    ensemble = pd.DataFrame(
        [
            {
                "cv": cv,
                "ensemble": name,
                **regression_metrics(y, pred),
                "prediction_variance": float(np.nanvar(pred)),
            }
            for name, pred in ensemble_preds.items()
        ]
    ).sort_values("spearman", ascending=False, na_position="last")
    best_ensemble_name = str(ensemble.iloc[0]["ensemble"])
    best_pred = ensemble_preds[best_ensemble_name]

    diversity_rows = []
    names = list(preds)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            diversity_rows.append(
                {
                    "cv": cv,
                    "model_a": left,
                    "model_b": right,
                    "prediction_pearson_corr": safe_pearson(preds[left], preds[right]),
                    "prediction_spearman_corr": safe_spearman(preds[left], preds[right]),
                    "residual_pearson_corr": safe_pearson(y - preds[left], y - preds[right]),
                    "mean_abs_prediction_gap": float(np.mean(np.abs(preds[left] - preds[right]))),
                }
            )
    diversity = pd.DataFrame(diversity_rows)

    scored = rows.copy()
    scored["ensemble_pred_ln_ic50"] = best_pred
    scored["ensemble_score"] = -best_pred
    group_cols = [
        "canonical_drug_id",
        "drug_name",
        "canonical_smiles",
        "target_genes",
        "PATHWAY_NAME_NORMALIZED",
        "classification",
    ]
    top = (
        scored.groupby(group_cols, dropna=False)
        .agg(
            mean_pred_ln_ic50=("ensemble_pred_ln_ic50", "mean"),
            ensemble_score=("ensemble_score", "mean"),
            screened_rows=("ensemble_score", "size"),
        )
        .reset_index()
        .sort_values("ensemble_score", ascending=False)
    )
    top["_drug_name_norm"] = top["drug_name"].map(lambda x: NAME_NORMALIZER.sub("", str(x).lower()))
    top = top.drop_duplicates("_drug_name_norm", keep="first").drop(columns=["_drug_name_norm"])
    top["rank"] = np.arange(1, len(top) + 1)
    top = top[["rank"] + [c for c in top.columns if c != "rank"]]

    out_dir = paths.results_dir / "ensemble" / cv / input_name
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, pred in ensemble_preds.items():
        np.save(out_dir / f"{name}_ensemble_oof.npy", pred)
    np.save(out_dir / "best_ensemble_oof.npy", best_pred)
    write_table(individual, out_dir / "individual_metrics.csv")
    write_table(ensemble, out_dir / "ensemble_metrics.csv")
    write_table(diversity, out_dir / "ensemble_diversity.csv")
    write_table(top.head(50), out_dir / "ensemble_top50_drugs.csv")
    result = {
        "cv": cv,
        "input_set": input_name,
        "members": list(preds),
        "best_ensemble": best_ensemble_name,
        "weights": {"spearman_weighted": weighted, "equal_weight": equal, "top4_spearman_weighted": top4},
        "individual_metrics": individual.to_dict(orient="records"),
        "ensemble_metrics": ensemble.to_dict(orient="records"),
        "diversity_rows": int(len(diversity)),
        "top50_path": str(out_dir / "ensemble_top50_drugs.csv"),
    }
    write_json(out_dir / "ensemble_results.json", result)
    return result


def survival_signal(values: pd.Series | None, clinical: pd.DataFrame | None) -> dict[str, Any]:
    if values is None or clinical is None or clinical.empty:
        return {"p_value": np.nan, "direction": "no_clinical_data"}
    endpoint = "OS_MONTHS" if "OS_MONTHS" in clinical.columns else None
    patient_col = "PATIENT_ID" if "PATIENT_ID" in clinical.columns else clinical.columns[0]
    if endpoint is None:
        return {"p_value": np.nan, "direction": "no_endpoint"}
    merged = clinical[[patient_col, endpoint]].copy()
    merged[endpoint] = pd.to_numeric(merged[endpoint], errors="coerce")
    expr = values.rename("target_expression").reset_index().rename(columns={"index": patient_col})
    merged = merged.merge(expr, on=patient_col, how="inner").dropna()
    if len(merged) < 8:
        return {"p_value": np.nan, "direction": "too_few_patients"}
    high = merged.loc[merged["target_expression"] >= merged["target_expression"].median(), endpoint]
    low = merged.loc[merged["target_expression"] < merged["target_expression"].median(), endpoint]
    if len(high) < 3 or len(low) < 3:
        return {"p_value": np.nan, "direction": "too_few_patients"}
    p = float(stats.mannwhitneyu(high, low, alternative="two-sided").pvalue)
    direction = "high_expression_longer_survival" if high.median() > low.median() else "high_expression_shorter_survival"
    return {"p_value": p, "direction": direction}


def clinical_trial_counts(cfg: dict[str, Any], drug_names: list[str]) -> dict[str, int]:
    ct_path = Path(cfg.get("external_sources", {}).get("clinical_trials", ""))
    if not ct_path.exists():
        ct_path = REPO_ROOT / "data" / "paad_source_staging" / "clinical_trials" / "clinicaltrials_paad_cancer_drug_20260421.json"
    text = ct_path.read_text(encoding="utf-8", errors="ignore").lower() if ct_path.exists() else ""
    return {name: int(text.count(str(name).lower())) if str(name).strip() else 0 for name in drug_names}


def prism_validation(staging: Path, candidates: pd.DataFrame) -> pd.DataFrame:
    depmap = staging / "depmap"
    model = pd.read_parquet(depmap / "depmap_model_basic_clean_20260406.parquet")
    pancreatic_models = set(
        model.loc[
            model["OncotreeLineage"].astype(str).str.contains("pancreas", case=False, na=False)
            | model["OncotreePrimaryDisease"].astype(str).str.contains("pancrea", case=False, na=False)
            | model["OncotreeSubtype"].astype(str).str.contains("pancrea", case=False, na=False),
            "ModelID",
        ].astype(str)
    )
    matrix = pd.read_parquet(depmap / "depmap_repurposing_response_matrix_20260406.parquet")
    catalog = pd.read_parquet(depmap / "depmap_repurposing_compound_catalog_table_20260406.parquet")
    catalog["_norm"] = catalog["DrugName"].map(SRC.norm_name)
    catalog["_syn_norm"] = catalog["Synonyms"].fillna("").map(SRC.norm_name)
    panc_cols = [c for c in matrix.columns if c in pancreatic_models]
    rows = []
    for cand in candidates.itertuples():
        norm = SRC.norm_name(cand.drug_name)
        hits = catalog.loc[catalog["_norm"].eq(norm) | catalog["_syn_norm"].str.contains(norm, regex=False, na=False)].copy()
        compound_ids = hits["CompoundID"].astype(str).dropna().unique().tolist()
        mat = matrix.loc[matrix["CompoundID"].astype(str).isin(compound_ids), panc_cols] if compound_ids and panc_cols else pd.DataFrame()
        values = pd.to_numeric(mat.drop(columns=["CompoundID"], errors="ignore").stack(), errors="coerce") if not mat.empty else pd.Series(dtype=float)
        rows.append(
            {
                "canonical_drug_id": cand.canonical_drug_id,
                "drug_name": cand.drug_name,
                "prism_compound_ids": ";".join(compound_ids[:5]),
                "prism_pancreatic_cell_lines": int(len(panc_cols)),
                "prism_observations": int(values.notna().sum()),
                "prism_mean_log2fc": float(values.mean()) if values.notna().any() else np.nan,
                "prism_sensitive_fraction_lt_minus1": float((values < -1.0).mean()) if values.notna().any() else np.nan,
                "prism_has_evidence": bool(values.notna().any()),
            }
        )
    return pd.DataFrame(rows)


def opentargets_evidence(staging: Path, candidates: pd.DataFrame) -> pd.DataFrame:
    ot = staging / "opentargets"
    disease = pd.read_parquet(ot / "opentargets_disease_basic_20260406.parquet", columns=["id", "name", "synonyms"])
    disease_text = disease["name"].astype(str).str.lower() + " " + disease["synonyms"].astype(str).str.lower()
    disease_ids = set(disease.loc[disease_text.str.contains("pancrea", na=False), "id"].astype(str))
    target = pd.read_parquet(ot / "opentargets_target_basic_20260406.parquet", columns=["id", "approvedSymbol"])
    symbol_to_id = dict(zip(target["approvedSymbol"].astype(str).str.upper(), target["id"].astype(str)))
    wanted_target_ids = set()
    drug_targets = {}
    for row in candidates.itertuples():
        genes = split_genes(getattr(row, "target_genes", ""))
        ids = [symbol_to_id[g] for g in genes if g in symbol_to_id]
        drug_targets[row.canonical_drug_id] = ids
        wanted_target_ids.update(ids)
    if not disease_ids or not wanted_target_ids:
        return pd.DataFrame(
            {
                "canonical_drug_id": candidates["canonical_drug_id"],
                "drug_name": candidates["drug_name"],
                "opentargets_max_score": np.nan,
                "opentargets_evidence_count": 0,
            }
        )
    assoc = pd.read_parquet(ot / "opentargets_association_overall_direct_basic_20260406.parquet")
    assoc = assoc.loc[assoc["diseaseId"].astype(str).isin(disease_ids) & assoc["targetId"].astype(str).isin(wanted_target_ids)].copy()
    rows = []
    for row in candidates.itertuples():
        ids = set(drug_targets.get(row.canonical_drug_id, []))
        sub = assoc.loc[assoc["targetId"].astype(str).isin(ids)] if ids else pd.DataFrame()
        rows.append(
            {
                "canonical_drug_id": row.canonical_drug_id,
                "drug_name": row.drug_name,
                "opentargets_target_ids": ";".join(sorted(ids)),
                "opentargets_max_score": float(sub["score"].max()) if not sub.empty else np.nan,
                "opentargets_evidence_count": int(sub["evidenceCount"].sum()) if not sub.empty else 0,
            }
        )
    return pd.DataFrame(rows)


def external_validate_and_rank(cfg: dict[str, Any], staging: Path, cv: str, input_name: str) -> dict[str, Any]:
    paths = pipeline_paths(cfg)
    cand_path = paths.results_dir / "ensemble" / cv / input_name / "ensemble_top50_drugs.csv"
    candidates = pd.read_csv(cand_path).sort_values("rank")
    candidates["canonical_drug_id"] = candidates["canonical_drug_id"].astype(str)
    expr = pd.read_csv(paths.external_dir / "paad_expression.csv")
    clin = pd.read_csv(paths.external_dir / "paad_clinical.csv")
    gene_col = "Hugo_Symbol" if "Hugo_Symbol" in expr.columns else expr.columns[0]
    expr[gene_col] = expr[gene_col].astype(str).str.upper()
    patient_cols = [c for c in expr.columns if c != gene_col]
    expr_values = expr.set_index(gene_col)[patient_cols].apply(pd.to_numeric, errors="coerce")
    global_median = float(np.nanmedian(expr_values.to_numpy())) if expr_values.size else np.nan
    expr_index = set(expr_values.index)
    known = {str(x).upper() for x in cfg.get("known_paad_drugs", [])}
    trial_counts = clinical_trial_counts(cfg, candidates["drug_name"].head(50).astype(str).tolist())
    ext_rows = []
    for row in candidates.head(50).itertuples():
        genes = split_genes(getattr(row, "target_genes", ""))
        matches = [g for g in genes if g in expr_index]
        if matches:
            values = expr_values.loc[matches].mean(axis=0)
            pct = float((values >= global_median).mean()) if np.isfinite(global_median) else np.nan
        else:
            values = None
            pct = np.nan
        surv = survival_signal(values, clin)
        ext_rows.append(
            {
                "rank": int(row.rank),
                "canonical_drug_id": row.canonical_drug_id,
                "drug_name": row.drug_name,
                "target_genes": getattr(row, "target_genes", ""),
                "target_match_genes": ";".join(matches),
                "target_gene_match_rate": float(len(matches) / len(genes)) if genes else np.nan,
                "target_expression_pct": pct,
                "target_expressed": bool(np.isfinite(pct) and pct >= 0.30),
                "survival_p_value": surv["p_value"],
                "survival_direction": surv["direction"],
                "known_paad_control": str(row.drug_name).upper() in known,
                "clinical_trial_mention_count": trial_counts.get(str(row.drug_name), 0),
                "ensemble_score": row.ensemble_score,
                "mean_pred_ln_ic50": row.mean_pred_ln_ic50,
            }
        )
    external = pd.DataFrame(ext_rows)
    external["canonical_drug_id"] = external["canonical_drug_id"].astype(str)
    prism = prism_validation(staging, candidates.head(50))
    prism["canonical_drug_id"] = prism["canonical_drug_id"].astype(str)
    ot = opentargets_evidence(staging, candidates.head(50))
    ot["canonical_drug_id"] = ot["canonical_drug_id"].astype(str)
    external = external.merge(prism, on=["canonical_drug_id", "drug_name"], how="left")
    external = external.merge(ot, on=["canonical_drug_id", "drug_name"], how="left")
    out_ext = paths.external_validation_dir / cv
    out_ext.mkdir(parents=True, exist_ok=True)
    write_table(external.head(15), out_ext / "top15_validated.csv")
    write_table(external, out_ext / "top50_external_validation.csv")

    admet = admet_assessment(cfg, external.head(15))
    admet["canonical_drug_id"] = admet["canonical_drug_id"].astype(str)
    final = external.head(15).merge(
        admet[
            [
                "canonical_drug_id",
                "admet_coverage",
                "toxicity_flags",
                "low_confidence_toxic_signals",
                "admet_category",
            ]
        ],
        on="canonical_drug_id",
        how="left",
    )
    final["model_rank_score"] = (len(final) - np.arange(len(final))) / max(1, len(final))
    final["target_expression_bonus"] = final["target_expressed"].astype(int) * 0.75
    final["survival_bonus"] = (
        pd.to_numeric(final["survival_p_value"], errors="coerce").lt(0.10)
        & final["survival_direction"].astype(str).str.contains("shorter", na=False)
    ).astype(int) * 0.5
    final["clinical_trial_bonus"] = np.minimum(pd.to_numeric(final["clinical_trial_mention_count"], errors="coerce").fillna(0), 5) / 5.0
    final["prism_bonus"] = (
        final["prism_has_evidence"].fillna(False).astype(bool).astype(float)
        * np.maximum(0, -pd.to_numeric(final["prism_mean_log2fc"], errors="coerce").fillna(0))
        .clip(upper=2)
        / 2.0
    )
    final["opentargets_bonus"] = pd.to_numeric(final["opentargets_max_score"], errors="coerce").fillna(0).clip(upper=1.0)
    final["admet_bonus"] = final["admet_category"].map({"Approved": 0.6, "Candidate": 0.4, "Caution": -0.4, "NO_SMILES": -0.8}).fillna(0.0)
    final["knowledge_score"] = (
        final["model_rank_score"]
        + final["target_expression_bonus"]
        + final["survival_bonus"]
        + final["clinical_trial_bonus"]
        + final["prism_bonus"]
        + final["opentargets_bonus"]
        + final["admet_bonus"]
    )
    final["tier"] = np.select(
        [
            (final["knowledge_score"] >= 3.0) & ~final["admet_category"].eq("Caution"),
            final["knowledge_score"] >= 2.1,
            final["admet_category"].eq("Caution"),
        ],
        ["Tier 1", "Tier 2", "Excluded"],
        default="Tier 3",
    )
    final["final_category"] = np.select(
        [
            final["known_paad_control"].fillna(False).astype(bool),
            pd.to_numeric(final["clinical_trial_mention_count"], errors="coerce").fillna(0).gt(0),
        ],
        ["Known pancreatic cancer control/support", "Pancreatic clinical-trial supported candidate"],
        default="Exploratory repurposing candidate",
    )
    tier_order = {"Tier 1": 1, "Tier 2": 2, "Tier 3": 3, "Excluded": 4}
    final["_tier_order"] = final["tier"].map(tier_order).fillna(99)
    final = (
        final.sort_values(["_tier_order", "knowledge_score", "ensemble_score"], ascending=[True, False, False])
        .drop(columns=["_tier_order"])
    )

    out_kg = paths.knowledge_validation_dir / cv
    out_phase = paths.phase5_dir / cv
    out_reports = paths.reports_dir / cv
    out_kg.mkdir(parents=True, exist_ok=True)
    out_phase.mkdir(parents=True, exist_ok=True)
    out_reports.mkdir(parents=True, exist_ok=True)
    write_table(final, out_kg / "validation_summary.csv")
    write_table(final, out_phase / "final_comprehensive_candidates.csv")
    write_table(final.loc[final["tier"].eq("Tier 1")], out_phase / "tier1_high_confidence.csv")
    qc = {
        "cv": cv,
        "input_set": input_name,
        "candidate_source": str(cand_path),
        "expression_gene_count": int(expr_values.shape[0]),
        "patient_count": int(expr_values.shape[1]),
        "target_expressed_count_top15": int(external.head(15)["target_expressed"].sum()),
        "clinical_trial_supported_top15": int(pd.to_numeric(external.head(15)["clinical_trial_mention_count"], errors="coerce").fillna(0).gt(0).sum()),
        "prism_evidence_top15": int(external.head(15)["prism_has_evidence"].fillna(False).astype(bool).sum()),
        "tier_counts": final["tier"].value_counts(dropna=False).to_dict(),
        "category_counts": final["final_category"].value_counts(dropna=False).to_dict(),
    }
    write_json(out_reports / "qc_paad_external_admet_knowledge_20260421.json", qc)
    return {"external": external, "final": final, "qc": qc}


def write_final_report(cfg: dict[str, Any], input_name: str, run_summary: dict[str, Any]) -> Path:
    paths = pipeline_paths(cfg)
    report = paths.root / "docs" / "PAAD_FULL_PIPELINE_20260421.md"
    random_ens = pd.read_csv(paths.results_dir / "ensemble" / "random4" / input_name / "ensemble_metrics.csv")
    group_ens = pd.read_csv(paths.results_dir / "ensemble" / "groupcv4_drug" / input_name / "ensemble_metrics.csv")
    random_ind = pd.read_csv(paths.results_dir / "ensemble" / "random4" / input_name / "individual_metrics.csv")
    group_ind = pd.read_csv(paths.results_dir / "ensemble" / "groupcv4_drug" / input_name / "individual_metrics.csv")
    top = pd.read_csv(paths.results_dir / "ensemble" / "groupcv4_drug" / input_name / "ensemble_top50_drugs.csv")
    external = pd.read_csv(paths.external_validation_dir / "groupcv4_drug" / "top50_external_validation.csv")
    final = pd.read_csv(paths.phase5_dir / "groupcv4_drug" / "final_comprehensive_candidates.csv")
    source_qc = read_json(paths.reports_dir / "qc_paad_source_to_model_ready_20260421.json")
    feature_qc = read_json(paths.reports_dir / "qc_paad_feature_matrices_20260421.json")
    text = f"""# PAAD Full Drug Repurposing Pipeline - 2026-04-21

## Run scope

- Cancer type: pancreatic cancer / TCGA-PAAD
- Primary input set: `{input_name}`
- CV checks: random sample 4-fold and drug GroupCV 4-fold
- Model families: ML, DL, and ML+DL ensembles
- Training label: GDSC2 PAAD `LN_IC50`
- Source prefix: `s3://say2-4team/PAAD_raw/`

## Input data

- Response rows after SMILES filter: `{source_qc.get("response_rows")}`
- PAAD cell lines: `{source_qc.get("response_cell_lines")}`
- Drugs: `{source_qc.get("response_drugs")}`
- Feature matrix shape: `{feature_qc.get("matrices", {}).get(input_name, {}).get("shape")}`
- LINCS numeric features: `{feature_qc.get("lincs_numeric_features")}`

## Random sample 4-fold ensemble

{_md_table(random_ens, ["ensemble", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20"])}

## Drug GroupCV 4-fold ensemble

{_md_table(group_ens, ["ensemble", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20"])}

## Random sample member models

{_md_table(random_ind, ["member", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20"], n=20)}

## Drug GroupCV member models

{_md_table(group_ind, ["member", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20"], n=20)}

## GroupCV ensemble Top 15

{_md_table(top, ["rank", "drug_name", "ensemble_score", "mean_pred_ln_ic50", "screened_rows", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"], n=15)}

## External validation snapshot

{_md_table(external, ["rank", "drug_name", "target_match_genes", "target_expression_pct", "target_expressed", "survival_p_value", "clinical_trial_mention_count", "prism_mean_log2fc", "opentargets_max_score"], n=15)}

## Final tiered candidates

{_md_table(final, ["drug_name", "tier", "knowledge_score", "final_category", "admet_category", "ensemble_score", "target_expressed", "clinical_trial_mention_count", "prism_mean_log2fc"], n=15)}

## Key outputs

- `results/paad/ml/random4/{input_name}/`
- `results/paad/ml/groupcv4_drug/{input_name}/`
- `results/paad/dl/random4/{input_name}/`
- `results/paad/dl/groupcv4_drug/{input_name}/`
- `results/paad/ensemble/random4/{input_name}/`
- `results/paad/ensemble/groupcv4_drug/{input_name}/`
- `external_validation/paad/groupcv4_drug/top50_external_validation.csv`
- `phase5_final_results/paad/groupcv4_drug/final_comprehensive_candidates.csv`
"""
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(text, encoding="utf-8")
    write_json(paths.reports_dir / "paad_full_pipeline_run_summary_20260421.json", run_summary)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the full PAAD drug repurposing pipeline")
    parser.add_argument("--config", default="config/paad_pipeline_config.json")
    parser.add_argument("--staging-dir", default="data/source_staging")
    parser.add_argument("--input-set", default=INPUT_NAME)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--skip-training", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = pipeline_paths(cfg)
    staging = Path(args.staging_dir)
    if not staging.is_absolute():
        staging = paths.root / staging

    run_summary: dict[str, Any] = {
        "status": "started",
        "config": str(Path(args.config).resolve()),
        "staging_dir": str(staging),
        "input_set": args.input_set,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    if not args.skip_build:
        log("[paad] build model-ready sources")
        run_summary["source_qc"] = build_model_ready(cfg, staging)
        log("[paad] build feature matrices")
        run_summary["feature_qc"] = write_feature_matrices(cfg)

    if not args.skip_training:
        cv_results = {}
        for cv in ["random4", "groupcv4_drug"]:
            log(f"[paad] train ML {cv}")
            ml_summary, _ = train_ml_cv(cfg, args.input_set, cv, force=args.force)
            log(f"[paad] train DL {cv}")
            dl_summary, _ = train_dl_cv(cfg, args.input_set, cv, force=args.force)
            log(f"[paad] build ensemble {cv}")
            ensemble = build_cv_ensemble(cfg, args.input_set, cv)
            cv_results[cv] = {
                "ml_best": ml_summary.iloc[0].to_dict() if not ml_summary.empty else {},
                "dl_best": dl_summary.iloc[0].to_dict() if not dl_summary.empty else {},
                "ensemble_best": ensemble.get("ensemble_metrics", [{}])[0],
            }
        run_summary["cv_results"] = cv_results

    if args.skip_training:
        run_summary["status"] = "completed_build_only"
        run_summary["completed_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    else:
        log("[paad] external validation and final ranking")
        run_summary["validation"] = external_validate_and_rank(cfg, staging, "groupcv4_drug", args.input_set)["qc"]
        report = write_final_report(cfg, args.input_set, run_summary)
        run_summary["status"] = "completed"
        run_summary["completed_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        run_summary["report"] = str(report)
    write_json(paths.reports_dir / "paad_full_pipeline_run_summary_20260421.json", run_summary)
    print(json.dumps(run_summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
