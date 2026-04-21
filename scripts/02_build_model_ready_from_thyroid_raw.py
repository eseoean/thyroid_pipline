#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import load_config, pipeline_paths, write_json, write_table

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except Exception:  # pragma: no cover - optional dependency
    Chem = None
    Descriptors = None
    rdMolDescriptors = None


TOKEN_SPLIT = re.compile(r"[;,|/]+")
NON_NAME = re.compile(r"[^a-z0-9]+")
GENE_WITH_ID = re.compile(r"^(.+?)\s+\(\d+\)$")
BRD_RE = re.compile(r"(BRD-[A-Z0-9]+)")

GENE_ALIASES = {
    "MEK1": ["MAP2K1"],
    "MEK2": ["MAP2K2"],
    "MTOR": ["MTOR"],
    "MTORC1": ["MTOR"],
    "MTORC2": ["MTOR"],
    "VEGFR": ["KDR", "FLT1", "FLT4"],
    "VEGFR1": ["FLT1"],
    "VEGFR2": ["KDR"],
    "VEGFR3": ["FLT4"],
    "FGFR": ["FGFR1", "FGFR2", "FGFR3", "FGFR4"],
    "NTRK": ["NTRK1", "NTRK2", "NTRK3"],
    "ERK1": ["MAPK3"],
    "ERK2": ["MAPK1"],
    "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
    "PARP": ["PARP1", "PARP2"],
    "CDK4/6": ["CDK4", "CDK6"],
}

NON_GENE_TERMS = {
    "OTHER",
    "UNKNOWN",
    "UNCLASSIFIED",
    "MITOSIS",
    "APOPTOSIS",
    "METABOLISM",
    "CHROMATIN",
    "CYTOSKELETON",
    "MICROTUBULE",
    "MICROTUBULE DESTABILISER",
    "MICROTUBULE STABILISER",
    "DNA REPLICATION",
    "DNA CROSSLINKER",
    "GENOME INTEGRITY",
    "PROTEIN STABILITY AND DEGRADATION",
    "EGFR SIGNALING",
    "ERK MAPK SIGNALING",
    "PI3K MTOR SIGNALING",
}


def norm_name(value: Any) -> str:
    return NON_NAME.sub("", "" if pd.isna(value) else str(value).lower())


def clean_gene_symbol(value: Any) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    match = GENE_WITH_ID.match(text)
    if match:
        text = match.group(1)
    return text.upper().strip()


def extract_brd(value: Any) -> str:
    match = BRD_RE.search("" if pd.isna(value) else str(value).upper())
    return match.group(1) if match else ""


def dedupe_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.loc[:, ~df.columns.duplicated()].copy()


def parse_targets(*values: Any) -> str:
    genes: set[str] = set()
    for value in values:
        if pd.isna(value):
            continue
        text = str(value).replace(" and ", ";").replace("+", ";")
        for raw in TOKEN_SPLIT.split(text):
            token = raw.strip()
            if not token:
                continue
            token = clean_gene_symbol(token)
            token = re.sub(r"\s+", " ", token)
            if token in GENE_ALIASES:
                genes.update(GENE_ALIASES[token])
                continue
            if token in NON_GENE_TERMS:
                continue
            if " " in token:
                continue
            if not re.match(r"^[A-Z][A-Z0-9.-]{1,14}$", token):
                continue
            genes.add(token)
    return ";".join(sorted(genes))


def canonicalize_smiles(smiles: Any) -> tuple[str, bool]:
    if pd.isna(smiles) or str(smiles).strip() == "":
        return "", False
    text = str(smiles).strip()
    if Chem is None:
        return text, True
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return "", False
    return Chem.MolToSmiles(mol, canonical=True), True


def morgan_bits(smiles: str, n_bits: int) -> list[int]:
    if Chem is None or rdMolDescriptors is None or not smiles:
        return [0] * n_bits
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [0] * n_bits
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=np.int8)
    for bit in fp.GetOnBits():
        arr[bit] = 1
    return arr.astype(int).tolist()


def smiles_descriptors(smiles: str) -> dict[str, float]:
    if Chem is None or Descriptors is None or not smiles:
        return {
            "drug_desc_hba": np.nan,
            "drug_desc_hbd": np.nan,
            "drug_desc_heavy_atoms": np.nan,
            "drug_desc_ring_count": np.nan,
            "drug_desc_rot_bonds": np.nan,
            "drug_desc_mol_wt": np.nan,
            "drug_desc_logp": np.nan,
        }
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "drug_desc_hba": np.nan,
            "drug_desc_hbd": np.nan,
            "drug_desc_heavy_atoms": np.nan,
            "drug_desc_ring_count": np.nan,
            "drug_desc_rot_bonds": np.nan,
            "drug_desc_mol_wt": np.nan,
            "drug_desc_logp": np.nan,
        }
    return {
        "drug_desc_hba": float(Descriptors.NumHAcceptors(mol)),
        "drug_desc_hbd": float(Descriptors.NumHDonors(mol)),
        "drug_desc_heavy_atoms": float(mol.GetNumHeavyAtoms()),
        "drug_desc_ring_count": float(Descriptors.RingCount(mol)),
        "drug_desc_rot_bonds": float(Descriptors.NumRotatableBonds(mol)),
        "drug_desc_mol_wt": float(Descriptors.MolWt(mol)),
        "drug_desc_logp": float(Descriptors.MolLogP(mol)),
    }


def top_variance_columns(df: pd.DataFrame, id_col: str, limit: int, priority_genes: set[str]) -> list[str]:
    numeric_cols = [c for c in df.columns if c != id_col and pd.api.types.is_numeric_dtype(df[c])]
    if not numeric_cols:
        return []
    variances = df[numeric_cols].var(axis=0, numeric_only=True).fillna(0).sort_values(ascending=False)
    selected = list(variances.head(limit).index)
    selected_set = set(selected)
    for col in numeric_cols:
        gene = clean_gene_symbol(col)
        if gene in priority_genes and col not in selected_set:
            selected.append(col)
            selected_set.add(col)
    return selected[: limit + len(priority_genes)]


def build_priority_gene_set(cfg: dict[str, Any]) -> set[str]:
    priority = {g.upper() for g in cfg.get("thyroid_biology_terms", []) if re.match(r"^[A-Za-z0-9.-]+$", g)}
    priority.update(
        {
            "BRAF",
            "RET",
            "NTRK1",
            "NTRK2",
            "NTRK3",
            "KRAS",
            "NRAS",
            "HRAS",
            "MAP2K1",
            "MAP2K2",
            "MAPK1",
            "MAPK3",
            "MTOR",
            "PIK3CA",
            "PIK3CB",
            "PIK3CD",
            "PIK3CG",
            "AKT1",
            "AKT2",
            "CDK4",
            "CDK6",
            "TP53",
            "TERT",
            "PTEN",
            "KDR",
            "FLT1",
            "FLT4",
        }
    )
    return priority


def gene_symbol_from_omics_col(col: str) -> str:
    text = str(col).strip()
    if text.lower() in {"modelid", "depmap_id", "profileid"}:
        return text
    text = re.sub(r"\s+\(\d+\)$", "", text)
    text = text.split(" ")[0] if text else text
    return clean_gene_symbol(text)


def select_omics_columns(df: pd.DataFrame, id_col: str, limit: int, priority_genes: set[str]) -> list[str]:
    numeric_cols = [c for c in df.columns if c != id_col and pd.api.types.is_numeric_dtype(df[c])]
    if not numeric_cols:
        return []
    variances = df[numeric_cols].var(axis=0, numeric_only=True).fillna(0).sort_values(ascending=False)
    selected = list(variances.head(limit).index)
    selected_set = set(selected)
    for col in numeric_cols:
        symbol = gene_symbol_from_omics_col(col)
        if symbol in priority_genes and col not in selected_set:
            selected.append(col)
            selected_set.add(col)
    return selected[: limit + len(priority_genes)]


def read_depmap_wide_omics(path: Path, model_ids: set[str], prefix: str, limit: int, priority_genes: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    if not path.exists():
        return pd.DataFrame({"ModelID": sorted(model_ids)}), {
            f"{prefix}_source_found": False,
            f"{prefix}_feature_rows": 0,
            f"{prefix}_selected_features": 0,
        }
    df = pd.read_csv(path)
    id_col = "ModelID" if "ModelID" in df.columns else df.columns[0]
    df[id_col] = df[id_col].astype(str)
    df = df.rename(columns={id_col: "ModelID"})
    selected = select_omics_columns(df, "ModelID", limit, priority_genes)
    small = df.loc[df["ModelID"].isin(model_ids), ["ModelID"] + selected].copy()
    rename: dict[str, str] = {}
    seen: set[str] = set()
    for col in selected:
        symbol = gene_symbol_from_omics_col(col)
        if not symbol or symbol in seen:
            continue
        rename[col] = f"sample__depmap_{prefix}__{symbol}"
        seen.add(symbol)
    small = small.rename(columns=rename)
    keep_cols = ["ModelID"] + list(rename.values())
    small = dedupe_columns(small[keep_cols])
    value_cols = [c for c in small.columns if c != "ModelID"]
    small[value_cols] = small[value_cols].apply(pd.to_numeric, errors="coerce")
    qc = {
        f"{prefix}_source_found": True,
        f"{prefix}_source_rows": int(len(df)),
        f"{prefix}_source_features": int(len(df.columns) - 1),
        f"{prefix}_feature_rows": int(small["ModelID"].nunique()),
        f"{prefix}_selected_features": int(len(value_cols)),
        f"{prefix}_matched_requested_models": int(small["ModelID"].isin(model_ids).sum()),
    }
    return small, qc


def read_depmap_mutation_features(path: Path, model_ids: set[str], priority_genes: set[str], top_limit: int) -> tuple[pd.DataFrame, dict[str, Any]]:
    base = pd.DataFrame({"ModelID": sorted(model_ids)})
    if not path.exists():
        return base, {"mutation_source_found": False, "mutation_feature_rows": 0, "mutation_selected_features": 0}

    header = pd.read_csv(path, nrows=0).columns.tolist()
    model_col = "ModelID" if "ModelID" in header else "DepMap_ID" if "DepMap_ID" in header else ""
    gene_col = "HugoSymbol" if "HugoSymbol" in header else "Hugo_Symbol" if "Hugo_Symbol" in header else "Gene" if "Gene" in header else ""
    if not model_col or not gene_col:
        return base, {
            "mutation_source_found": True,
            "mutation_feature_rows": 0,
            "mutation_selected_features": 0,
            "mutation_error": f"Required columns not found in {path.name}",
        }

    counts: dict[str, int] = defaultdict(int)
    model_gene: dict[str, set[str]] = {m: set() for m in model_ids}
    for chunk in pd.read_csv(path, usecols=[model_col, gene_col], chunksize=250_000):
        chunk[model_col] = chunk[model_col].astype(str)
        chunk = chunk.loc[chunk[model_col].isin(model_ids)].copy()
        if chunk.empty:
            continue
        chunk[gene_col] = chunk[gene_col].map(clean_gene_symbol)
        chunk = chunk.loc[chunk[gene_col].ne("")]
        for gene, count in chunk[gene_col].value_counts().items():
            counts[gene] += int(count)
        for model_id, genes in chunk.groupby(model_col)[gene_col]:
            model_gene.setdefault(str(model_id), set()).update(set(genes))

    top_genes = [g for g, _ in sorted(counts.items(), key=lambda x: (-x[1], x[0]))[:top_limit]]
    selected = sorted(set(top_genes).union(priority_genes))
    selected = [g for g in selected if g]
    feature_data = {
        f"sample__depmap_mut__{gene}": [int(gene in model_gene.get(mid, set())) for mid in base["ModelID"]]
        for gene in selected
    }
    feature_data["sample_has_depmap_mutation"] = [int(bool(model_gene.get(mid, set()))) for mid in base["ModelID"]]
    base = pd.concat([base, pd.DataFrame(feature_data)], axis=1)
    qc = {
        "mutation_source_found": True,
        "mutation_feature_rows": int(sum(bool(model_gene.get(mid, set())) for mid in model_ids)),
        "mutation_selected_features": int(len(selected)),
        "mutation_top_genes_from_thyroid_models": top_genes[:50],
    }
    return base, qc


def read_required(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def build_gdsc_response_and_base(staging: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gdsc_dir = staging / "gdsc"
    gdsc = pd.read_parquet(gdsc_dir / "gdsc2_annotation_normalized_20260406.parquet")
    thca = gdsc.loc[gdsc["TCGA_DESC"].astype(str).str.upper().eq("THCA")].copy()
    thca["sample_id"] = thca["SANGER_MODEL_ID"].astype(str)
    thca["canonical_drug_id"] = thca["DRUG_ID"].astype(str)
    thca["cell_line_name"] = thca["CELL_LINE_NAME"].astype(str)
    thca["drug_name"] = thca["DRUG_NAME"].astype(str)
    thca["disease_label"] = "thyroid cancer"
    response = thca[
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
    cell = cell.loc[cell["TCGA_DESC"].astype(str).str.upper().eq("THCA")].copy()
    drug_ann = pd.read_parquet(gdsc_dir / "gdsc2_drug_annotation_table_20260406.parquet").copy()
    drug_ann["canonical_drug_id"] = drug_ann["DRUG_ID"].astype(str)
    return response, cell, drug_ann


def build_sample_features(staging: Path, cell: pd.DataFrame, cfg: dict[str, Any], feature_limit: int) -> tuple[pd.DataFrame, dict[str, Any]]:
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
        "Age",
    ]
    model = model[[c for c in keep_model_cols if c in model.columns]].copy()

    cell = cell.merge(model, left_on="SANGER_MODEL_ID", right_on="SangerModelID", how="left")
    sample = pd.DataFrame(
        {
            "sample_id": cell["SANGER_MODEL_ID"].astype(str),
            "cell_line_name": cell["CELL_LINE_NAME"].astype(str),
            "disease_label": cell["OncotreeSubtype"].fillna(cell["TCGA_DESC"]).astype(str) + " thyroid cancer",
            "sample_cosmic_id": pd.to_numeric(cell["COSMIC_ID"], errors="coerce"),
            "sample_has_depmap_model": cell["ModelID"].notna().astype(int),
            "sample__depmap_age": pd.to_numeric(cell.get("Age"), errors="coerce") if "Age" in cell.columns else np.nan,
        }
    )
    subtype_text = cell.get("OncotreeSubtype", pd.Series([""] * len(cell))).fillna("").astype(str)
    sample["sample__is_papillary"] = subtype_text.str.contains("Papillary", case=False).astype(int)
    sample["sample__is_follicular"] = subtype_text.str.contains("Follicular", case=False).astype(int)
    sample["sample__is_anaplastic"] = subtype_text.str.contains("Anaplastic", case=False).astype(int)
    sample["sample__is_medullary"] = subtype_text.str.contains("Medullary", case=False).astype(int)
    sample["sample__is_poorly_differentiated"] = subtype_text.str.contains("Poorly", case=False).astype(int)

    crispr = pd.read_parquet(depmap_dir / "depmap_crispr_gene_dependency_basic_clean_20260406.parquet")
    priority = build_priority_gene_set(cfg)
    selected = top_variance_columns(crispr, "ModelID", feature_limit, priority)
    crispr_small = crispr[["ModelID"] + selected].copy()
    rename = {}
    for col in selected:
        gene = clean_gene_symbol(col)
        rename[col] = f"sample__depmap_dependency__{gene}"
    crispr_small = crispr_small.rename(columns=rename)
    sample = sample.merge(cell[["SANGER_MODEL_ID", "ModelID"]], left_on="sample_id", right_on="SANGER_MODEL_ID", how="left")
    sample = sample.merge(crispr_small, on="ModelID", how="left")
    sample["sample_has_crispr"] = sample["ModelID"].isin(set(crispr["ModelID"])).astype(int)
    model_ids = set(sample["ModelID"].dropna().astype(str))

    enhancement_cfg = cfg.get("sample_feature_enhancement", {})
    expr_limit = int(enhancement_cfg.get("expression_feature_limit", 512))
    cnv_limit = int(enhancement_cfg.get("cnv_feature_limit", 512))
    mut_limit = int(enhancement_cfg.get("mutation_top_gene_limit", 128))
    expr_small, expr_qc = read_depmap_wide_omics(
        depmap_dir / "OmicsExpressionProteinCodingGenesTPMLogp1_24Q2.csv",
        model_ids,
        "expr",
        expr_limit,
        priority,
    )
    cnv_small, cnv_qc = read_depmap_wide_omics(
        depmap_dir / "OmicsCNGene_24Q2.csv",
        model_ids,
        "cnv",
        cnv_limit,
        priority,
    )
    mut_small, mut_qc = read_depmap_mutation_features(
        depmap_dir / "OmicsSomaticMutations_24Q2.csv",
        model_ids,
        priority,
        mut_limit,
    )
    sample = sample.merge(expr_small, on="ModelID", how="left")
    sample = sample.merge(cnv_small, on="ModelID", how="left")
    sample = sample.merge(mut_small, on="ModelID", how="left")
    expr_feature_cols = [c for c in sample.columns if c.startswith("sample__depmap_expr__")]
    cnv_feature_cols = [c for c in sample.columns if c.startswith("sample__depmap_cnv__")]
    mut_feature_cols = [c for c in sample.columns if c.startswith("sample__depmap_mut__")]
    sample["sample_has_depmap_expression"] = sample[expr_feature_cols].notna().any(axis=1).astype(int) if expr_feature_cols else 0
    sample["sample_has_depmap_cnv"] = sample[cnv_feature_cols].notna().any(axis=1).astype(int) if cnv_feature_cols else 0
    if "sample_has_depmap_mutation" not in sample.columns:
        sample["sample_has_depmap_mutation"] = 0
    sample["sample_has_non_crispr_omics_fallback"] = (
        sample[["sample_has_depmap_expression", "sample_has_depmap_cnv", "sample_has_depmap_mutation"]].sum(axis=1).gt(0).astype(int)
    )
    missing_crispr = sample["sample_has_crispr"].eq(0)
    fallback_covered = int((missing_crispr & sample["sample_has_non_crispr_omics_fallback"].eq(1)).sum())
    sample = sample.drop(columns=[c for c in ["SANGER_MODEL_ID", "ModelID"] if c in sample.columns])

    qc = {
        "thca_cell_lines": int(cell["SANGER_MODEL_ID"].nunique()),
        "depmap_model_matches": int(cell["ModelID"].notna().sum()),
        "crispr_feature_rows": int(sample["sample_has_crispr"].sum()),
        "selected_crispr_features": int(len(selected)),
        **expr_qc,
        **cnv_qc,
        **mut_qc,
        "missing_crispr_cell_lines": int(missing_crispr.sum()),
        "missing_crispr_with_non_crispr_omics_fallback": fallback_covered,
        "selected_non_crispr_omics_features": int(len(expr_feature_cols) + len(cnv_feature_cols) + len(mut_feature_cols)),
    }
    return sample, qc


def load_drug_reference_maps(staging: Path) -> dict[str, Any]:
    refs: dict[str, Any] = {}
    gdsc_dir = staging / "gdsc"
    catalog_path = gdsc_dir / "drug_features_catalog_brca_reference_20260420.parquet"
    if catalog_path.exists():
        cat = pd.read_parquet(catalog_path)
        cat["canonical_drug_id"] = cat["DRUG_ID"].astype(str)
        refs["catalog_by_id"] = cat.drop_duplicates("canonical_drug_id").set_index("canonical_drug_id")

    drugbank = pd.read_parquet(
        staging / "drugbank" / "drugbank_drug_master_basic_20260406.parquet",
        columns=["drugbank_id", "name", "smiles", "indication", "mechanism_of_action"],
    )
    drugbank["_norm"] = drugbank["name"].map(norm_name)
    refs["drugbank_by_norm"] = drugbank.drop_duplicates("_norm").set_index("_norm")

    target = pd.read_parquet(staging / "drugbank" / "drugbank_target_table_basic_20260406.parquet")
    target_genes = (
        target.dropna(subset=["gene_name"])
        .groupby("drugbank_id")["gene_name"]
        .apply(lambda s: ";".join(sorted(set(clean_gene_symbol(x) for x in s if clean_gene_symbol(x)))))
        .to_dict()
    )
    refs["drugbank_targets"] = target_genes

    chembl = pd.read_parquet(
        staging / "chembl" / "chembl_compound_master_basic_20260406.parquet",
        columns=["chembl_id", "pref_name", "canonical_smiles", "max_phase", "first_approval"],
    )
    chembl = chembl.loc[chembl["pref_name"].notna()].copy()
    chembl["_norm"] = chembl["pref_name"].map(norm_name)
    refs["chembl_by_norm"] = chembl.drop_duplicates("_norm").set_index("_norm")
    return refs


def build_lincs_bridge_candidates(drugs: pd.DataFrame, lincs: pd.DataFrame, staging: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    pert_path = staging / "lincs" / "lincs_pert_info_basic_20260406.parquet"
    if not pert_path.exists():
        return pd.DataFrame(), {"lincs_bridge_metadata_found": False, "lincs_bridge_candidates": 0}

    pert = pd.read_parquet(pert_path)
    pert = pert.loc[pert.get("pert_type", "").astype(str).eq("trt_cp")].copy()
    if pert.empty:
        return pd.DataFrame(), {"lincs_bridge_metadata_found": True, "lincs_bridge_candidates": 0}
    pert["_norm"] = pert["pert_iname"].map(norm_name)
    pert["_brd_id"] = pert["pert_id"].map(extract_brd)
    pert["_canonical_smiles_norm"] = pert["canonical_smiles"].map(lambda x: canonicalize_smiles(x)[0])
    pert["_touchstone_rank"] = pd.to_numeric(pert.get("is_touchstone", 0), errors="coerce").fillna(0)
    pert = pert.sort_values(["_touchstone_rank", "pert_id"], ascending=[False, True])

    direct_ids = set(lincs["canonical_drug_id"].astype(str)) if "canonical_drug_id" in lincs.columns else set()
    missing = drugs.loc[~drugs["canonical_drug_id"].astype(str).isin(direct_ids)].copy()
    rows: list[dict[str, Any]] = []
    used_ids: set[str] = set()
    for row in missing.to_dict(orient="records"):
        drug_id = str(row["canonical_drug_id"])
        norm = norm_name(row.get("drug_name", ""))
        can = canonicalize_smiles(row.get("canonical_smiles", ""))[0]
        hit = pd.DataFrame()
        match_type = ""
        if norm:
            hit = pert.loc[pert["_norm"].eq(norm)].head(1)
            match_type = "name"
        if hit.empty and can:
            hit = pert.loc[pert["_canonical_smiles_norm"].eq(can)].head(1)
            match_type = "smiles"
        if hit.empty:
            continue
        h = hit.iloc[0]
        brd_id = str(h.get("_brd_id", ""))
        if not brd_id or drug_id in used_ids:
            continue
        used_ids.add(drug_id)
        rows.append(
            {
                "canonical_drug_id": drug_id,
                "drug_name": row.get("drug_name", ""),
                "pert_id": h.get("pert_id", ""),
                "brd_id": brd_id,
                "pert_iname": h.get("pert_iname", ""),
                "match_type": match_type,
                "is_touchstone": int(h.get("is_touchstone", 0) or 0),
            }
        )
    bridge = pd.DataFrame(rows)
    return bridge, {
        "lincs_bridge_metadata_found": True,
        "lincs_bridge_candidates": int(len(bridge)),
        "lincs_bridge_match_type_counts": bridge["match_type"].value_counts().to_dict() if not bridge.empty else {},
    }


def recover_lincs_signature_from_mcf7(
    mcf7_path: Path,
    bridge: pd.DataFrame,
    selected_lincs_cols: list[str],
    rename_map: dict[str, str],
    gene_info_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    if bridge.empty:
        return pd.DataFrame(), bridge, {"lincs_mcf7_source_found": bool(mcf7_path.exists()), "lincs_recovered_drugs": 0}
    if not mcf7_path.exists():
        bridge = bridge.copy()
        bridge["recovery_status"] = "mcf7_source_missing"
        return pd.DataFrame(), bridge, {"lincs_mcf7_source_found": False, "lincs_recovered_drugs": 0}

    gene_info = pd.read_parquet(gene_info_path)
    symbol_to_entrez = {
        clean_gene_symbol(row.pr_gene_symbol): str(row.pr_gene_id)
        for row in gene_info[["pr_gene_id", "pr_gene_symbol"]].dropna().itertuples(index=False)
    }
    schema_cols = set(pq.read_schema(mcf7_path).names)
    raw_to_out: dict[str, str] = {}
    for selected in selected_lincs_cols:
        out_col = rename_map[selected]
        if selected in schema_cols:
            raw_to_out[selected] = out_col
            continue
        symbol = clean_gene_symbol(selected.replace("crispr__", ""))
        entrez = symbol_to_entrez.get(symbol, "")
        if entrez and entrez in schema_cols:
            raw_to_out[entrez] = out_col

    if not raw_to_out or "sig_id" not in schema_cols:
        bridge = bridge.copy()
        bridge["recovery_status"] = "mcf7_schema_incompatible"
        return pd.DataFrame(), bridge, {
            "lincs_mcf7_source_found": True,
            "lincs_recovered_drugs": 0,
            "lincs_mcf7_selected_raw_features": 0,
        }

    brd_to_drug = dict(zip(bridge["brd_id"].astype(str), bridge["canonical_drug_id"].astype(str)))
    wanted_brd = set(brd_to_drug)
    read_cols = ["sig_id"] + list(raw_to_out)
    pieces: list[pd.DataFrame] = []
    pf = pq.ParquetFile(mcf7_path)
    for i in range(pf.num_row_groups):
        rg = pf.read_row_group(i, columns=read_cols).to_pandas()
        rg["brd_id"] = rg["sig_id"].map(extract_brd)
        rg = rg.loc[rg["brd_id"].isin(wanted_brd)].copy()
        if rg.empty:
            continue
        rg["canonical_drug_id"] = rg["brd_id"].map(brd_to_drug)
        value_cols = list(raw_to_out)
        rg[value_cols] = rg[value_cols].apply(pd.to_numeric, errors="coerce")
        pieces.append(rg[["canonical_drug_id", "brd_id"] + value_cols])

    bridge = bridge.copy()
    if not pieces:
        bridge["recovery_status"] = "no_mcf7_signature_rows"
        return pd.DataFrame(), bridge, {
            "lincs_mcf7_source_found": True,
            "lincs_recovered_drugs": 0,
            "lincs_mcf7_selected_raw_features": int(len(raw_to_out)),
            "lincs_mcf7_row_groups_scanned": int(pf.num_row_groups),
        }

    matched = pd.concat(pieces, ignore_index=True)
    grouped = matched.groupby("canonical_drug_id", as_index=False)[list(raw_to_out)].mean()
    grouped = grouped.rename(columns=raw_to_out)
    source_counts = matched.groupby("canonical_drug_id").size().rename("mcf7_signature_rows").reset_index()
    bridge = bridge.merge(source_counts, on="canonical_drug_id", how="left")
    bridge["recovery_status"] = np.where(bridge["mcf7_signature_rows"].fillna(0).gt(0), "recovered_from_mcf7", "no_mcf7_signature_rows")
    qc = {
        "lincs_mcf7_source_found": True,
        "lincs_recovered_drugs": int(grouped["canonical_drug_id"].nunique()),
        "lincs_mcf7_signature_rows_used": int(len(matched)),
        "lincs_mcf7_selected_raw_features": int(len(raw_to_out)),
        "lincs_mcf7_row_groups_scanned": int(pf.num_row_groups),
    }
    return grouped, bridge, qc


def build_drug_features(staging: Path, drug_ann: pd.DataFrame, cfg: dict[str, Any], lincs_limit: int, morgan_bits_n: int) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    refs = load_drug_reference_maps(staging)
    known = {norm_name(x) for x in cfg.get("known_thyroid_drugs", [])}
    clinical_trials_text = ""
    ct_path = staging / "clinical_trials" / "clinicaltrials_thyroid_cancer_drug_20260420.json"
    if ct_path.exists():
        clinical_trials_text = ct_path.read_text(encoding="utf-8", errors="ignore").lower()

    records = []
    for row in drug_ann.itertuples(index=False):
        drug_id = str(row.canonical_drug_id)
        drug_name = str(row.DRUG_NAME)
        norm = norm_name(drug_name)
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
                source = "brca_gdsc_catalog"
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
        can, ok = canonicalize_smiles(smiles)
        db_targets = refs["drugbank_targets"].get(drugbank_id, "")
        target_genes = parse_targets(db_targets, row.PUTATIVE_TARGET_NORMALIZED, row.PUTATIVE_TARGET)
        indication_low = f"{indication} {clinical_trials_text if norm and norm in clinical_trials_text else ''}".lower()
        if norm in known:
            classification = "approved"
        elif "thyroid" in indication_low:
            classification = "indication_expansion"
        else:
            classification = "screened_candidate"
        bridge_bits = [
            ok,
            bool(target_genes),
            norm in refs["drugbank_by_norm"].index,
            norm in refs["chembl_by_norm"].index,
        ]
        bridge_strength = "multi_source" if sum(bridge_bits) >= 3 else "single_source" if sum(bridge_bits) >= 1 else "weak"
        desc = smiles_descriptors(can)
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
            "TCGA_DESC": cfg["project"].get("tcga_code", "THCA"),
            "drug__has_smiles": int(ok),
            "drug__has_target": int(bool(target_genes)),
            "drug__is_known_thyroid_control": int(norm in known),
            "drug__is_clinical_trial_mentioned": int(norm in clinical_trials_text) if norm else 0,
        }
        rec.update(desc)
        for i, bit in enumerate(morgan_bits(can, morgan_bits_n)):
            rec[f"drug_morgan_{i:04d}"] = bit
        records.append(rec)

    drugs = pd.DataFrame(records)

    lincs_path = staging / "lincs" / "lincs_drug_signature_normalized.parquet"
    lincs_qc = {"lincs_rows": 0, "selected_lincs_features": 0, "drug_lincs_matches": 0, "drug_lincs_recovered_matches": 0}
    if lincs_path.exists():
        lincs = pd.read_parquet(lincs_path)
        lincs["canonical_drug_id"] = lincs["canonical_drug_id"].astype(str)
        selected = top_variance_columns(lincs, "canonical_drug_id", lincs_limit, set())
        rename = {c: f"drug__lincs__{clean_gene_symbol(c.replace('crispr__', ''))}" for c in selected}
        lincs_small = lincs[["canonical_drug_id"] + selected].rename(columns=rename)
        lincs_small["drug__has_lincs_signature"] = 1
        lincs_small["drug__lincs_signature_source_direct"] = 1
        lincs_small["drug__lincs_signature_source_recovered"] = 0

        bridge, bridge_qc = build_lincs_bridge_candidates(drugs, lincs, staging)
        recovered, bridge_review, recovery_qc = recover_lincs_signature_from_mcf7(
            staging / "lincs" / "lincs_mcf7.parquet",
            bridge,
            selected,
            rename,
            staging / "lincs" / "lincs_gene_info_basic_20260406.parquet",
        )
        if not bridge_review.empty:
            write_table(bridge_review, REPO_ROOT / "reports" / "missingness" / "lincs_signature_recovery_review.csv")
        recovered_count = 0
        if not recovered.empty:
            recovered["drug__has_lincs_signature"] = 1
            recovered["drug__lincs_signature_source_direct"] = 0
            recovered["drug__lincs_signature_source_recovered"] = 1
            recovered_count = int(recovered["canonical_drug_id"].nunique())
            lincs_small = pd.concat([lincs_small, recovered], ignore_index=True, sort=False)
            lincs_small = lincs_small.drop_duplicates("canonical_drug_id", keep="first")

        drugs = drugs.merge(lincs_small, on="canonical_drug_id", how="left")
        drugs["drug__has_lincs_signature"] = drugs["drug__has_lincs_signature"].fillna(0).astype(int)
        drugs["drug__lincs_signature_source_direct"] = drugs["drug__lincs_signature_source_direct"].fillna(0).astype(int)
        drugs["drug__lincs_signature_source_recovered"] = drugs["drug__lincs_signature_source_recovered"].fillna(0).astype(int)
        lincs_qc = {
            "lincs_rows": int(len(lincs)),
            "selected_lincs_features": int(len(selected)),
            "drug_lincs_direct_matches": int(drugs["canonical_drug_id"].isin(set(lincs["canonical_drug_id"])).sum()),
            "drug_lincs_recovered_matches": recovered_count,
            "drug_lincs_matches": int(drugs["drug__has_lincs_signature"].sum()),
            **bridge_qc,
            **recovery_qc,
        }
    else:
        drugs["drug__has_lincs_signature"] = 0
        drugs["drug__lincs_signature_source_direct"] = 0
        drugs["drug__lincs_signature_source_recovered"] = 0

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
        "morgan_bits": int(morgan_bits_n),
        **lincs_qc,
    }
    return drugs, annotations, qc


def build_external_validation_files(staging: Path, output_dir: Path, target_genes: set[str]) -> dict[str, Any]:
    tcga = staging / "tcga_thca"
    gencode = staging / "gencode" / "gencode.v36.annotation.gtf.gz"
    expression = tcga / "TCGA-THCA.star_tpm.tsv.gz"
    survival = tcga / "TCGA-THCA.survival.tsv.gz"
    clinical = tcga / "TCGA-THCA.clinical.tsv.gz"
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_map: dict[str, str] = {}
    with gzip.open(gencode, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            gid = re.search(r'gene_id "([^"]+)"', attrs)
            gname = re.search(r'gene_name "([^"]+)"', attrs)
            if gid and gname:
                gene_map[gid.group(1).split(".")[0]] = gname.group(1).upper()

    rows = []
    header = None
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
    write_table(expr_df, output_dir / "thyroid_expression.csv")

    surv = pd.read_csv(survival, sep="\t", compression="gzip")
    clin = pd.read_csv(clinical, sep="\t", compression="gzip")
    clinical_out = surv.rename(columns={"sample": "PATIENT_ID", "OS.time": "OS_DAYS", "OS": "OS_STATUS"}).copy()
    clinical_out["OS_MONTHS"] = pd.to_numeric(clinical_out["OS_DAYS"], errors="coerce") / 30.4375
    if "sample" in clin.columns:
        clin_small = clin[["sample"] + [c for c in ["gender.demographic", "age_at_index.demographic", "tumor_stage.diagnoses"] if c in clin.columns]].copy()
        clin_small = clin_small.rename(columns={"sample": "PATIENT_ID"})
        clinical_out = clinical_out.merge(clin_small, on="PATIENT_ID", how="left")
    write_table(clinical_out, output_dir / "thyroid_clinical.csv")

    return {
        "target_gene_universe": int(len(target_upper)),
        "external_expression_genes_written": int(expr_df["Hugo_Symbol"].nunique()) if not expr_df.empty else 0,
        "external_expression_samples": int(max(0, expr_df.shape[1] - 1)),
        "clinical_rows": int(len(clinical_out)),
    }


def build_admet_csvs(staging: Path, admet_out: Path) -> dict[str, Any]:
    admet_out = admet_out / "tdc"
    admet_out.mkdir(parents=True, exist_ok=True)
    source = staging / "admet"
    assays = defaultdict(list)
    for path in source.glob("*/*_basic_clean_20260406.parquet"):
        assays[path.parent.name].append(path)

    written = []
    for assay, paths in sorted(assays.items()):
        frames = []
        for path in sorted(paths):
            df = pd.read_parquet(path)
            if "Drug" not in df.columns or "Y" not in df.columns:
                continue
            tmp = df[["Drug", "Y"]].rename(columns={"Drug": "smiles", "Y": "label"})
            frames.append(tmp)
        if not frames:
            continue
        out = pd.concat(frames, ignore_index=True).dropna().drop_duplicates()
        out.to_csv(admet_out / f"{assay}.csv", index=False)
        written.append({"assay": assay, "rows": int(len(out))})
    return {"admet_assays_written": len(written), "admet_rows_by_assay": written}


def apply_primary_candidate_filters(
    response: pd.DataFrame,
    drug_features: pd.DataFrame,
    drug_annotations: pd.DataFrame,
    cfg: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    filters = cfg.get("data_filters", {})
    require_smiles = bool(filters.get("require_valid_smiles_drugs", False))
    before = {
        "response_rows": int(len(response)),
        "response_cell_lines": int(response["sample_id"].nunique()),
        "response_drugs": int(response["canonical_drug_id"].nunique()),
        "drug_feature_rows": int(len(drug_features)),
    }
    qc: dict[str, Any] = {
        "require_valid_smiles_drugs": require_smiles,
        "before": before,
    }
    if not require_smiles:
        qc["after"] = before
        qc["removed"] = {"response_rows": 0, "drugs": 0}
        return response, drug_features, drug_annotations, qc

    smiles_ok = (
        drug_features["canonical_smiles"].fillna("").astype(str).str.len().gt(0)
        & drug_features.get("smiles_parse_ok", pd.Series(0, index=drug_features.index)).fillna(0).astype(int).eq(1)
    )
    keep_ids = set(drug_features.loc[smiles_ok, "canonical_drug_id"].astype(str))
    missing = drug_features.loc[~smiles_ok, ["canonical_drug_id", "drug_name", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"]].copy()

    response = response.loc[response["canonical_drug_id"].astype(str).isin(keep_ids)].copy()
    drug_features = drug_features.loc[drug_features["canonical_drug_id"].astype(str).isin(keep_ids)].copy()
    drug_annotations = drug_annotations.loc[drug_annotations["canonical_drug_id"].astype(str).isin(keep_ids)].copy()

    after = {
        "response_rows": int(len(response)),
        "response_cell_lines": int(response["sample_id"].nunique()),
        "response_drugs": int(response["canonical_drug_id"].nunique()),
        "drug_feature_rows": int(len(drug_features)),
    }
    qc.update(
        {
            "after": after,
            "removed": {
                "response_rows": int(before["response_rows"] - after["response_rows"]),
                "drugs": int(before["response_drugs"] - after["response_drugs"]),
            },
            "removed_drugs": missing.to_dict(orient="records"),
        }
    )
    return response, drug_features, drug_annotations, qc


def main() -> int:
    parser = argparse.ArgumentParser(description="Build model-ready thyroid pipeline inputs from s3 thyroid_raw source staging")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--staging-dir", default="data/source_staging")
    parser.add_argument("--sample-feature-limit", type=int, default=4096)
    parser.add_argument("--lincs-feature-limit", type=int, default=1024)
    parser.add_argument("--morgan-bits", type=int, default=512)
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = pipeline_paths(cfg)
    staging = Path(args.staging_dir)
    if not staging.is_absolute():
        staging = paths.root / staging

    response, cell, drug_ann = build_gdsc_response_and_base(staging)
    sample_features, sample_qc = build_sample_features(staging, cell, cfg, args.sample_feature_limit)
    drug_features, drug_annotations, drug_qc = build_drug_features(staging, drug_ann, cfg, args.lincs_feature_limit, args.morgan_bits)
    response, drug_features, drug_annotations, filter_qc = apply_primary_candidate_filters(response, drug_features, drug_annotations, cfg)

    target_genes = set()
    for value in drug_features["target_genes"].dropna():
        target_genes.update([g for g in str(value).split(";") if g])
    target_genes.update({"BRAF", "RET", "NTRK1", "NTRK2", "NTRK3", "KRAS", "NRAS", "HRAS", "MAP2K1", "MAP2K2", "MTOR", "PIK3CA", "CDK4", "CDK6", "KDR", "FLT1", "FLT4"})
    external_qc = build_external_validation_files(staging, paths.external_dir, target_genes)
    admet_qc = build_admet_csvs(staging, paths.admet_source_dir)

    write_table(response, paths.raw_dir / "thyroid_response_pairs.csv")
    write_table(sample_features, paths.raw_dir / "sample_features.csv")
    write_table(drug_features, paths.raw_dir / "drug_features.csv")
    write_table(drug_annotations, paths.raw_dir / "drug_annotations.csv")

    qc = {
        "step": "source_to_model_ready",
        "source_staging_dir": str(staging),
        "response_rows": int(len(response)),
        "response_cell_lines": int(response["sample_id"].nunique()),
        "response_drugs": int(response["canonical_drug_id"].nunique()),
        "sample_features_shape": [int(sample_features.shape[0]), int(sample_features.shape[1])],
        "drug_features_shape": [int(drug_features.shape[0]), int(drug_features.shape[1])],
        "drug_annotations_shape": [int(drug_annotations.shape[0]), int(drug_annotations.shape[1])],
        "sample_qc": sample_qc,
        "drug_qc": drug_qc,
        "primary_filter_qc": filter_qc,
        "external_qc": external_qc,
        "admet_qc": admet_qc,
        "outputs": {
            "response_table": str(paths.raw_dir / "thyroid_response_pairs.csv"),
            "sample_features": str(paths.raw_dir / "sample_features.csv"),
            "drug_features": str(paths.raw_dir / "drug_features.csv"),
            "drug_annotations": str(paths.raw_dir / "drug_annotations.csv"),
            "external_expression": str(paths.external_dir / "thyroid_expression.csv"),
            "external_clinical": str(paths.external_dir / "thyroid_clinical.csv"),
            "admet_dir": str(paths.admet_source_dir),
        },
    }
    write_json(paths.reports_dir / "qc_source_to_model_ready_20260420.json", qc)
    print(json.dumps(qc, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
