#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import admet_assessment, load_config, pipeline_paths, split_genes, write_json  # noqa: E402


INPUT_NAME = "numeric_strong_context_smiles_pan_lincs"
CV_NAME = "groupcv4_drug"
VARIANT = "groupcv4_drug_v2_sources"

PAAD_SEED_GENES = {
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
    "AURKA",
    "AURKB",
    "WEE1",
}

PAAD_HALLMARKS = {
    "HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_G2M_CHECKPOINT",
}

SERIOUS_SIDE_EFFECT_KEYWORDS = {
    "death",
    "fatal",
    "cardiac",
    "heart",
    "myocardial",
    "arrhythmia",
    "neutropenia",
    "leukopenia",
    "thrombocytopenia",
    "myelosuppression",
    "sepsis",
    "infection",
    "renal",
    "kidney",
    "hepatic",
    "liver",
    "pancreatitis",
    "neuropathy",
    "pulmonary",
    "embolism",
    "hemorrhage",
}


def norm_name(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def clean_symbol(value: Any) -> str:
    text = str(value).strip().strip('"').upper()
    text = text.split("|")[0]
    text = re.sub(r"[^A-Z0-9_.-]", "", text)
    return text


def finite_mean(values: list[float]) -> float:
    arr = np.asarray([v for v in values if pd.notna(v) and np.isfinite(v)], dtype=float)
    return float(arr.mean()) if arr.size else np.nan


def clip01(value: Any) -> float:
    try:
        v = float(value)
    except Exception:
        return 0.0
    if not np.isfinite(v):
        return 0.0
    return float(min(1.0, max(0.0, v)))


def safe_mannwhitney(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan
    try:
        return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except Exception:
        return np.nan


def signature_table(
    source: str,
    tumor: dict[str, np.ndarray],
    normal: dict[str, np.ndarray],
    raw_tpm: bool = False,
) -> pd.DataFrame:
    rows = []
    for gene in sorted(set(tumor) & set(normal)):
        t = np.asarray(tumor[gene], dtype=float)
        n = np.asarray(normal[gene], dtype=float)
        if raw_tpm:
            t = np.log2(np.maximum(t, 0) + 1.0)
            n = np.log2(np.maximum(n, 0) + 1.0)
        t = t[np.isfinite(t)]
        n = n[np.isfinite(n)]
        if len(t) == 0 or len(n) == 0:
            continue
        delta = float(np.nanmedian(t) - np.nanmedian(n))
        rows.append(
            {
                "source": source,
                "gene": gene,
                "tumor_n": int(len(t)),
                "normal_n": int(len(n)),
                "tumor_median": float(np.nanmedian(t)),
                "normal_median": float(np.nanmedian(n)),
                "log2fc_or_delta": delta,
                "p_value": safe_mannwhitney(t, n),
                "support_up": bool(delta > 0.25),
            }
        )
    return pd.DataFrame(rows)


def load_candidate_tables(paths: Any, input_name: str) -> pd.DataFrame:
    model_path = paths.results_dir / "ensemble" / CV_NAME / input_name / "ensemble_top50_drugs.csv"
    v1_path = paths.external_validation_dir / CV_NAME / "top50_external_validation.csv"
    if not model_path.exists():
        raise FileNotFoundError(f"Missing model top50 candidates: {model_path}")
    model = pd.read_csv(model_path)
    model["canonical_drug_id"] = model["canonical_drug_id"].astype(str)
    if v1_path.exists():
        v1 = pd.read_csv(v1_path)
        v1["canonical_drug_id"] = v1["canonical_drug_id"].astype(str)
        missing_cols = [c for c in model.columns if c not in v1.columns]
        merged = v1.merge(model[["canonical_drug_id"] + missing_cols], on="canonical_drug_id", how="left")
    else:
        merged = model.copy()
    merged = merged.sort_values("rank").reset_index(drop=True)
    return merged


def target_universe(candidates: pd.DataFrame, cfg: dict[str, Any]) -> set[str]:
    genes: set[str] = set(PAAD_SEED_GENES)
    for value in candidates["target_genes"].fillna(""):
        genes.update(split_genes(value))
    for value in cfg.get("paad_biology_terms", []):
        if isinstance(value, str):
            gene = clean_symbol(value)
            if 2 <= len(gene) <= 16:
                genes.add(gene)
    return {g for g in genes if g}


def load_paad_external_expression(path: Path, wanted: set[str]) -> dict[str, np.ndarray]:
    df = pd.read_csv(path)
    gene_col = "Hugo_Symbol" if "Hugo_Symbol" in df.columns else df.columns[0]
    df[gene_col] = df[gene_col].map(clean_symbol)
    df = df.loc[df[gene_col].isin(wanted)].copy()
    value_cols = [c for c in df.columns if c != gene_col]
    out: dict[str, np.ndarray] = {}
    for gene, sub in df.groupby(gene_col):
        vals = sub[value_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        out[gene] = np.nanmean(vals, axis=0)
    return out


def load_gtex_pancreas(path: Path, wanted: set[str]) -> dict[str, np.ndarray]:
    df = pd.read_csv(path, sep="\t", compression="gzip", skiprows=2)
    symbol_col = "Description" if "Description" in df.columns else df.columns[2]
    df[symbol_col] = df[symbol_col].map(clean_symbol)
    df = df.loc[df[symbol_col].isin(wanted)].copy()
    value_cols = [c for c in df.columns if c not in {"id", "Name", "Description"}]
    out: dict[str, np.ndarray] = {}
    for gene, sub in df.groupby(symbol_col):
        vals = sub[value_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        out[gene] = np.nanmean(vals, axis=0)
    return out


def parse_geo_series_matrix(path: Path, wanted: set[str] | None = None) -> tuple[dict[str, list[str]], pd.DataFrame]:
    metadata: dict[str, list[str]] = {}
    rows: list[list[str]] = []
    header: list[str] | None = None
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parsed = next(csv.reader([line], delimiter="\t"))
            if not in_table:
                if parsed and parsed[0].startswith("!Sample_"):
                    metadata[parsed[0]] = parsed[1:]
                continue
            if header is None:
                header = [x.strip('"') for x in parsed]
                continue
            gene = clean_symbol(parsed[0])
            if wanted is not None and gene not in wanted:
                continue
            rows.append([gene] + parsed[1:])
    if header is None:
        return metadata, pd.DataFrame()
    df = pd.DataFrame(rows, columns=["gene"] + header[1:])
    if not df.empty:
        value_cols = [c for c in df.columns if c != "gene"]
        df[value_cols] = df[value_cols].apply(pd.to_numeric, errors="coerce")
        df = df.groupby("gene", as_index=False)[value_cols].mean()
    return metadata, df


def load_gpl6244_mapping(path: Path, wanted: set[str]) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith(("^", "!", "#"))), delimiter="\t")
        for row in reader:
            probe = str(row.get("ID", "")).strip()
            symbols = []
            for part in re.split(r"\s*///\s*|\s*//\s*|;\s*", str(row.get("Gene symbol", ""))):
                sym = clean_symbol(part)
                if sym in wanted:
                    symbols.append(sym)
            if probe and symbols:
                mapping[probe] = sorted(set(symbols))
    return mapping


def parse_gse62452(path: Path, annot_path: Path, wanted: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    probe_to_genes = load_gpl6244_mapping(annot_path, wanted)
    metadata: dict[str, list[str]] = {}
    tissue_labels: list[str] = []
    header: list[str] | None = None
    gene_rows: dict[str, list[np.ndarray]] = defaultdict(list)
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parsed = next(csv.reader([line], delimiter="\t"))
            if not in_table:
                if parsed and parsed[0].startswith("!Sample_"):
                    metadata[parsed[0]] = parsed[1:]
                    if parsed[0] == "!Sample_characteristics_ch1" and any("tissue:" in str(x).lower() for x in parsed[1:]):
                        tissue_labels = parsed[1:]
                continue
            if header is None:
                header = [x.strip('"') for x in parsed]
                continue
            probe = str(parsed[0]).strip('"')
            genes = probe_to_genes.get(probe)
            if not genes:
                continue
            values = pd.to_numeric(pd.Series(parsed[1:]), errors="coerce").to_numpy(dtype=float)
            for gene in genes:
                gene_rows[gene].append(values)
    if header is None:
        return pd.DataFrame(), {"platform_genes_matched": 0}
    sample_ids = header[1:]
    supp = metadata.get("!Sample_supplementary_file", [""] * len(sample_ids))
    tumor_idx = {i for i, text in enumerate(supp) if re.search(r"_T\.", str(text), re.IGNORECASE)}
    normal_idx = {i for i, text in enumerate(supp) if re.search(r"_N\.", str(text), re.IGNORECASE)}
    for i, text in enumerate(tissue_labels[: len(sample_ids)]):
        lower = str(text).lower()
        if "non-tumor" in lower or "normal" in lower:
            normal_idx.add(i)
        elif "tumor" in lower:
            tumor_idx.add(i)
    tumor_idx -= normal_idx
    tumor_cols = [sample_ids[i] for i in sorted(tumor_idx)]
    normal_cols = [sample_ids[i] for i in sorted(normal_idx)]
    rows = []
    for gene, arrays in gene_rows.items():
        vals = np.nanmean(np.vstack(arrays), axis=0)
        rows.append([gene] + vals.tolist())
    df = pd.DataFrame(rows, columns=["gene"] + sample_ids)
    df = df.groupby("gene", as_index=False)[sample_ids].mean() if not df.empty else df
    sig = signature_table(
        "GSE62452_tumor_vs_adjacent",
        {r.gene: pd.to_numeric(pd.Series([getattr(r, c) for c in tumor_cols]), errors="coerce").to_numpy(dtype=float) for r in df.itertuples(index=False)},
        {r.gene: pd.to_numeric(pd.Series([getattr(r, c) for c in normal_cols]), errors="coerce").to_numpy(dtype=float) for r in df.itertuples(index=False)},
        raw_tpm=False,
    )
    qc = {
        "platform_mapped_target_genes": int(df["gene"].nunique()) if not df.empty else 0,
        "tumor_samples": int(len(tumor_cols)),
        "normal_samples": int(len(normal_cols)),
    }
    return sig, qc


def parse_gse71729(path: Path, wanted: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    metadata, df = parse_geo_series_matrix(path, wanted)
    source = metadata.get("!Sample_source_name_ch2") or metadata.get("!Sample_title") or []
    sample_ids = [c for c in df.columns if c != "gene"]
    tumor_cols: list[str] = []
    normal_cols: list[str] = []
    cellline_cols: list[str] = []
    for idx, sample in enumerate(sample_ids):
        text = str(source[idx] if idx < len(source) else sample).lower()
        if "cellline" in text or "cell line" in text:
            cellline_cols.append(sample)
        elif "normal" in text:
            normal_cols.append(sample)
        elif "metastasis" in text or "primary" in text or "tumor" in text:
            tumor_cols.append(sample)
    sig = signature_table(
        "GSE71729_tumor_met_vs_normal",
        {r.gene: pd.to_numeric(pd.Series([getattr(r, c) for c in tumor_cols]), errors="coerce").to_numpy(dtype=float) for r in df.itertuples(index=False)},
        {r.gene: pd.to_numeric(pd.Series([getattr(r, c) for c in normal_cols]), errors="coerce").to_numpy(dtype=float) for r in df.itertuples(index=False)},
        raw_tpm=False,
    )
    qc = {
        "gene_rows_matched": int(df["gene"].nunique()) if not df.empty else 0,
        "tumor_or_metastasis_samples": int(len(tumor_cols)),
        "normal_samples": int(len(normal_cols)),
        "cell_line_samples_excluded": int(len(cellline_cols)),
    }
    return sig, qc


def parse_cptac_proteomics(stage: Path, wanted: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    tumor_path = stage / "cptac_paad" / "proteomics_gene_level_MD_abundance_tumor.cct"
    normal_path = stage / "cptac_paad" / "proteomics_gene_level_MD_abundance_normal.cct"
    tumor = pd.read_csv(tumor_path, sep="\t", index_col=0, na_values=["NA", ""])
    normal = pd.read_csv(normal_path, sep="\t", index_col=0, na_values=["NA", ""])
    tumor.index = [clean_symbol(x) for x in tumor.index]
    normal.index = [clean_symbol(x) for x in normal.index]
    tumor = tumor.loc[tumor.index.isin(wanted)].groupby(level=0).mean()
    normal = normal.loc[normal.index.isin(wanted)].groupby(level=0).mean()
    sig = signature_table(
        "CPTAC_proteomics_tumor_vs_normal",
        {g: pd.to_numeric(tumor.loc[g], errors="coerce").to_numpy(dtype=float) for g in tumor.index},
        {g: pd.to_numeric(normal.loc[g], errors="coerce").to_numpy(dtype=float) for g in normal.index},
        raw_tpm=False,
    )
    qc = {
        "matched_genes": int(sig["gene"].nunique()) if not sig.empty else 0,
        "tumor_samples": int(tumor.shape[1]),
        "normal_samples": int(normal.shape[1]),
    }
    return sig, qc


def parse_rppa(stage: Path, wanted: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    path = stage / "cbioportal" / "paad_tcga_pan_can_atlas_2018" / "data_rppa_zscores.txt"
    rppa = pd.read_csv(path, sep="\t")
    first = rppa.columns[0]
    rppa["gene"] = rppa[first].map(clean_symbol)
    rppa = rppa.loc[rppa["gene"].isin(wanted)].copy()
    value_cols = [c for c in rppa.columns if c not in {first, "gene"}]
    rows = []
    for gene, sub in rppa.groupby("gene"):
        values = sub[value_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        flat = values[np.isfinite(values)]
        rows.append(
            {
                "gene": gene,
                "rppa_antibodies": int(len(sub)),
                "rppa_mean_z": float(np.nanmean(flat)) if flat.size else np.nan,
                "rppa_positive_fraction": float((flat > 0).mean()) if flat.size else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    qc = {"matched_genes": int(out["gene"].nunique()) if not out.empty else 0, "rows": int(len(out))}
    return out, qc


def parse_mutations(stage: Path, wanted: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    path = stage / "cbioportal" / "paad_tcga_pan_can_atlas_2018" / "data_mutations.txt"
    muts = pd.read_csv(path, sep="\t", usecols=["Hugo_Symbol", "Tumor_Sample_Barcode"], low_memory=False)
    muts["gene"] = muts["Hugo_Symbol"].map(clean_symbol)
    muts = muts.loc[muts["gene"].isin(wanted)].copy()
    out = (
        muts.groupby("gene")
        .agg(
            cbio_mutation_count=("Tumor_Sample_Barcode", "size"),
            cbio_mutated_samples=("Tumor_Sample_Barcode", lambda x: int(pd.Series(x).nunique())),
        )
        .reset_index()
    )
    qc = {"matched_genes": int(out["gene"].nunique()) if not out.empty else 0, "mutation_rows": int(len(muts))}
    return out, qc


def parse_string(stage: Path, wanted: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    info_path = stage / "string" / "9606.protein.info.v12.0.txt.gz"
    links_path = stage / "string" / "9606.protein.links.v12.0.txt.gz"
    wanted_all = set(wanted) | PAAD_SEED_GENES
    protein_to_gene: dict[str, str] = {}
    gene_to_proteins: dict[str, set[str]] = defaultdict(set)
    with gzip.open(info_path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = clean_symbol(row.get("preferred_name", ""))
            if gene in wanted_all:
                pid = str(row.get("#string_protein_id", ""))
                protein_to_gene[pid] = gene
                gene_to_proteins[gene].add(pid)
    interesting = set(protein_to_gene)
    seed_genes = PAAD_SEED_GENES & set(gene_to_proteins)
    best_by_gene = {g: 1.0 for g in wanted if g in PAAD_SEED_GENES}
    seed_hits: dict[str, set[str]] = defaultdict(set)
    with gzip.open(links_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        for line in handle:
            p1, p2, score_raw = line.rstrip("\n").split()[:3]
            if p1 not in interesting or p2 not in interesting:
                continue
            g1 = protein_to_gene[p1]
            g2 = protein_to_gene[p2]
            if g1 == g2:
                continue
            score = float(score_raw) / 1000.0
            if g1 in wanted and g2 in seed_genes:
                best_by_gene[g1] = max(best_by_gene.get(g1, 0.0), score)
                seed_hits[g1].add(g2)
            if g2 in wanted and g1 in seed_genes:
                best_by_gene[g2] = max(best_by_gene.get(g2, 0.0), score)
                seed_hits[g2].add(g1)
    rows = [
        {
            "gene": gene,
            "string_paad_seed_max_score": float(best_by_gene.get(gene, np.nan)),
            "string_seed_hit_genes": ";".join(sorted(seed_hits.get(gene, set()))),
        }
        for gene in sorted(wanted)
        if gene in best_by_gene or gene in seed_hits
    ]
    out = pd.DataFrame(rows)
    qc = {
        "mapped_wanted_genes": int(sum(1 for g in wanted if g in gene_to_proteins)),
        "mapped_seed_genes": int(len(seed_genes)),
        "genes_with_seed_edge": int(len(out)),
    }
    return out, qc


def parse_msigdb(stage: Path, wanted: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    path = stage / "msigdb" / "h.all.v7.5.symbols.gmt"
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] not in PAAD_HALLMARKS:
                continue
            pathway = parts[0]
            for gene in {clean_symbol(x) for x in parts[2:]} & wanted:
                rows.append({"gene": gene, "msigdb_paad_hallmark": pathway})
    out = pd.DataFrame(rows)
    qc = {
        "hallmarks_used": int(out["msigdb_paad_hallmark"].nunique()) if not out.empty else 0,
        "matched_genes": int(out["gene"].nunique()) if not out.empty else 0,
    }
    return out, qc


def parse_sider(stage: Path, candidates: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    names_path = stage / "sider" / "drug_names.tsv"
    se_path = stage / "sider" / "meddra_all_se.tsv.gz"
    name_to_cids: dict[str, set[str]] = defaultdict(set)
    with names_path.open("r", encoding="utf-8", errors="replace") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2:
                name_to_cids[norm_name(row[1])].add(row[0])
    wanted_cids: set[str] = set()
    cand_cids: dict[str, set[str]] = {}
    for row in candidates.itertuples():
        ids = set(name_to_cids.get(norm_name(row.drug_name), set()))
        cand_cids[str(row.canonical_drug_id)] = ids
        wanted_cids.update(ids)
    effects_by_cid: dict[str, set[str]] = defaultdict(set)
    serious_by_cid: dict[str, set[str]] = defaultdict(set)
    if wanted_cids:
        with gzip.open(se_path, "rt", encoding="utf-8", errors="replace") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) < 6:
                    continue
                cids = {row[0], row[1]} & wanted_cids
                if not cids:
                    continue
                effect = row[5]
                low = effect.lower()
                for cid in cids:
                    effects_by_cid[cid].add(effect)
                    if any(k in low for k in SERIOUS_SIDE_EFFECT_KEYWORDS):
                        serious_by_cid[cid].add(effect)
    rows = []
    for row in candidates.itertuples():
        cid_key = str(row.canonical_drug_id)
        ids = cand_cids.get(cid_key, set())
        effects = set().union(*(effects_by_cid.get(cid, set()) for cid in ids)) if ids else set()
        serious = set().union(*(serious_by_cid.get(cid, set()) for cid in ids)) if ids else set()
        rows.append(
            {
                "canonical_drug_id": cid_key,
                "drug_name": row.drug_name,
                "sider_cids": ";".join(sorted(ids)),
                "sider_has_match": bool(ids),
                "sider_side_effect_count": int(len(effects)),
                "sider_serious_keyword_count": int(len(serious)),
                "sider_serious_examples": ";".join(sorted(serious)[:8]),
            }
        )
    out = pd.DataFrame(rows)
    qc = {
        "candidate_matches": int(out["sider_has_match"].sum()) if not out.empty else 0,
        "candidate_count": int(len(out)),
    }
    return out, qc


def aggregate_gene_signature(candidates: pd.DataFrame, sig: pd.DataFrame, prefix: str) -> pd.DataFrame:
    if sig.empty:
        return pd.DataFrame(
            {
                "canonical_drug_id": candidates["canonical_drug_id"].astype(str),
                f"{prefix}_matched_genes": "",
                f"{prefix}_mean_delta": np.nan,
                f"{prefix}_support_fraction": np.nan,
            }
        )
    by_gene = sig.set_index("gene")
    rows = []
    for row in candidates.itertuples():
        genes = split_genes(getattr(row, "target_genes", ""))
        matched = [g for g in genes if g in by_gene.index]
        deltas = [float(by_gene.loc[g, "log2fc_or_delta"]) for g in matched]
        supports = [bool(by_gene.loc[g, "support_up"]) for g in matched]
        rows.append(
            {
                "canonical_drug_id": str(row.canonical_drug_id),
                f"{prefix}_matched_genes": ";".join(matched),
                f"{prefix}_mean_delta": finite_mean(deltas),
                f"{prefix}_support_fraction": float(np.mean(supports)) if supports else np.nan,
            }
        )
    return pd.DataFrame(rows)


def aggregate_gene_table(candidates: pd.DataFrame, table: pd.DataFrame, prefix: str, value_cols: list[str]) -> pd.DataFrame:
    if table.empty:
        return pd.DataFrame({"canonical_drug_id": candidates["canonical_drug_id"].astype(str)})
    by_gene = table.set_index("gene")
    rows = []
    for row in candidates.itertuples():
        genes = split_genes(getattr(row, "target_genes", ""))
        matched = [g for g in genes if g in by_gene.index]
        record: dict[str, Any] = {"canonical_drug_id": str(row.canonical_drug_id), f"{prefix}_matched_genes": ";".join(matched)}
        for col in value_cols:
            vals = [float(by_gene.loc[g, col]) for g in matched if col in by_gene.columns and pd.notna(by_gene.loc[g, col])]
            record[f"{prefix}_{col}"] = finite_mean(vals)
        if prefix == "msigdb":
            pathways = sorted({str(by_gene.loc[g, "msigdb_paad_hallmark"]) for g in matched if "msigdb_paad_hallmark" in by_gene.columns})
            record["msigdb_paad_hallmark_count"] = len(pathways)
            record["msigdb_paad_hallmarks"] = ";".join(pathways)
        rows.append(record)
    return pd.DataFrame(rows)


def build_v2_scores(df: pd.DataFrame) -> pd.DataFrame:
    n = len(df)
    out = df.copy()
    out["model_rank_bonus"] = (n - pd.to_numeric(out["rank"], errors="coerce") + 1) / max(1, n)
    out["tcga_gtex_bonus"] = pd.to_numeric(out["tcga_gtex_mean_delta"], errors="coerce").fillna(0).clip(lower=0, upper=2) / 2
    geo_vals = pd.concat(
        [
            pd.to_numeric(out.get("gse62452_mean_delta"), errors="coerce"),
            pd.to_numeric(out.get("gse71729_mean_delta"), errors="coerce"),
        ],
        axis=1,
    )
    out["geo_reproducibility_bonus"] = geo_vals.clip(lower=0, upper=1.5).mean(axis=1, skipna=True).fillna(0) / 1.5
    out["cptac_proteomics_bonus"] = pd.to_numeric(out["cptac_proteomics_mean_delta"], errors="coerce").fillna(0).clip(lower=0, upper=1.5) / 1.5
    out["rppa_bonus"] = pd.to_numeric(out.get("rppa_rppa_mean_z"), errors="coerce").fillna(0).clip(lower=0, upper=1.0)
    out["string_bonus"] = pd.to_numeric(out.get("string_string_paad_seed_max_score"), errors="coerce").fillna(0).clip(lower=0, upper=1.0)
    out["msigdb_bonus"] = pd.to_numeric(out.get("msigdb_paad_hallmark_count"), errors="coerce").fillna(0).clip(upper=3) / 3
    out["clinical_trial_bonus"] = pd.to_numeric(out.get("clinical_trial_mention_count"), errors="coerce").fillna(0).clip(upper=5) / 5
    out["prism_bonus"] = (
        out.get("prism_has_evidence", False).fillna(False).astype(bool).astype(float)
        * (-pd.to_numeric(out.get("prism_mean_log2fc"), errors="coerce").fillna(0)).clip(lower=0, upper=2)
        / 2
    )
    out["opentargets_bonus"] = pd.to_numeric(out.get("opentargets_max_score"), errors="coerce").fillna(0).clip(lower=0, upper=1)
    out["admet_bonus"] = out.get("admet_category", pd.Series([""] * len(out))).map(
        {"Approved": 0.60, "Candidate": 0.40, "Caution": -0.40, "NO_SMILES": -0.80}
    ).fillna(0.0)
    out["sider_penalty"] = pd.to_numeric(out.get("sider_serious_keyword_count"), errors="coerce").fillna(0).clip(upper=80) / 160
    support_cols = [
        "tcga_gtex_mean_delta",
        "gse62452_mean_delta",
        "gse71729_mean_delta",
        "cptac_proteomics_mean_delta",
        "rppa_rppa_mean_z",
        "string_string_paad_seed_max_score",
    ]
    support = pd.DataFrame(index=out.index)
    for col in support_cols:
        support[col] = pd.to_numeric(out.get(col), errors="coerce").fillna(0) > (0.15 if "string" not in col else 0.50)
    support["msigdb"] = pd.to_numeric(out.get("msigdb_paad_hallmark_count"), errors="coerce").fillna(0) > 0
    support["prism"] = out.get("prism_has_evidence", False).fillna(False).astype(bool)
    support["opentargets"] = pd.to_numeric(out.get("opentargets_evidence_count"), errors="coerce").fillna(0) > 0
    out["external_source_support_count"] = support.sum(axis=1).astype(int)
    out["external_v2_score"] = (
        0.80 * out["model_rank_bonus"]
        + 0.80 * out["tcga_gtex_bonus"]
        + 0.70 * out["geo_reproducibility_bonus"]
        + 0.70 * out["cptac_proteomics_bonus"]
        + 0.35 * out["rppa_bonus"]
        + 0.45 * out["string_bonus"]
        + 0.35 * out["msigdb_bonus"]
        + 0.45 * out["clinical_trial_bonus"]
        + 0.55 * out["prism_bonus"]
        + 0.45 * out["opentargets_bonus"]
        + out["admet_bonus"]
        - out["sider_penalty"]
    )
    out["tier_v2"] = np.select(
        [
            (out["external_v2_score"] >= 3.0) & (out["external_source_support_count"] >= 5) & ~out["admet_category"].eq("Caution"),
            (out["external_v2_score"] >= 2.2) & (out["external_source_support_count"] >= 4),
            out["admet_category"].eq("Caution"),
        ],
        ["Tier 1", "Tier 2", "Excluded"],
        default="Tier 3",
    )
    out["external_v2_category"] = np.select(
        [
            out.get("known_paad_control", False).fillna(False).astype(bool),
            pd.to_numeric(out.get("clinical_trial_mention_count"), errors="coerce").fillna(0).gt(0),
            out["external_source_support_count"].ge(5),
        ],
        [
            "Known PAAD control/support",
            "PAAD clinical-trial supported candidate",
            "Multi-source external evidence candidate",
        ],
        default="Exploratory candidate",
    )
    tier_order = {"Tier 1": 1, "Tier 2": 2, "Tier 3": 3, "Excluded": 4}
    out["_tier_order"] = out["tier_v2"].map(tier_order).fillna(99)
    return out.sort_values(["_tier_order", "external_v2_score", "ensemble_score"], ascending=[True, False, False]).drop(columns=["_tier_order"])


def write_report(out_dir: Path, report_dir: Path, final: pd.DataFrame, qc: dict[str, Any]) -> Path:
    report = report_dir / "PAAD_EXTERNAL_VALIDATION_V2_20260423.md"
    cols = [
        "rank",
        "drug_name",
        "tier_v2",
        "external_v2_score",
        "external_source_support_count",
        "tcga_gtex_mean_delta",
        "gse62452_mean_delta",
        "gse71729_mean_delta",
        "cptac_proteomics_mean_delta",
        "string_string_paad_seed_max_score",
        "msigdb_paad_hallmark_count",
        "sider_serious_keyword_count",
    ]
    top = final[[c for c in cols if c in final.columns]].head(20).copy()
    lines = ["| " + " | ".join(top.columns) + " |", "| " + " | ".join(["---"] * len(top.columns)) + " |"]
    for row in top.itertuples(index=False):
        vals = []
        for val in row:
            if isinstance(val, float):
                vals.append("" if not np.isfinite(val) else f"{val:.4f}")
            else:
                vals.append(str(val))
        lines.append("| " + " | ".join(vals) + " |")
    validation_summary_path = Path(qc.get("outputs", {}).get("knowledge_validation", out_dir)) / "validation_summary_v2.csv"
    text = f"""# PAAD External Validation v2 - 2026-04-23

## Scope

This run revalidates the PAAD GroupCV Top50 candidates with the newly staged external sources:

- GTEx Pancreas vs TCGA-PAAD target expression
- GEO GSE62452 and GSE71729 tumor-normal reproducibility
- CPTAC-PDAC proteomics tumor-normal support
- cBioPortal PAAD RPPA/mutation support
- STRING PAAD seed network proximity
- MSigDB Hallmark pathway overlap
- SIDER side-effect evidence
- Existing PRISM/DepMap, OpenTargets, ClinicalTrials, and ADMET evidence

## QC

```json
{json.dumps(qc, ensure_ascii=False, indent=2)}
```

## Top validated candidates

{chr(10).join(lines)}

## Outputs

- `{out_dir / "top50_external_validation_v2.csv"}`
- `{out_dir / "top15_validated_v2.csv"}`
- `{out_dir / "source_level_gene_evidence.csv"}`
- `{validation_summary_path}`
"""
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(text, encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Run PAAD external validation v2 with GTEx/GEO/CPTAC/SIDER/STRING/MSigDB")
    parser.add_argument("--config", default="config/paad_pipeline_config.json")
    parser.add_argument("--input-set", default=INPUT_NAME)
    parser.add_argument("--variant", default=VARIANT)
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = pipeline_paths(cfg)
    stage = paths.root / "data" / "paad_source_staging"
    candidates = load_candidate_tables(paths, args.input_set)
    candidates["canonical_drug_id"] = candidates["canonical_drug_id"].astype(str)
    wanted = target_universe(candidates, cfg)

    qc: dict[str, Any] = {"candidate_count": int(len(candidates)), "target_universe_genes": int(len(wanted))}

    tcga = load_paad_external_expression(paths.root / "data" / "paad_external" / "paad_expression.csv", wanted)
    gtex = load_gtex_pancreas(stage / "gtex" / "gene_tpm_2017-06-05_v8_pancreas.gct.gz", wanted)
    sig_tcga_gtex = signature_table("TCGA_PAAD_vs_GTEx_pancreas", tcga, gtex, raw_tpm=True)
    qc["tcga_gtex"] = {
        "matched_genes": int(sig_tcga_gtex["gene"].nunique()) if not sig_tcga_gtex.empty else 0,
        "tcga_samples": int(len(next(iter(tcga.values())))) if tcga else 0,
        "gtex_samples": int(len(next(iter(gtex.values())))) if gtex else 0,
    }

    sig_gse62452, qc_gse62452 = parse_gse62452(
        stage / "geo" / "GSE62452" / "GSE62452_series_matrix.txt.gz",
        stage / "geo" / "GPL6244" / "GPL6244.annot.gz",
        wanted,
    )
    sig_gse71729, qc_gse71729 = parse_gse71729(stage / "geo" / "GSE71729" / "GSE71729_series_matrix.txt.gz", wanted)
    qc["gse62452"] = qc_gse62452
    qc["gse71729"] = qc_gse71729

    sig_cptac, qc_cptac = parse_cptac_proteomics(stage, wanted)
    rppa, qc_rppa = parse_rppa(stage, wanted)
    mutations, qc_mut = parse_mutations(stage, wanted)
    string, qc_string = parse_string(stage, wanted)
    msigdb, qc_msigdb = parse_msigdb(stage, wanted)
    sider, qc_sider = parse_sider(stage, candidates)
    qc["cptac_proteomics"] = qc_cptac
    qc["cbioportal_rppa"] = qc_rppa
    qc["cbioportal_mutations"] = qc_mut
    qc["string"] = qc_string
    qc["msigdb"] = qc_msigdb
    qc["sider"] = qc_sider

    merged = candidates.copy()
    for table in [
        aggregate_gene_signature(candidates, sig_tcga_gtex, "tcga_gtex"),
        aggregate_gene_signature(candidates, sig_gse62452, "gse62452"),
        aggregate_gene_signature(candidates, sig_gse71729, "gse71729"),
        aggregate_gene_signature(candidates, sig_cptac, "cptac_proteomics"),
        aggregate_gene_table(candidates, rppa, "rppa", ["rppa_mean_z", "rppa_positive_fraction"]),
        aggregate_gene_table(candidates, mutations, "cbio", ["cbio_mutation_count", "cbio_mutated_samples"]),
        aggregate_gene_table(candidates, string, "string", ["string_paad_seed_max_score"]),
        aggregate_gene_table(candidates, msigdb, "msigdb", []),
        sider,
    ]:
        duplicate_cols = [c for c in table.columns if c != "canonical_drug_id" and c in merged.columns]
        if duplicate_cols:
            table = table.drop(columns=duplicate_cols)
        merged = merged.merge(table, on="canonical_drug_id", how="left")

    cfg_admet = json.loads(json.dumps(cfg))
    cfg_admet["paths"]["admet_output_dir"] = f"admet/paad/{args.variant}"
    cfg_admet["paths"]["reports_dir"] = f"reports/paad/{args.variant}"
    admet_input = merged.head(50).copy()
    admet = admet_assessment(cfg_admet, admet_input)
    admet["canonical_drug_id"] = admet["canonical_drug_id"].astype(str)
    admet_cols = ["canonical_drug_id", "admet_coverage", "toxicity_flags", "low_confidence_toxic_signals", "admet_category"]
    merged = merged.drop(columns=[c for c in admet_cols if c in merged.columns and c != "canonical_drug_id"], errors="ignore")
    merged = merged.merge(admet[[c for c in admet_cols if c in admet.columns]], on="canonical_drug_id", how="left")

    final = build_v2_scores(merged)
    out_dir = paths.external_validation_dir / args.variant
    kg_dir = paths.knowledge_validation_dir / args.variant
    phase_dir = paths.phase5_dir / args.variant
    report_dir = paths.reports_dir / args.variant
    for directory in [out_dir, kg_dir, phase_dir, report_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    gene_evidence = pd.concat(
        [sig_tcga_gtex, sig_gse62452, sig_gse71729, sig_cptac],
        ignore_index=True,
    )
    final.to_csv(out_dir / "top50_external_validation_v2.csv", index=False)
    final.head(15).to_csv(out_dir / "top15_validated_v2.csv", index=False)
    gene_evidence.to_csv(out_dir / "source_level_gene_evidence.csv", index=False)
    final.to_csv(kg_dir / "validation_summary_v2.csv", index=False)
    final.to_csv(phase_dir / "final_comprehensive_candidates.csv", index=False)
    final.loc[final["tier_v2"].eq("Tier 1")].to_csv(phase_dir / "tier1_high_confidence.csv", index=False)

    qc["tier_counts"] = final["tier_v2"].value_counts(dropna=False).to_dict()
    qc["top15_support_mean"] = float(final.head(15)["external_source_support_count"].mean())
    qc["top15_tier1_count"] = int(final.head(15)["tier_v2"].eq("Tier 1").sum())
    qc["outputs"] = {
        "external_validation": str(out_dir),
        "knowledge_validation": str(kg_dir),
        "phase5": str(phase_dir),
        "reports": str(report_dir),
    }
    write_json(report_dir / "qc_paad_external_validation_v2_20260423.json", qc)
    report = write_report(out_dir, report_dir, final, qc)

    summary = {
        "status": "completed",
        "variant": args.variant,
        "report": str(report),
        "top15": final.head(15)[
            ["rank", "drug_name", "tier_v2", "external_v2_score", "external_source_support_count"]
        ].to_dict(orient="records"),
        "qc": qc,
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
