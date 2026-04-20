#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import load_config, pipeline_paths, write_json, write_table

try:
    from rdkit import Chem
except Exception:  # pragma: no cover - optional dependency
    Chem = None


NON_NAME = re.compile(r"[^a-z0-9]+")


def norm_name(value: Any) -> str:
    return NON_NAME.sub("", "" if pd.isna(value) else str(value).lower())


def canonical_smiles(value: Any) -> str:
    if pd.isna(value) or str(value).strip() == "":
        return ""
    text = str(value).strip()
    if Chem is None:
        return text
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def read_parquet(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_parquet(path, columns=columns)


def audit_removed_smiles(staging: Path, reports_dir: Path, source_qc: dict[str, Any]) -> pd.DataFrame:
    removed = pd.DataFrame(source_qc.get("primary_filter_qc", {}).get("removed_drugs", []))
    if removed.empty:
        return removed
    removed["_norm"] = removed["drug_name"].map(norm_name)

    drugbank = read_parquet(
        staging / "drugbank" / "drugbank_drug_master_basic_20260406.parquet",
        columns=["drugbank_id", "name", "smiles", "inchikey"],
    )
    chembl = read_parquet(
        staging / "chembl" / "chembl_compound_master_basic_20260406.parquet",
        columns=["chembl_id", "pref_name", "canonical_smiles", "standard_inchi_key"],
    )
    lincs = read_parquet(
        staging / "lincs" / "lincs_pert_info_basic_20260406.parquet",
        columns=["pert_id", "pert_iname", "canonical_smiles", "inchi_key"],
    )

    for df, col in [(drugbank, "name"), (chembl, "pref_name"), (lincs, "pert_iname")]:
        if not df.empty and col in df.columns:
            df["_norm"] = df[col].map(norm_name)

    rows = []
    for row in removed.to_dict(orient="records"):
        norm = row["_norm"]
        candidates = []
        if not drugbank.empty:
            hit = drugbank.loc[drugbank["_norm"].eq(norm)].head(1)
            if not hit.empty:
                r = hit.iloc[0]
                candidates.append(("drugbank_name", r.get("drugbank_id", ""), r.get("name", ""), r.get("smiles", ""), r.get("inchikey", "")))
        if not chembl.empty:
            hit = chembl.loc[chembl["_norm"].eq(norm)].head(1)
            if not hit.empty:
                r = hit.iloc[0]
                candidates.append(("chembl_pref_name", r.get("chembl_id", ""), r.get("pref_name", ""), r.get("canonical_smiles", ""), r.get("standard_inchi_key", "")))
        if not lincs.empty:
            hit = lincs.loc[lincs["_norm"].eq(norm)].head(1)
            if not hit.empty:
                r = hit.iloc[0]
                candidates.append(("lincs_pert_name", r.get("pert_id", ""), r.get("pert_iname", ""), r.get("canonical_smiles", ""), r.get("inchi_key", "")))
        if not candidates:
            candidates.append(("no_local_name_match", "", "", "", ""))
        for source, source_id, source_name, smiles, inchikey in candidates:
            can = canonical_smiles(smiles)
            rows.append(
                {
                    "canonical_drug_id": row.get("canonical_drug_id", ""),
                    "drug_name": row.get("drug_name", ""),
                    "target_genes": row.get("target_genes", ""),
                    "pathway": row.get("PATHWAY_NAME_NORMALIZED", ""),
                    "candidate_source": source,
                    "candidate_source_id": source_id,
                    "candidate_name": source_name,
                    "candidate_has_smiles": bool(can),
                    "candidate_canonical_smiles": can,
                    "candidate_inchikey": inchikey,
                }
            )
    out = pd.DataFrame(rows)
    write_table(out, reports_dir / "missingness" / "smiles_removed_drug_recovery_review.csv")
    return out


def audit_lincs(drug_features: pd.DataFrame, staging: Path, reports_dir: Path) -> pd.DataFrame:
    lincs_cols = [c for c in drug_features.columns if c.startswith("drug__lincs__")]
    has_lincs = drug_features[lincs_cols].notna().any(axis=1) if lincs_cols else pd.Series(False, index=drug_features.index)
    review = drug_features[["canonical_drug_id", "drug_name", "canonical_smiles", "target_genes", "classification"]].copy()
    review["has_lincs_signature"] = has_lincs.values
    review["_drug_norm"] = review["drug_name"].map(norm_name)
    review["_canonical_smiles_norm"] = review["canonical_smiles"].map(canonical_smiles)

    lincs = read_parquet(
        staging / "lincs" / "lincs_pert_info_basic_20260406.parquet",
        columns=["pert_id", "pert_iname", "canonical_smiles", "inchi_key", "inchi_key_prefix", "pubchem_cid"],
    )
    if not lincs.empty:
        lincs["_drug_norm"] = lincs["pert_iname"].map(norm_name)
        lincs["_canonical_smiles_norm"] = lincs["canonical_smiles"].map(canonical_smiles)
        name_hits = set(lincs["_drug_norm"].dropna())
        smiles_hits = set(lincs["_canonical_smiles_norm"].dropna())
        review["local_lincs_name_match_available"] = review["_drug_norm"].isin(name_hits)
        review["local_lincs_smiles_match_available"] = review["_canonical_smiles_norm"].isin(smiles_hits) & review["_canonical_smiles_norm"].ne("")
    else:
        review["local_lincs_name_match_available"] = False
        review["local_lincs_smiles_match_available"] = False

    review["suggested_action"] = "keep_existing_lincs"
    needs = ~review["has_lincs_signature"]
    review.loc[needs & review["local_lincs_name_match_available"], "suggested_action"] = "recover_by_lincs_name_bridge"
    review.loc[needs & review["local_lincs_smiles_match_available"], "suggested_action"] = "recover_by_lincs_smiles_bridge"
    review.loc[needs & ~review["local_lincs_name_match_available"] & ~review["local_lincs_smiles_match_available"], "suggested_action"] = "needs_external_lincs_or_analog_imputation"
    review = review.drop(columns=["_drug_norm", "_canonical_smiles_norm"])
    write_table(review, reports_dir / "missingness" / "lincs_drug_mapping_review.csv")
    return review


def audit_cell_lines(sample_features: pd.DataFrame, staging: Path, reports_dir: Path) -> pd.DataFrame:
    model = read_parquet(
        staging / "depmap" / "depmap_model_basic_clean_20260406.parquet",
        columns=["ModelID", "CellLineName", "StrippedCellLineName", "SangerModelID", "COSMICID", "OncotreeSubtype"],
    )
    crispr = read_parquet(staging / "depmap" / "depmap_crispr_gene_dependency_basic_clean_20260406.parquet", columns=["ModelID"])
    crispr_ids = set(crispr["ModelID"]) if not crispr.empty else set()

    review = sample_features[["sample_id", "cell_line_name", "sample_has_crispr"]].copy()
    if not model.empty:
        model = model.drop_duplicates("SangerModelID")
        review = review.merge(
            model,
            left_on="sample_id",
            right_on="SangerModelID",
            how="left",
        )
        review["depmap_model_available"] = review["ModelID"].notna()
        review["depmap_crispr_model_available"] = review["ModelID"].isin(crispr_ids)
    else:
        review["depmap_model_available"] = False
        review["depmap_crispr_model_available"] = False

    review["suggested_action"] = "keep_crispr_features"
    missing = ~review["sample_has_crispr"].astype(bool)
    review.loc[missing & review["depmap_model_available"], "suggested_action"] = "add_expression_cnv_mutation_fallback_or_latest_depmap_check"
    review.loc[missing & ~review["depmap_model_available"], "suggested_action"] = "resolve_alias_with_sanger_cosmic_cell_model_passports"
    write_table(review, reports_dir / "missingness" / "cellline_feature_coverage_review.csv")
    return review


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit thyroid missingness and local source recovery options")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--staging-dir", default="data/source_staging")
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = pipeline_paths(cfg)
    staging = Path(args.staging_dir)
    if not staging.is_absolute():
        staging = paths.root / staging

    reports_missing = paths.reports_dir / "missingness"
    reports_missing.mkdir(parents=True, exist_ok=True)
    drug_features = pd.read_csv(paths.raw_dir / "drug_features.csv")
    sample_features = pd.read_csv(paths.raw_dir / "sample_features.csv")
    source_qc_path = paths.reports_dir / "qc_source_to_model_ready_20260420.json"
    source_qc = json.loads(source_qc_path.read_text()) if source_qc_path.exists() else {}

    smiles_review = audit_removed_smiles(staging, paths.reports_dir, source_qc)
    lincs_review = audit_lincs(drug_features, staging, paths.reports_dir)
    cell_review = audit_cell_lines(sample_features, staging, paths.reports_dir)

    summary = {
        "step": "missingness_source_audit",
        "smiles_removed_drugs": int(source_qc.get("primary_filter_qc", {}).get("removed", {}).get("drugs", 0)),
        "smiles_removed_recovery_candidates_with_local_smiles": int(smiles_review.get("candidate_has_smiles", pd.Series(dtype=bool)).sum()) if not smiles_review.empty else 0,
        "lincs_drugs_total": int(len(lincs_review)),
        "lincs_direct_signature_drugs": int(lincs_review["has_lincs_signature"].sum()),
        "lincs_local_name_or_smiles_bridge_candidates": int(
            ((~lincs_review["has_lincs_signature"]) & (lincs_review["local_lincs_name_match_available"] | lincs_review["local_lincs_smiles_match_available"])).sum()
        ),
        "cell_lines_total": int(len(cell_review)),
        "cell_lines_with_crispr": int(cell_review["sample_has_crispr"].sum()),
        "cell_lines_missing_crispr_but_depmap_model_available": int(
            ((~cell_review["sample_has_crispr"].astype(bool)) & cell_review["depmap_model_available"]).sum()
        ),
        "outputs": {
            "smiles_removed_drug_recovery_review": str(reports_missing / "smiles_removed_drug_recovery_review.csv"),
            "lincs_drug_mapping_review": str(reports_missing / "lincs_drug_mapping_review.csv"),
            "cellline_feature_coverage_review": str(reports_missing / "cellline_feature_coverage_review.csv"),
        },
    }
    write_json(reports_missing / "missingness_source_audit_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
