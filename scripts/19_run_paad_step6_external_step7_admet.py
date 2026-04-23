#!/usr/bin/env python3
"""PAAD Step 6 external-cohort validation and separate Step 7 ADMET gate.

This mirrors the BRCA METABRIC flow:
  Step 6: validate model-recommended Top30 drugs in external PAAD cohorts.
  Step 7: run ADMET only after Step 6 top15 selection.

Unlike scripts/18_run_paad_external_validation_v2.py, Step 6 here does not
mix ADMET/SIDER into the external validation counts.
"""

from __future__ import annotations

import argparse
import json
import sys
from copy import deepcopy
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from thyroid_pipeline.core import admet_assessment, load_config, pipeline_paths


DEFAULT_INPUT_SET = "numeric_strong_context_smiles_pan_lincs"
DEFAULT_CV = "groupcv4_drug"
DEFAULT_V2_VARIANT = "groupcv4_drug_v2_sources"
DEFAULT_STEP6_VARIANT = "groupcv4_drug_step6_external_cohort"
DEFAULT_STEP7_VARIANT = "groupcv4_drug_step7_admet_after_external"
STAMP = "20260423"


def safe_float(value: Any) -> float | None:
    try:
        out = float(value)
    except Exception:
        return None
    return out if np.isfinite(out) else None


def json_safe(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    if isinstance(obj, tuple):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return safe_float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if pd.isna(obj) if not isinstance(obj, (list, dict, tuple)) else False:
        return None
    return obj


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(payload), ensure_ascii=False, indent=2), encoding="utf-8")


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows_"
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for row in df.itertuples(index=False):
        values = []
        for value in row:
            if isinstance(value, float):
                values.append("" if not np.isfinite(value) else f"{value:.4f}")
            else:
                text = "" if pd.isna(value) else str(value)
                values.append(text.replace("|", "/"))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def positive(series: pd.Series, threshold: float = 0.0) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(-np.inf).gt(threshold)


def nonempty(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().ne("").astype(bool)


def pathway_relevant(row: pd.Series, biology_terms: list[str]) -> bool:
    haystack = " ".join(
        str(row.get(col, ""))
        for col in ["drug_name", "target_genes", "target_match_genes", "PATHWAY_NAME_NORMALIZED", "classification"]
    ).upper()
    return any(term.upper() in haystack for term in biology_terms)


def known_precision(df: pd.DataFrame, known_names: set[str]) -> dict[str, dict[str, float | int]]:
    out: dict[str, dict[str, float | int]] = {}
    names = df.sort_values("rank")["drug_name"].astype(str).str.upper().tolist()
    for k in [5, 10, 15, 20, 25, 30]:
        top_k = names[: min(k, len(names))]
        hits = int(sum(name in known_names for name in top_k))
        total = int(len(top_k))
        out[f"P@{k}"] = {
            "precision": float(hits / total) if total else 0.0,
            "hits": hits,
            "total": total,
        }
    return out


def build_step6(
    cfg: dict[str, Any],
    top50_path: Path,
    out_dir: Path,
    report_dir: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    df = pd.read_csv(top50_path)
    top30 = df.sort_values("rank").head(30).copy()
    top30["canonical_drug_id"] = top30["canonical_drug_id"].astype(str)
    top30["drug_id"] = top30["canonical_drug_id"]

    top30["target_mapped"] = nonempty(top30.get("target_match_genes", pd.Series("", index=top30.index)))
    top30["tcga_gtex_target_up"] = positive(top30.get("tcga_gtex_mean_delta", pd.Series(np.nan, index=top30.index)))
    top30["gse62452_target_up"] = positive(top30.get("gse62452_mean_delta", pd.Series(np.nan, index=top30.index)))
    top30["gse71729_target_up"] = positive(top30.get("gse71729_mean_delta", pd.Series(np.nan, index=top30.index)))
    top30["cptac_protein_up"] = positive(top30.get("cptac_proteomics_mean_delta", pd.Series(np.nan, index=top30.index)))
    top30["geo_external_validated"] = top30["gse62452_target_up"] | top30["gse71729_target_up"]
    top30["both_geo_validated"] = top30["gse62452_target_up"] & top30["gse71729_target_up"]
    top30["external_cohort_validated"] = top30["geo_external_validated"] | top30["cptac_protein_up"]
    top30["any_disease_reference_validated"] = top30["external_cohort_validated"] | top30["tcga_gtex_target_up"]

    top30["survival_p"] = pd.to_numeric(top30.get("survival_p_value"), errors="coerce")
    top30["survival_sig"] = top30["survival_p"].lt(0.05).fillna(False)
    top30["survival_risk_direction"] = top30.get("survival_direction", "").astype(str).str.contains("shorter", case=False, na=False)
    top30["survival_risk_sig"] = top30["survival_sig"] & top30["survival_risk_direction"]

    known_names = {str(x).upper() for x in cfg.get("known_paad_drugs", [])}
    top30["known_paad"] = top30["drug_name"].astype(str).str.upper().isin(known_names) | top30.get(
        "known_paad_control", False
    ).fillna(False).astype(bool)
    top30["paad_biology_relevant"] = top30.apply(lambda row: pathway_relevant(row, cfg.get("paad_biology_terms", [])), axis=1)
    top30["prism_evidence"] = top30.get("prism_has_evidence", False).fillna(False).astype(bool)
    top30["high_prism_sensitivity"] = (
        pd.to_numeric(top30.get("prism_sensitive_fraction_lt_minus1", 0), errors="coerce").fillna(0).ge(0.25)
    )

    top30["validation_score"] = (
        top30["external_cohort_validated"].astype(float) * 2.0
        + top30["both_geo_validated"].astype(float) * 0.75
        + top30["cptac_protein_up"].astype(float) * 1.0
        + top30["survival_sig"].astype(float) * 2.0
        + top30["known_paad"].astype(float) * 2.0
        + top30["paad_biology_relevant"].astype(float) * 1.0
        + top30["prism_evidence"].astype(float) * 1.0
        + top30["high_prism_sensitivity"].astype(float) * 0.5
        - pd.to_numeric(top30["rank"], errors="coerce").fillna(30) * 0.03
    )

    top15 = (
        top30.sort_values(["validation_score", "rank"], ascending=[False, True])
        .drop_duplicates("drug_name", keep="first")
        .head(15)
        .sort_values("rank")
        .copy()
    )
    top15["final_rank"] = range(1, len(top15) + 1)

    cols = [
        "final_rank",
        "rank",
        "drug_id",
        "canonical_drug_id",
        "drug_name",
        "canonical_smiles",
        "target_genes",
        "target_match_genes",
        "PATHWAY_NAME_NORMALIZED",
        "classification",
        "mean_pred_ln_ic50",
        "ensemble_score",
        "target_mapped",
        "tcga_gtex_target_up",
        "gse62452_target_up",
        "gse71729_target_up",
        "cptac_protein_up",
        "geo_external_validated",
        "both_geo_validated",
        "external_cohort_validated",
        "any_disease_reference_validated",
        "survival_sig",
        "survival_p",
        "survival_direction",
        "known_paad",
        "paad_biology_relevant",
        "prism_evidence",
        "high_prism_sensitivity",
        "validation_score",
    ]
    score_cols = [c for c in cols if c in top30.columns]

    method_a_details = top30[
        [
            "rank",
            "drug_name",
            "target_match_genes",
            "target_mapped",
            "tcga_gtex_target_up",
            "gse62452_target_up",
            "gse71729_target_up",
            "cptac_protein_up",
            "geo_external_validated",
            "external_cohort_validated",
            "any_disease_reference_validated",
            "tcga_gtex_mean_delta",
            "gse62452_mean_delta",
            "gse71729_mean_delta",
            "cptac_proteomics_mean_delta",
        ]
    ].copy()
    method_b_details = top30[
        [
            "rank",
            "drug_name",
            "survival_sig",
            "survival_risk_sig",
            "survival_p",
            "survival_direction",
        ]
    ].copy()

    p_at_k = known_precision(top30, known_names)
    summary = {
        "step": 6,
        "description": "PAAD external-cohort validation (BRCA/METABRIC-style A+B+C)",
        "source_input": str(top50_path),
        "n_total": int(len(top30)),
        "method_a": {
            "name": "External target expression/protein validation",
            "primary_definition": "external_cohort_validated = GSE62452_up OR GSE71729_up OR CPTAC_protein_up",
            "n_target_mapped": int(top30["target_mapped"].sum()),
            "n_external_cohort_validated": int(top30["external_cohort_validated"].sum()),
            "n_geo_external_validated": int(top30["geo_external_validated"].sum()),
            "n_both_geo_validated": int(top30["both_geo_validated"].sum()),
            "n_any_disease_reference_validated": int(top30["any_disease_reference_validated"].sum()),
            "per_source": {
                "tcga_gtex_target_up": int(top30["tcga_gtex_target_up"].sum()),
                "gse62452_target_up": int(top30["gse62452_target_up"].sum()),
                "gse71729_target_up": int(top30["gse71729_target_up"].sum()),
                "cptac_protein_up": int(top30["cptac_protein_up"].sum()),
            },
            "details": method_a_details.to_dict(orient="records"),
        },
        "method_b": {
            "name": "Target-expression survival stratification",
            "significance_cutoff": 0.05,
            "n_significant": int(top30["survival_sig"].sum()),
            "n_risk_direction_significant": int(top30["survival_risk_sig"].sum()),
            "details": method_b_details.to_dict(orient="records"),
        },
        "method_c": {
            "name": "Known PAAD drug precision (P@K)",
            "precision_at_k": p_at_k,
            "known_paad_drugs": sorted(known_names),
        },
        "prism_crosscheck": {
            "n_prism_evidence": int(top30["prism_evidence"].sum()),
            "n_high_prism_sensitivity": int(top30["high_prism_sensitivity"].sum()),
        },
        "top15_validated": top15[score_cols].to_dict(orient="records"),
        "all_30_scores": top30[score_cols].to_dict(orient="records"),
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    top30[score_cols].to_csv(out_dir / "all_30_scores.csv", index=False)
    top15[score_cols].to_csv(out_dir / "top15_validated.csv", index=False)
    method_a_details.to_csv(out_dir / "method_a_target_expression_validation.csv", index=False)
    method_b_details.to_csv(out_dir / "method_b_survival_validation.csv", index=False)
    pd.DataFrame(
        [{"metric": key, **value} for key, value in p_at_k.items()]
    ).to_csv(out_dir / "method_c_known_paad_precision.csv", index=False)
    write_json(out_dir / "step6_external_cohort_results.json", summary)

    write_step6_report(report_dir / f"PAAD_STEP6_EXTERNAL_COHORT_VALIDATION_{STAMP}.md", summary, top15[score_cols])
    return top15[score_cols], summary


def write_step6_report(path: Path, summary: dict[str, Any], top15: pd.DataFrame) -> None:
    method_a = summary["method_a"]
    method_b = summary["method_b"]
    method_c = summary["method_c"]
    lines = [
        f"# PAAD Step 6 External Cohort Validation - {date.today().isoformat()}",
        "",
        "This is the BRCA/METABRIC-style validation: count how many model Top30 candidates are supported in external PAAD cohorts before ADMET.",
        "",
        "## Summary",
        "",
        f"- Primary external target/protein validation: {method_a['n_external_cohort_validated']}/{summary['n_total']}",
        f"- GEO external validation: {method_a['n_geo_external_validated']}/{summary['n_total']}",
        f"- Both GEO cohorts: {method_a['n_both_geo_validated']}/{summary['n_total']}",
        f"- Any disease reference including TCGA-GTEx: {method_a['n_any_disease_reference_validated']}/{summary['n_total']}",
        f"- Survival significant, p<0.05: {method_b['n_significant']}/{summary['n_total']}",
        f"- Survival significant with high-expression risk direction: {method_b['n_risk_direction_significant']}/{summary['n_total']}",
        f"- PRISM evidence: {summary['prism_crosscheck']['n_prism_evidence']}/{summary['n_total']}",
        "",
        "## Per-Source Method A",
        "",
    ]
    for key, value in method_a["per_source"].items():
        lines.append(f"- {key}: {value}/{summary['n_total']}")
    lines.extend(["", "## Method C P@K", ""])
    for key, value in method_c["precision_at_k"].items():
        lines.append(f"- {key}: {value['hits']}/{value['total']} ({value['precision']:.2%})")
    lines.extend(["", "## Top15 for Step 7 ADMET", ""])
    display = [
        "final_rank",
        "rank",
        "drug_name",
        "external_cohort_validated",
        "geo_external_validated",
        "cptac_protein_up",
        "survival_sig",
        "known_paad",
        "validation_score",
    ]
    lines.append(markdown_table(top15[[c for c in display if c in top15.columns]]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_step7_admet(
    cfg: dict[str, Any],
    step6_top15: pd.DataFrame,
    step6_variant: str,
    step7_variant: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    cfg2 = deepcopy(cfg)
    cfg2["paths"]["admet_output_dir"] = f"admet/paad/{step7_variant}"
    cfg2["paths"]["reports_dir"] = f"reports/paad/{step7_variant}"
    paths2 = pipeline_paths(cfg2)
    admet = admet_assessment(cfg2, step6_top15)
    admet["canonical_drug_id"] = admet["canonical_drug_id"].astype(str)
    merged = step6_top15.merge(
        admet[
            [
                "canonical_drug_id",
                "admet_coverage",
                "admet_no_match_assays",
                "toxicity_flags",
                "low_confidence_toxic_signals",
                "admet_category",
            ]
        ],
        on="canonical_drug_id",
        how="left",
    )
    merged["safety_score"] = (
        10.0
        * pd.to_numeric(merged["admet_coverage"], errors="coerce").fillna(0)
        - merged["toxicity_flags"].fillna("").astype(str).str.len().gt(0).astype(float) * 3.0
        - merged["low_confidence_toxic_signals"].fillna("").astype(str).str.len().gt(0).astype(float) * 1.0
    ).round(3)
    merged["step7_combined_score"] = (
        pd.to_numeric(merged["validation_score"], errors="coerce").fillna(0)
        + merged["admet_category"].map({"Approved": 2.0, "Candidate": 1.0, "Caution": -2.0, "NO_SMILES": -3.0}).fillna(0)
        + merged["safety_score"] / 10.0
    ).round(3)
    merged = merged.sort_values(["step7_combined_score", "validation_score"], ascending=[False, False]).copy()
    merged["step7_rank"] = range(1, len(merged) + 1)
    merged.to_csv(paths2.admet_output_dir / "step7_admet_ranked_candidates.csv", index=False)

    admet_summary_path = paths2.admet_output_dir / "admet_summary.json"
    admet_summary = json.loads(admet_summary_path.read_text(encoding="utf-8")) if admet_summary_path.exists() else {}
    step7_summary = {
        "step": 7,
        "description": "PAAD ADMET gate after Step 6 external-cohort validation",
        "step6_variant": step6_variant,
        "step7_variant": step7_variant,
        "n_drugs_input": int(len(step6_top15)),
        "n_drugs_output": int(len(merged)),
        "admet_summary": admet_summary,
        "final_candidates": merged.to_dict(orient="records"),
        "outputs": {
            "admet_dir": str(paths2.admet_output_dir),
            "reports_dir": str(paths2.reports_dir),
        },
    }
    write_json(paths2.admet_output_dir / "step7_admet_results.json", step7_summary)
    write_step7_report(paths2.reports_dir / f"PAAD_STEP7_ADMET_AFTER_EXTERNAL_{STAMP}.md", step7_summary, merged)
    return merged, step7_summary


def write_step7_report(path: Path, summary: dict[str, Any], ranked: pd.DataFrame) -> None:
    counts = summary.get("admet_summary", {}).get("category_counts", {})
    lines = [
        f"# PAAD Step 7 ADMET After External Validation - {date.today().isoformat()}",
        "",
        "This step mirrors BRCA Step 7: ADMET is applied after the Step 6 external-cohort Top15 selection.",
        "",
        "## ADMET Summary",
        "",
        f"- Input drugs: {summary['n_drugs_input']}",
        f"- Output drugs: {summary['n_drugs_output']}",
        f"- Category counts: {counts}",
        "",
        "## Ranked Candidates",
        "",
    ]
    display = [
        "step7_rank",
        "rank",
        "drug_name",
        "validation_score",
        "admet_coverage",
        "toxicity_flags",
        "low_confidence_toxic_signals",
        "admet_category",
        "safety_score",
        "step7_combined_score",
    ]
    lines.append(markdown_table(ranked[[c for c in display if c in ranked.columns]]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run PAAD BRCA-style Step6 external validation and separate Step7 ADMET")
    parser.add_argument("--config", default="config/paad_pipeline_config.json")
    parser.add_argument("--input-set", default=DEFAULT_INPUT_SET)
    parser.add_argument("--cv", default=DEFAULT_CV)
    parser.add_argument("--v2-variant", default=DEFAULT_V2_VARIANT)
    parser.add_argument("--step6-variant", default=DEFAULT_STEP6_VARIANT)
    parser.add_argument("--step7-variant", default=DEFAULT_STEP7_VARIANT)
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = pipeline_paths(cfg)
    top50_path = paths.external_validation_dir / args.v2_variant / "top50_external_validation_v2.csv"
    if not top50_path.exists():
        raise FileNotFoundError(f"Missing v2 external validation file: {top50_path}")

    step6_dir = paths.external_validation_dir / args.step6_variant
    step6_report_dir = paths.reports_dir / args.step6_variant
    step6_top15, step6_summary = build_step6(cfg, top50_path, step6_dir, step6_report_dir)
    step7_ranked, step7_summary = run_step7_admet(cfg, step6_top15, args.step6_variant, args.step7_variant)

    run_summary = {
        "status": "completed",
        "step6_variant": args.step6_variant,
        "step7_variant": args.step7_variant,
        "step6_primary_external_validated": f"{step6_summary['method_a']['n_external_cohort_validated']}/{step6_summary['n_total']}",
        "step6_any_disease_reference_validated": f"{step6_summary['method_a']['n_any_disease_reference_validated']}/{step6_summary['n_total']}",
        "step7_category_counts": step7_summary.get("admet_summary", {}).get("category_counts", {}),
        "outputs": {
            "step6": str(step6_dir),
            "step6_report": str(step6_report_dir),
            "step7": step7_summary["outputs"]["admet_dir"],
            "step7_report": step7_summary["outputs"]["reports_dir"],
        },
        "top_step7": step7_ranked[["step7_rank", "drug_name", "admet_category", "step7_combined_score"]]
        .head(10)
        .to_dict(orient="records"),
    }
    write_json(paths.reports_dir / f"paad_step6_step7_separated_run_summary_{STAMP}.json", run_summary)
    print(json.dumps(json_safe(run_summary), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
