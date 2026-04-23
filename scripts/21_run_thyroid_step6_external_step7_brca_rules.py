#!/usr/bin/env python3
"""Run thyroid Step 6 external validation and Step 7 ADMET with BRCA rules.

This mirrors the separated PAAD/BRCA-style flow:
  Step 6: validate Top30 thyroid model outputs in the external THCA cohort
  Step 7: run ADMET only after Step 6 top15 selection, using BRCA rules
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

from thyroid_pipeline.core import (  # noqa: E402
    _fingerprint,
    _load_assay,
    _nearest_admet_match,
    canonicalize_smiles,
    load_config,
    pipeline_paths,
)


DEFAULT_SOURCE_VARIANT = "mixed_groupcv_pan_lincs"
DEFAULT_TOP30 = "results/mixed_groupcv/crossattention_residualmlp_lightgbm_groupcv/mixed_groupcv_ensemble_top30_drugs.csv"
DEFAULT_STEP6_VARIANT = "mixed_groupcv_pan_lincs_step6_external_cohort"
DEFAULT_STEP7_VARIANT = "mixed_groupcv_pan_lincs_step7_admet_brca_rules"
STAMP = "20260423"

ADMET_ASSAYS = {
    "caco2_wang": {"category": "Absorption", "name": "Caco-2 Permeability", "type": "regression", "good_direction": "high", "unit": "log(cm/s)", "threshold": -5.15},
    "hia_hou": {"category": "Absorption", "name": "HIA (Human Intestinal Absorption)", "type": "binary", "good_value": 1},
    "pgp_broccatelli": {"category": "Absorption", "name": "P-gp Inhibitor", "type": "binary", "good_value": 0},
    "bioavailability_ma": {"category": "Absorption", "name": "Oral Bioavailability (F>20%)", "type": "binary", "good_value": 1},
    "bbb_martins": {"category": "Distribution", "name": "BBB Penetration", "type": "binary", "good_value": None},
    "ppbr_az": {"category": "Distribution", "name": "Plasma Protein Binding Rate", "type": "regression", "good_direction": "low", "unit": "%", "threshold": 90},
    "vdss_lombardo": {"category": "Distribution", "name": "Volume of Distribution", "type": "regression", "good_direction": None, "unit": "L/kg"},
    "cyp2c9_veith": {"category": "Metabolism", "name": "CYP2C9 Inhibitor", "type": "binary", "good_value": 0},
    "cyp2d6_veith": {"category": "Metabolism", "name": "CYP2D6 Inhibitor", "type": "binary", "good_value": 0},
    "cyp3a4_veith": {"category": "Metabolism", "name": "CYP3A4 Inhibitor", "type": "binary", "good_value": 0},
    "cyp2c9_substrate_carbonmangels": {"category": "Metabolism", "name": "CYP2C9 Substrate", "type": "binary", "good_value": None},
    "cyp2d6_substrate_carbonmangels": {"category": "Metabolism", "name": "CYP2D6 Substrate", "type": "binary", "good_value": None},
    "cyp3a4_substrate_carbonmangels": {"category": "Metabolism", "name": "CYP3A4 Substrate", "type": "binary", "good_value": None},
    "clearance_hepatocyte_az": {"category": "Excretion", "name": "Hepatocyte Clearance", "type": "regression", "good_direction": None, "unit": "uL/min/10^6 cells"},
    "clearance_microsome_az": {"category": "Excretion", "name": "Microsome Clearance", "type": "regression", "good_direction": None, "unit": "mL/min/g"},
    "half_life_obach": {"category": "Excretion", "name": "Half-Life", "type": "regression", "good_direction": "high", "unit": "hr", "threshold": 3},
    "ames": {"category": "Toxicity", "name": "Ames Mutagenicity", "type": "binary", "good_value": 0},
    "dili": {"category": "Toxicity", "name": "DILI (Drug-Induced Liver Injury)", "type": "binary", "good_value": 0},
    "herg": {"category": "Toxicity", "name": "hERG Cardiotoxicity", "type": "binary", "good_value": 0},
    "ld50_zhu": {"category": "Toxicity", "name": "Acute Toxicity (LD50)", "type": "regression", "good_direction": "high", "unit": "log(mol/kg)"},
    "lipophilicity_astrazeneca": {"category": "Properties", "name": "Lipophilicity (logD)", "type": "regression", "good_direction": None, "unit": "logD", "ideal_range": (-0.4, 5.6)},
    "solubility_aqsoldb": {"category": "Properties", "name": "Aqueous Solubility", "type": "regression", "good_direction": "high", "unit": "logS"},
}

KNOWN_APPROVED = {
    "Docetaxel",
    "Paclitaxel",
    "Vinblastine",
    "Vinorelbine",
    "Rapamycin",
    "Bortezomib",
    "Romidepsin",
    "Dactinomycin",
}


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
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    try:
        if pd.isna(obj):
            return None
    except Exception:
        pass
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
            elif isinstance(value, (list, tuple, set, np.ndarray)):
                items = [str(item) for item in value if not pd.isna(item)]
                values.append(", ".join(items).replace("|", "/"))
            else:
                is_missing = bool(pd.isna(value)) if np.isscalar(value) else False
                text = "" if is_missing else str(value)
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


def split_tokens(value: Any) -> list[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    return [token for token in text.split(";") if token]


def gene_match_rate(target_genes: Any, matched_genes: Any) -> float | None:
    targets = split_tokens(target_genes)
    matched = split_tokens(matched_genes)
    if not targets:
        return None
    return float(len(matched) / len(targets))


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


def load_top30(top30_path: Path) -> pd.DataFrame:
    df = pd.read_csv(top30_path).copy()
    df["canonical_drug_id"] = df["canonical_drug_id"].astype(str)
    if "ensemble_score" not in df.columns and "mixed_groupcv_score" in df.columns:
        df["ensemble_score"] = pd.to_numeric(df["mixed_groupcv_score"], errors="coerce")
    if "drug_id" not in df.columns:
        df["drug_id"] = df["canonical_drug_id"]
    return df


def load_external_tables(cfg: dict[str, Any], source_variant: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = pipeline_paths(cfg)
    expr_path = paths.external_validation_dir / source_variant / "thyroid_target_expression.csv"
    surv_path = paths.external_validation_dir / source_variant / "thyroid_survival_validation.csv"
    expr_df = pd.read_csv(expr_path) if expr_path.exists() else pd.DataFrame()
    surv_df = pd.read_csv(surv_path) if surv_path.exists() else pd.DataFrame()
    for df in [expr_df, surv_df]:
        if not df.empty and "canonical_drug_id" in df.columns:
            df["canonical_drug_id"] = df["canonical_drug_id"].astype(str)
    return expr_df, surv_df


def build_step6(
    cfg: dict[str, Any],
    top30_path: Path,
    source_variant: str,
    step6_variant: str,
) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Path]]:
    paths = pipeline_paths(cfg)
    out_dir = paths.external_validation_dir / step6_variant
    report_dir = paths.reports_dir / step6_variant
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    top30 = load_top30(top30_path)
    expr_df, surv_df = load_external_tables(cfg, source_variant)

    if not expr_df.empty:
        expr_keep = [c for c in ["canonical_drug_id", "target_match_genes", "target_expression_pct", "target_expressed"] if c in expr_df.columns]
        top30 = top30.merge(expr_df[expr_keep], on="canonical_drug_id", how="left")
    if not surv_df.empty:
        surv_keep = [c for c in ["canonical_drug_id", "survival_p_value", "survival_direction"] if c in surv_df.columns]
        top30 = top30.merge(surv_df[surv_keep], on="canonical_drug_id", how="left")

    top30["target_mapped"] = nonempty(top30.get("target_match_genes", pd.Series("", index=top30.index)))
    top30["target_gene_match_rate"] = [
        gene_match_rate(targets, matched)
        for targets, matched in zip(top30.get("target_genes", pd.Series("", index=top30.index)), top30.get("target_match_genes", pd.Series("", index=top30.index)))
    ]
    top30["target_expression_pct"] = pd.to_numeric(top30.get("target_expression_pct"), errors="coerce")
    top30["target_expressed"] = top30.get("target_expressed", False).fillna(False).astype(bool)
    top30["external_cohort_validated"] = top30["target_expressed"]

    top30["survival_p"] = pd.to_numeric(top30.get("survival_p_value"), errors="coerce")
    top30["survival_sig"] = top30["survival_p"].lt(0.05).fillna(False)
    top30["survival_risk_direction"] = top30.get("survival_direction", "").astype(str).str.contains("shorter", case=False, na=False)
    top30["survival_risk_sig"] = top30["survival_sig"] & top30["survival_risk_direction"]

    known_names = {str(x).upper() for x in cfg.get("known_thyroid_drugs", [])}
    top30["known_thyroid"] = top30["drug_name"].astype(str).str.upper().isin(known_names)
    top30["thyroid_biology_relevant"] = top30.apply(lambda row: pathway_relevant(row, cfg.get("thyroid_biology_terms", [])), axis=1)

    top30["validation_score"] = (
        top30["external_cohort_validated"].astype(float) * 2.0
        + top30["target_expression_pct"].fillna(0.0).clip(lower=0.0, upper=1.0) * 1.0
        + top30["survival_sig"].astype(float) * 2.0
        + top30["survival_risk_sig"].astype(float) * 0.5
        + top30["known_thyroid"].astype(float) * 2.0
        + top30["thyroid_biology_relevant"].astype(float) * 1.0
        - pd.to_numeric(top30["rank"], errors="coerce").fillna(30) * 0.03
    )

    top15 = (
        top30.sort_values(["validation_score", "rank"], ascending=[False, True])
        .drop_duplicates("drug_name", keep="first")
        .head(15)
        .copy()
    )
    top15["final_rank"] = range(1, len(top15) + 1)

    score_cols = [
        "rank",
        "drug_id",
        "canonical_drug_id",
        "drug_name",
        "canonical_smiles",
        "target_genes",
        "target_match_genes",
        "target_gene_match_rate",
        "PATHWAY_NAME_NORMALIZED",
        "classification",
        "mean_pred_ln_ic50",
        "ensemble_score",
        "target_expression_pct",
        "target_expressed",
        "external_cohort_validated",
        "survival_sig",
        "survival_risk_sig",
        "survival_p",
        "survival_direction",
        "known_thyroid",
        "thyroid_biology_relevant",
        "validation_score",
    ]
    score_cols = [c for c in score_cols if c in top30.columns]

    method_a = top30[
        [
            "rank",
            "drug_name",
            "target_match_genes",
            "target_gene_match_rate",
            "target_expression_pct",
            "target_expressed",
            "external_cohort_validated",
        ]
    ].copy()
    method_b = top30[
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
        "description": "Thyroid external-cohort validation (BRCA-style A+B+C)",
        "source_variant": source_variant,
        "source_input": str(top30_path),
        "n_total": int(len(top30)),
        "method_a": {
            "name": "External target expression validation",
            "primary_definition": "external_cohort_validated = target_expressed in TCGA-THCA",
            "n_target_mapped": int(top30["target_mapped"].sum()),
            "n_external_cohort_validated": int(top30["external_cohort_validated"].sum()),
            "details": method_a.to_dict(orient="records"),
        },
        "method_b": {
            "name": "Survival stratification",
            "n_survival_significant": int(top30["survival_sig"].sum()),
            "n_survival_risk_significant": int(top30["survival_risk_sig"].sum()),
            "details": method_b.to_dict(orient="records"),
        },
        "method_c": {
            "name": "Known thyroid drug precision",
            "known_set_size": int(len(known_names)),
            "precision_at_k": p_at_k,
        },
        "top15_validated": top15[
            [c for c in ["final_rank"] + score_cols if c in top15.columns]
        ].to_dict(orient="records"),
    }

    all_30_path = out_dir / "all_30_scores.csv"
    method_a_path = out_dir / "method_a_target_expression_validation.csv"
    method_b_path = out_dir / "method_b_survival_validation.csv"
    method_c_path = out_dir / "method_c_known_thyroid_precision.csv"
    top15_path = out_dir / "top15_validated.csv"
    json_path = out_dir / "step6_external_cohort_results.json"
    report_path = report_dir / f"THYROID_STEP6_EXTERNAL_COHORT_VALIDATION_{STAMP}.md"

    top30[score_cols].sort_values("rank").to_csv(all_30_path, index=False)
    method_a.sort_values("rank").to_csv(method_a_path, index=False)
    method_b.sort_values("rank").to_csv(method_b_path, index=False)
    pd.DataFrame(
        [
            {
                "metric": metric,
                "precision": values["precision"],
                "hits": values["hits"],
                "total": values["total"],
            }
            for metric, values in p_at_k.items()
        ]
    ).to_csv(method_c_path, index=False)
    top15[[c for c in ["final_rank"] + score_cols if c in top15.columns]].to_csv(top15_path, index=False)
    write_json(json_path, summary)

    report_lines = [
        f"# Thyroid Step 6 External Cohort Validation - {date.today().isoformat()}",
        "",
        "This run mirrors the BRCA/PAAD separated Step 6 style for the final thyroid model output.",
        "",
        "## Step 6 counts",
        "",
        f"- Top30 model outputs reviewed: `{len(top30)}`",
        f"- External THCA target expression support: `{int(top30['external_cohort_validated'].sum())}/{len(top30)}`",
        f"- Survival significant (p < 0.05): `{int(top30['survival_sig'].sum())}/{len(top30)}`",
        f"- Known thyroid controls in Top30: `{sum(values['hits'] for values in p_at_k.values() if values['total'] == 30) if 'P@30' in p_at_k else 0}/{p_at_k.get('P@30', {}).get('total', 0)}`",
        "",
        "## Method A. Target Expression",
        "",
        markdown_table(method_a.head(15)),
        "",
        "## Method B. Survival",
        "",
        markdown_table(method_b.head(15)),
        "",
        "## Method C. Known Thyroid Drug Precision",
        "",
        markdown_table(pd.DataFrame(
            [
                {"metric": metric, "precision": values["precision"], "hits": values["hits"], "total": values["total"]}
                for metric, values in p_at_k.items()
            ]
        )),
        "",
        "## Top15 for Step 7",
        "",
        markdown_table(top15[[c for c in ["final_rank", "rank", "drug_name", "validation_score", "target_expressed", "survival_p", "known_thyroid"] if c in top15.columns]]),
        "",
    ]
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    return top15, summary, {
        "external_dir": out_dir,
        "report_dir": report_dir,
        "report_path": report_path,
    }


def load_assays(cfg: dict[str, Any]) -> dict[str, pd.DataFrame]:
    paths = pipeline_paths(cfg)
    assay_cache: dict[str, pd.DataFrame] = {}
    for assay_name in ADMET_ASSAYS:
        path = paths.admet_source_dir / "tdc" / f"{assay_name}.csv"
        if not path.exists():
            path = paths.admet_source_dir / f"{assay_name}.csv"
        if path.exists():
            _, assay_df = _load_assay(path, cfg)
            assay_cache[assay_name] = assay_df
        else:
            assay_cache[assay_name] = pd.DataFrame()
    return assay_cache


def run_brca_style_admet(cfg: dict[str, Any], top15: pd.DataFrame) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    assay_cache = load_assays(cfg)
    profiles: list[dict[str, Any]] = []
    detailed_rows: list[dict[str, Any]] = []

    for _, row in top15.iterrows():
        drug_id = str(row["canonical_drug_id"])
        drug_name = str(row["drug_name"])
        smiles = str(row.get("canonical_smiles", "") or "")
        can, ok = canonicalize_smiles(smiles)
        candidate_fp = _fingerprint(can, cfg) if ok else None

        n_assays_tested = 0
        n_pass = 0
        n_caution = 0
        n_nodata = 0
        flags: list[str] = []
        assay_details: dict[str, Any] = {}

        for assay_name, assay_info in ADMET_ASSAYS.items():
            match = _nearest_admet_match(candidate_fp, assay_cache.get(assay_name, pd.DataFrame()), cfg)
            predicted = match.get("predicted_label")
            match_type = str(match.get("match_type", "no_match")).lower()
            similarity = match.get("similarity", np.nan)

            if predicted is None or (isinstance(predicted, float) and not np.isfinite(predicted)) or match_type in {"no_match", "no_smiles"}:
                status = "no_data"
                n_nodata += 1
            elif assay_info["type"] == "binary":
                n_assays_tested += 1
                good_val = assay_info.get("good_value")
                if good_val is None:
                    status = "info"
                elif int(float(predicted)) == int(good_val):
                    status = "pass"
                    n_pass += 1
                else:
                    if assay_name in {"ames", "dili", "herg"}:
                        status = "caution"
                        n_caution += 1
                        flags.append(f"{assay_info['name']}(+)")
                    else:
                        status = "minor"
                        n_caution += 1
            else:
                n_assays_tested += 1
                status = "measured"
                n_pass += 1

            assay_details[assay_name] = {
                "value": predicted,
                "status": status,
                "match_type": match.get("match_type"),
                "similarity": similarity,
            }
            detailed_rows.append(
                {
                    "canonical_drug_id": drug_id,
                    "drug_name": drug_name,
                    "assay": assay_name,
                    "category": assay_info["category"],
                    "assay_name": assay_info["name"],
                    "assay_type": assay_info["type"],
                    "match_type": match.get("match_type"),
                    "similarity": similarity,
                    "predicted_label": predicted,
                    "status": status,
                }
            )

        if n_assays_tested > 0:
            safety_score = n_pass / max(n_assays_tested, 1) * 10
            for flag in flags:
                if "hERG" in flag:
                    safety_score -= 2
                elif "Ames" in flag:
                    safety_score -= 1.5
                elif "DILI" in flag:
                    safety_score -= 1
        else:
            safety_score = 5.0

        known_approved = drug_name in KNOWN_APPROVED
        if known_approved:
            safety_score += 2.0

        profiles.append(
            {
                "drug_id": drug_id,
                "canonical_drug_id": drug_id,
                "drug_name": drug_name,
                "target": row.get("target_genes", ""),
                "pathway": row.get("PATHWAY_NAME_NORMALIZED", ""),
                "pred_ic50": float(row.get("mean_pred_ln_ic50", np.nan)),
                "validation_score": float(row.get("validation_score", np.nan)),
                "step6_rank": int(row.get("final_rank", np.nan)),
                "model_rank": int(row.get("rank", np.nan)),
                "n_assays_tested": int(n_assays_tested),
                "n_pass": int(n_pass),
                "n_caution": int(n_caution),
                "n_nodata": int(n_nodata),
                "safety_score": round(float(safety_score), 2),
                "flags": flags,
                "known_approved": known_approved,
                "assay_details": assay_details,
            }
        )

    for profile in profiles:
        ic50_rank = sorted(profiles, key=lambda x: x["pred_ic50"]).index(profile)
        profile["efficacy_rank"] = ic50_rank + 1
        profile["combined_score"] = round(profile["safety_score"] + (15 - ic50_rank) * 0.5, 2)

    profiles.sort(key=lambda x: -x["combined_score"])

    final: list[dict[str, Any]] = []
    for idx, profile in enumerate(profiles, start=1):
        if profile["known_approved"]:
            category = "Approved"
        elif profile["safety_score"] >= 4.0:
            category = "Candidate"
        else:
            category = "Caution"
        final.append(
            {
                "final_rank": idx,
                "drug_id": profile["drug_id"],
                "canonical_drug_id": profile["canonical_drug_id"],
                "drug_name": profile["drug_name"],
                "target": profile["target"],
                "pathway": profile["pathway"],
                "pred_ic50": profile["pred_ic50"],
                "validation_score": profile["validation_score"],
                "step6_rank": profile["step6_rank"],
                "model_rank": profile["model_rank"],
                "safety_score": profile["safety_score"],
                "category": category,
                "flags": profile["flags"],
                "combined_score": profile["combined_score"],
                "n_assays_tested": profile["n_assays_tested"],
                "known_approved": profile["known_approved"],
            }
        )

    return profiles, final, detailed_rows


def write_step7_report(path: Path, final_df: pd.DataFrame) -> None:
    lines = [
        f"# Thyroid Step 7 ADMET With BRCA Rules - {date.today().isoformat()}",
        "",
        "This run applies the BRCA Step7 ADMET rules unchanged to the thyroid Step6 top15.",
        "",
        "## Rules copied from BRCA",
        "",
        "- Same 22 ADMET assays",
        "- Same safety score formula",
        "- Same toxicity penalties: hERG -2, Ames -1.5, DILI -1",
        "- Same known-approved override set",
        "- Final category: `Approved` if known-approved, else `Candidate` if safety_score >= 4, else `Caution`",
        "",
        "## Final Candidates",
        "",
        markdown_table(final_df[[c for c in ["final_rank", "drug_name", "step6_rank", "model_rank", "validation_score", "safety_score", "category", "n_assays_tested", "combined_score", "flags"] if c in final_df.columns]]),
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_step7(
    cfg: dict[str, Any],
    top15: pd.DataFrame,
    step7_variant: str,
) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Path]]:
    cfg2 = deepcopy(cfg)
    cfg2["paths"]["admet_output_dir"] = f"admet/{step7_variant}"
    cfg2["paths"]["reports_dir"] = f"reports/{step7_variant}"
    paths2 = pipeline_paths(cfg2)
    paths2.admet_output_dir.mkdir(parents=True, exist_ok=True)
    paths2.reports_dir.mkdir(parents=True, exist_ok=True)

    profiles, final_rows, detailed_rows = run_brca_style_admet(cfg, top15)
    detailed_df = pd.DataFrame(detailed_rows)
    final_df = pd.DataFrame(final_rows)

    detailed_df.to_csv(paths2.admet_output_dir / "admet_detailed_candidates.csv", index=False)
    final_df.to_csv(paths2.admet_output_dir / "final_drug_candidates.csv", index=False)

    result_payload = {
        "status": "completed",
        "category_counts": final_df["category"].value_counts(dropna=False).to_dict(),
        "top10": final_rows[:10],
        "profiles": profiles,
    }
    write_json(paths2.admet_output_dir / "step7_admet_results.json", result_payload)
    write_json(paths2.reports_dir / "qc_step7_admet_brca_rules.json", {"category_counts": result_payload["category_counts"], "n_final": int(len(final_df))})
    report_path = paths2.reports_dir / f"THYROID_STEP7_ADMET_BRCA_RULES_{STAMP}.md"
    write_step7_report(report_path, final_df)

    return final_df, result_payload, {
        "admet_dir": paths2.admet_output_dir,
        "report_dir": paths2.reports_dir,
        "report_path": report_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run thyroid Step6/Step7 separated validation with BRCA Step7 ADMET rules")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--source-variant", default=DEFAULT_SOURCE_VARIANT)
    parser.add_argument("--top30", default=DEFAULT_TOP30)
    parser.add_argument("--step6-variant", default=DEFAULT_STEP6_VARIANT)
    parser.add_argument("--step7-variant", default=DEFAULT_STEP7_VARIANT)
    args = parser.parse_args()

    cfg = load_config(args.config)
    top30_path = Path(args.top30)
    if not top30_path.is_absolute():
        top30_path = ROOT / args.top30

    top15, step6_summary, step6_paths = build_step6(cfg, top30_path, args.source_variant, args.step6_variant)
    final_df, step7_summary, step7_paths = build_step7(cfg, top15, args.step7_variant)

    run_summary = {
        "status": "completed",
        "source_variant": args.source_variant,
        "step6_variant": args.step6_variant,
        "step7_variant": args.step7_variant,
        "step6_counts": {
            "external_target_expression": int(step6_summary["method_a"]["n_external_cohort_validated"]),
            "top30": int(step6_summary["n_total"]),
            "survival_sig": int(step6_summary["method_b"]["n_survival_significant"]),
            "known_p_at_5": step6_summary["method_c"]["precision_at_k"]["P@5"],
            "known_p_at_10": step6_summary["method_c"]["precision_at_k"]["P@10"],
            "known_p_at_20": step6_summary["method_c"]["precision_at_k"]["P@20"],
        },
        "step7_category_counts": step7_summary["category_counts"],
        "outputs": {
            "step6_external_dir": str(step6_paths["external_dir"]),
            "step6_report_dir": str(step6_paths["report_dir"]),
            "step7_admet_dir": str(step7_paths["admet_dir"]),
            "step7_report_dir": str(step7_paths["report_dir"]),
        },
    }
    write_json(ROOT / "reports" / f"thyroid_step6_step7_brca_run_summary_{STAMP}.json", run_summary)
    print(json.dumps(json_safe(run_summary), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
