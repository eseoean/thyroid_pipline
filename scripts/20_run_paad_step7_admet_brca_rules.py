#!/usr/bin/env python3
"""Run PAAD Step 7 ADMET with the exact BRCA Step 7 scoring/category rules."""

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
    _load_assay,
    _nearest_admet_match,
    _fingerprint,
    canonicalize_smiles,
    load_config,
    pipeline_paths,
)


DEFAULT_STEP6_VARIANT = "groupcv4_drug_step6_external_cohort"
DEFAULT_STEP7_VARIANT = "groupcv4_drug_step7_admet_brca_rules"
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


def load_step6_top15(cfg: dict[str, Any], step6_variant: str) -> pd.DataFrame:
    paths = pipeline_paths(cfg)
    path = paths.external_validation_dir / step6_variant / "top15_validated.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing Step6 top15 file: {path}")
    df = pd.read_csv(path)
    df["canonical_drug_id"] = df["canonical_drug_id"].astype(str)
    if "drug_id" not in df.columns:
        df["drug_id"] = df["canonical_drug_id"]
    return df


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
                "step6_rank": int(row.get("rank", np.nan)),
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
    for i, profile in enumerate(profiles, start=1):
        category = "Approved" if profile["known_approved"] else ("Candidate" if profile["safety_score"] >= 4 else "Caution")
        final.append(
            {
                "final_rank": i,
                "drug_id": profile["drug_id"],
                "canonical_drug_id": profile["canonical_drug_id"],
                "drug_name": profile["drug_name"],
                "target": profile["target"],
                "pathway": profile["pathway"],
                "pred_ic50": profile["pred_ic50"],
                "validation_score": profile["validation_score"],
                "step6_rank": profile["step6_rank"],
                "safety_score": profile["safety_score"],
                "category": category,
                "flags": profile["flags"],
                "combined_score": profile["combined_score"],
                "n_assays_tested": profile["n_assays_tested"],
                "known_approved": profile["known_approved"],
            }
        )
    return profiles, final, detailed_rows


def write_report(path: Path, final_df: pd.DataFrame) -> None:
    lines = [
        f"# PAAD Step 7 ADMET With BRCA Rules - {date.today().isoformat()}",
        "",
        "This run applies the BRCA Step7 ADMET rules unchanged to the PAAD Step6 top15.",
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
    ]
    display = [
        "final_rank",
        "drug_name",
        "step6_rank",
        "validation_score",
        "safety_score",
        "category",
        "n_assays_tested",
        "combined_score",
        "flags",
    ]
    lines.append(markdown_table(final_df[[c for c in display if c in final_df.columns]]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run PAAD Step7 ADMET using BRCA Step7 rules")
    parser.add_argument("--config", default="config/paad_pipeline_config.json")
    parser.add_argument("--step6-variant", default=DEFAULT_STEP6_VARIANT)
    parser.add_argument("--step7-variant", default=DEFAULT_STEP7_VARIANT)
    args = parser.parse_args()

    cfg = load_config(args.config)
    top15 = load_step6_top15(cfg, args.step6_variant)

    profiles, final, detailed_rows = run_brca_style_admet(cfg, top15)

    cfg2 = deepcopy(cfg)
    cfg2["paths"]["admet_output_dir"] = f"admet/paad/{args.step7_variant}"
    cfg2["paths"]["reports_dir"] = f"reports/paad/{args.step7_variant}"
    paths2 = pipeline_paths(cfg2)

    final_df = pd.DataFrame(final)
    detailed_df = pd.DataFrame(detailed_rows)

    summary = {
        "step": 7,
        "description": "PAAD Step7 ADMET with BRCA Step7 rules",
        "source_step6_variant": args.step6_variant,
        "n_assays": len(ADMET_ASSAYS),
        "n_drugs_input": len(profiles),
        "n_drugs_output": len(final),
        "assay_list": {k: v["name"] for k, v in ADMET_ASSAYS.items()},
        "drug_profiles": [{k: v for k, v in p.items() if k != "assay_details"} for p in profiles],
        "final_candidates": final,
        "detailed_profiles": profiles,
        "category_counts": final_df["category"].value_counts(dropna=False).to_dict() if not final_df.empty else {},
    }

    final_df.to_csv(paths2.admet_output_dir / "final_drug_candidates.csv", index=False)
    detailed_df.to_csv(paths2.admet_output_dir / "admet_detailed_candidates.csv", index=False)
    write_json(paths2.admet_output_dir / "step7_admet_results.json", summary)
    write_report(paths2.reports_dir / f"PAAD_STEP7_ADMET_BRCA_RULES_{STAMP}.md", final_df)
    write_json(paths2.reports_dir / "qc_step7_admet_brca_rules.json", {"category_counts": summary["category_counts"], "n_assays": len(ADMET_ASSAYS), "n_drugs_input": len(profiles)})

    run_summary = {
        "status": "completed",
        "step6_variant": args.step6_variant,
        "step7_variant": args.step7_variant,
        "category_counts": summary["category_counts"],
        "top10": final_df[["final_rank", "drug_name", "category", "safety_score", "combined_score"]].head(10).to_dict(orient="records"),
        "outputs": {
            "admet_dir": str(paths2.admet_output_dir),
            "report_dir": str(paths2.reports_dir),
        },
    }
    write_json(paths2.root / "reports" / "paad" / f"paad_step7_admet_brca_rules_summary_{STAMP}.json", run_summary)
    print(json.dumps(json_safe(run_summary), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
