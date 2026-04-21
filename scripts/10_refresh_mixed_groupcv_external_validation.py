from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import (  # noqa: E402
    admet_assessment,
    external_validation,
    knowledge_validation,
    load_config,
    pipeline_paths,
    write_json,
)


DEFAULT_CANDIDATES = "results/mixed_groupcv/crossattention_residualmlp_lightgbm_groupcv/mixed_groupcv_ensemble_top30_drugs.csv"
VARIANT_NAME = "mixed_groupcv_pan_lincs"


def _fmt(value: Any) -> str:
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.4f}"
    except Exception:
        return str(value)


def _md_table(df: pd.DataFrame, cols: list[str], n: int = 15) -> str:
    sub = df[[c for c in cols if c in df.columns]].head(n).copy()
    if sub.empty:
        return "_No rows._"
    lines = ["| " + " | ".join(sub.columns) + " |", "| " + " | ".join(["---"] * len(sub.columns)) + " |"]
    for row in sub.itertuples(index=False):
        lines.append("| " + " | ".join(_fmt(v) for v in row) + " |")
    return "\n".join(lines)


def _prepare_candidates(path: Path) -> pd.DataFrame:
    candidates = pd.read_csv(path)
    if "mixed_groupcv_score" in candidates.columns and "ensemble_score" not in candidates.columns:
        candidates = candidates.rename(columns={"mixed_groupcv_score": "ensemble_score"})
    if "mixed_score" in candidates.columns and "ensemble_score" not in candidates.columns:
        candidates = candidates.rename(columns={"mixed_score": "ensemble_score"})
    required = {"rank", "canonical_drug_id", "drug_name", "canonical_smiles", "target_genes", "ensemble_score", "mean_pred_ln_ic50"}
    missing = sorted(required - set(candidates.columns))
    if missing:
        raise ValueError(f"Candidate file is missing required columns: {missing}")
    candidates["rank"] = candidates["rank"].astype(int)
    candidates = candidates.sort_values("rank").reset_index(drop=True)
    return candidates


def _variant_cfg(cfg: dict[str, Any], variant: str) -> dict[str, Any]:
    cfg = json.loads(json.dumps(cfg))
    cfg["paths"]["external_validation_dir"] = f"external_validation/{variant}"
    cfg["paths"]["admet_output_dir"] = f"admet/{variant}"
    cfg["paths"]["knowledge_validation_dir"] = f"knowledge_validation/{variant}"
    cfg["paths"]["phase5_dir"] = f"phase5_final_results/{variant}"
    cfg["paths"]["reports_dir"] = f"reports/{variant}"
    return cfg


def write_report(
    cfg: dict[str, Any],
    variant: str,
    source_path: Path,
    candidates: pd.DataFrame,
    external: pd.DataFrame,
    admet: pd.DataFrame,
    knowledge: pd.DataFrame,
) -> Path:
    paths = pipeline_paths(cfg)
    report = paths.root / "docs" / "THYROID_MIXED_GROUPCV_EXTERNAL_VALIDATION_20260421.md"
    ext_qc_path = paths.reports_dir / "qc_step8_external_validation.json"
    admet_qc_path = paths.reports_dir / "qc_step9_admet.json"
    kg_qc_path = paths.reports_dir / "qc_step10_knowledge_validation.json"
    ext_qc = json.loads(ext_qc_path.read_text(encoding="utf-8")) if ext_qc_path.exists() else {}
    admet_qc = json.loads(admet_qc_path.read_text(encoding="utf-8")) if admet_qc_path.exists() else {}
    kg_qc = json.loads(kg_qc_path.read_text(encoding="utf-8")) if kg_qc_path.exists() else {}

    candidate_cols = ["rank", "drug_name", "ensemble_score", "mean_pred_ln_ic50", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"]
    external_cols = ["rank", "drug_name", "target_match_genes", "target_expression_pct", "target_expressed", "survival_p_value", "survival_direction", "known_thyroid_control"]
    admet_cols = ["drug_name", "target_expressed", "admet_coverage", "toxicity_flags", "low_confidence_toxic_signals", "admet_category"]
    knowledge_cols = ["drug_name", "ensemble_score", "admet_category", "knowledge_score", "tier", "final_category"]

    text = f"""# Thyroid Mixed GroupCV External Validation - 2026-04-21

## 목적

최신 모델 후보인 `pan-cancer LINCS + CrossAttention + ResidualMLP + LightGBM` GroupCV 앙상블 Top30을 기준으로 외부검증, ADMET, KG/Tier 검증을 다시 수행했다.

## 입력 후보

- Variant: `{variant}`
- Source: `{source_path}`
- Model candidate source: `CrossAttention + ResidualMLP + LightGBM`
- Input features: `numeric + strong context + smiles + pan-cancer LINCS`
- CV basis: `GroupCV by canonical_drug_id`

## 최신 앙상블 Top15

{_md_table(candidates, candidate_cols)}

## Step 8. TCGA-THCA 외부검증 요약

- Expression genes: `{ext_qc.get("expression_gene_count", "")}`
- Patient count: `{ext_qc.get("patient_count", "")}`
- Mean target gene match rate: `{_fmt(ext_qc.get("target_gene_match_rate_mean"))}`
- Target expressed in Top15: `{ext_qc.get("target_expressed_count_top15", "")} / 15`
- Known thyroid positive-control precision: `{json.dumps(ext_qc.get("topk_precision", {}), ensure_ascii=False)}`

{_md_table(external, external_cols)}

## Step 9. ADMET 요약

- Candidate count: `{admet_qc.get("candidate_count", "")}`
- Assay count: `{admet_qc.get("assay_count", "")}`
- Category counts: `{json.dumps(admet_qc.get("category_counts", {}), ensure_ascii=False)}`
- Toxicity flag count: `{admet_qc.get("toxicity_flag_count", "")}`
- Low-confidence toxic signal count: `{admet_qc.get("low_confidence_toxic_signal_count", "")}`

{_md_table(admet, admet_cols)}

## Step 10. KG/Tier 검증 요약

- Candidate count: `{kg_qc.get("candidate_count", "")}`
- Tier counts: `{json.dumps(kg_qc.get("tier_counts", {}), ensure_ascii=False)}`
- Category counts: `{json.dumps(kg_qc.get("category_counts", {}), ensure_ascii=False)}`

{_md_table(knowledge, knowledge_cols)}

## 산출물

- `external_validation/{variant}/top15_validated.csv`
- `external_validation/{variant}/thyroid_target_expression.csv`
- `external_validation/{variant}/thyroid_survival_validation.csv`
- `admet/{variant}/final_drug_candidates.csv`
- `knowledge_validation/{variant}/validation_summary.csv`
- `phase5_final_results/{variant}/final_comprehensive_candidates.csv`
- `phase5_final_results/{variant}/tier1_high_confidence.csv`
- `reports/{variant}/qc_step8_external_validation.json`
- `reports/{variant}/qc_step9_admet.json`
- `reports/{variant}/qc_step10_knowledge_validation.json`
"""
    report.write_text(text, encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Refresh external/ADMET/KG validation for mixed GroupCV thyroid ensemble")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--candidates", default=DEFAULT_CANDIDATES)
    parser.add_argument("--variant", default=VARIANT_NAME)
    args = parser.parse_args()

    base_cfg = load_config(args.config)
    cfg = _variant_cfg(base_cfg, args.variant)
    paths = pipeline_paths(cfg)
    source_path = (paths.root / args.candidates).resolve() if not Path(args.candidates).is_absolute() else Path(args.candidates)
    candidates = _prepare_candidates(source_path)

    external = external_validation(cfg, candidates)
    admet = admet_assessment(cfg, external)
    knowledge = knowledge_validation(cfg, admet)
    report = write_report(cfg, args.variant, source_path, candidates, external, admet, knowledge)

    result = {
        "status": "completed",
        "variant": args.variant,
        "candidate_source": str(source_path),
        "report": str(report),
        "outputs": {
            "external_validation": str(paths.external_validation_dir),
            "admet": str(paths.admet_output_dir),
            "knowledge_validation": str(paths.knowledge_validation_dir),
            "phase5": str(paths.phase5_dir),
            "reports": str(paths.reports_dir),
        },
        "top_tier_rows": knowledge.head(5).to_dict(orient="records"),
    }
    write_json(paths.reports_dir / "mixed_groupcv_external_validation_summary.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
