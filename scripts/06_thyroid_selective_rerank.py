from __future__ import annotations

import argparse
import copy
import json
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import (  # noqa: E402
    NAME_NORMALIZER,
    admet_assessment,
    external_validation,
    knowledge_validation,
    load_config,
    pipeline_paths,
    read_table,
    split_genes,
    write_json,
    write_table,
)


THYROID_CONTEXT_GENES = {
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
    "EGFR",
    "MET",
    "KDR",
    "FLT1",
    "FLT4",
    "FGFR1",
    "FGFR2",
    "FGFR3",
    "FGFR4",
    "PIK3CA",
    "PIK3CB",
    "PIK3CD",
    "PIK3CG",
    "AKT1",
    "AKT2",
    "MTOR",
    "CDK4",
    "CDK6",
    "PTEN",
    "TP53",
    "TERT",
}

THYROID_CONTEXT_TERMS = {
    "BRAF",
    "RET",
    "NTRK",
    "RAS",
    "MAPK",
    "ERK",
    "RTK",
    "EGFR",
    "VEGFR",
    "FGFR",
    "MET",
    "PI3K",
    "AKT",
    "MTOR",
    "KINASE",
    "ANGIOGENESIS",
}

BROAD_CYTOTOXIC_TERMS = {
    "MITOSIS",
    "MICROTUBULE",
    "DNA REPLICATION",
    "DNA DAMAGE",
    "DNA CROSSLINKER",
    "TOPOISOMERASE",
    "PROTEIN STABILITY",
    "PROTEASOME",
    "APOPTOSIS",
    "CYTOSKELETON",
    "CHROMATIN",
}


def _norm_name(value: Any) -> str:
    return NAME_NORMALIZER.sub("", "" if pd.isna(value) else str(value).lower())


def _percentile_high_good(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().sum() == 0:
        return pd.Series(0.5, index=values.index, dtype=float)
    filled = numeric.fillna(numeric.median())
    return filled.rank(method="average", pct=True).astype(float)


def _safe_minmax(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().sum() == 0:
        return pd.Series(0.0, index=values.index, dtype=float)
    filled = numeric.fillna(numeric.median())
    lo = float(filled.min())
    hi = float(filled.max())
    if np.isclose(lo, hi):
        return pd.Series(0.5, index=values.index, dtype=float)
    return ((filled - lo) / (hi - lo)).astype(float)


def _contains_any(text: Any, terms: set[str]) -> bool:
    upper = "" if pd.isna(text) else str(text).upper()
    return any(term in upper for term in terms)


def _context_score(row: pd.Series, known_upper: set[str]) -> float:
    genes = set(split_genes(row.get("target_genes", "")))
    pathway = str(row.get("PATHWAY_NAME_NORMALIZED", "")).upper()
    classification = str(row.get("classification", "")).lower()
    drug_name = str(row.get("drug_name", "")).upper()

    score = 0.0
    if genes & THYROID_CONTEXT_GENES:
        score += 0.55
    if _contains_any(pathway, THYROID_CONTEXT_TERMS):
        score += 0.30
    if classification in {"approved", "indication_expansion"}:
        score += 0.10
    if drug_name in known_upper:
        score += 0.15
    return float(min(score, 1.0))


def _broad_cytotoxic_flag(row: pd.Series) -> bool:
    text = f"{row.get('target_genes', '')} {row.get('PATHWAY_NAME_NORMALIZED', '')}"
    return _contains_any(text, BROAD_CYTOTOXIC_TERMS)


def _model_drug_scores(cfg: dict[str, Any], input_name: str) -> pd.DataFrame:
    paths = pipeline_paths(cfg)
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    if rows is None:
        raise FileNotFoundError("data/row_metadata.csv is missing")
    pred_path = paths.results_dir / "pancancer_lincs" / "ensemble" / f"{input_name}_weighted_ensemble_oof.npy"
    if not pred_path.exists():
        raise FileNotFoundError(pred_path)

    pred = np.load(pred_path)
    if len(pred) != len(rows):
        raise ValueError(f"Prediction length {len(pred)} does not match row_metadata rows {len(rows)}")

    candidate_rows = rows.copy()
    candidate_rows["ensemble_pred_ln_ic50"] = pred.astype(float)
    candidate_rows["ensemble_score"] = -candidate_rows["ensemble_pred_ln_ic50"]
    group_cols = [
        "canonical_drug_id",
        "drug_name",
        "canonical_smiles",
        "target_genes",
        "PATHWAY_NAME_NORMALIZED",
        "classification",
    ]
    available = [c for c in group_cols if c in candidate_rows.columns]
    out = (
        candidate_rows.groupby(available, dropna=False)
        .agg(
            mean_pred_ln_ic50=("ensemble_pred_ln_ic50", "mean"),
            ensemble_score=("ensemble_score", "mean"),
            pred_sd_ln_ic50=("ensemble_pred_ln_ic50", "std"),
            screened_rows=("ensemble_score", "size"),
            thyroid_cell_lines=("sample_id", "nunique"),
        )
        .reset_index()
    )
    out["canonical_drug_id"] = out["canonical_drug_id"].astype(str)
    out["_drug_name_norm"] = out["drug_name"].map(_norm_name)
    out = out.sort_values("ensemble_score", ascending=False).drop_duplicates("_drug_name_norm", keep="first")
    return out.drop(columns=["_drug_name_norm"])


def _thyroid_response_summary(cfg: dict[str, Any]) -> pd.DataFrame:
    paths = pipeline_paths(cfg)
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    if rows is None:
        raise FileNotFoundError("data/row_metadata.csv is missing")
    rows = rows.copy()
    rows["canonical_drug_id"] = rows["canonical_drug_id"].astype(str)
    global_sensitive_cut = float(pd.to_numeric(rows["LN_IC50"], errors="coerce").median())
    summary = (
        rows.groupby("canonical_drug_id", dropna=False)
        .agg(
            thyroid_mean_ln_ic50=("LN_IC50", "mean"),
            thyroid_median_ln_ic50=("LN_IC50", "median"),
            thyroid_q25_ln_ic50=("LN_IC50", lambda s: float(pd.to_numeric(s, errors="coerce").quantile(0.25))),
            thyroid_observed_sd=("LN_IC50", "std"),
            thyroid_observed_rows=("LN_IC50", "size"),
            thyroid_observed_cell_lines=("sample_id", "nunique"),
            thyroid_sensitive_rate=("LN_IC50", lambda s: float((pd.to_numeric(s, errors="coerce") <= global_sensitive_cut).mean())),
        )
        .reset_index()
    )
    summary["thyroid_observed_score"] = -summary["thyroid_mean_ln_ic50"]
    return summary


def _pancancer_response_summary(cfg: dict[str, Any]) -> pd.DataFrame:
    paths = pipeline_paths(cfg)
    source = paths.processed_dir / "source_staging" / "gdsc" / "gdsc2_annotation_normalized_20260406.parquet"
    if not source.exists():
        raise FileNotFoundError(source)
    gdsc = pd.read_parquet(
        source,
        columns=["DRUG_ID", "DRUG_NAME", "TCGA_DESC", "SANGER_MODEL_ID", "LN_IC50"],
    )
    gdsc = gdsc.copy()
    gdsc["canonical_drug_id"] = gdsc["DRUG_ID"].astype(str)
    gdsc["LN_IC50"] = pd.to_numeric(gdsc["LN_IC50"], errors="coerce")
    gdsc = gdsc.loc[gdsc["LN_IC50"].notna()].copy()
    non_thyroid = gdsc.loc[~gdsc["TCGA_DESC"].astype(str).str.upper().eq("THCA")].copy()
    summary = (
        non_thyroid.groupby("canonical_drug_id", dropna=False)
        .agg(
            pancancer_mean_ln_ic50=("LN_IC50", "mean"),
            pancancer_median_ln_ic50=("LN_IC50", "median"),
            pancancer_q25_ln_ic50=("LN_IC50", lambda s: float(pd.to_numeric(s, errors="coerce").quantile(0.25))),
            pancancer_rows=("LN_IC50", "size"),
            pancancer_cell_lines=("SANGER_MODEL_ID", "nunique"),
            pancancer_tcga_types=("TCGA_DESC", "nunique"),
        )
        .reset_index()
    )
    summary["canonical_drug_id"] = summary["canonical_drug_id"].astype(str)
    summary["pancancer_activity_score"] = -summary["pancancer_mean_ln_ic50"]
    return summary


def _compare_brca_overlap(cfg: dict[str, Any], ranked: pd.DataFrame, brca_top30_path: str | None) -> dict[str, Any]:
    if brca_top30_path is None:
        default = Path("/Users/skku_aws2_18/team4_project/pre_project/say2_preproject/models/ensemble_results_random3_strong_context_smiles/top30_drugs.csv")
        brca_path = default if default.exists() else None
    else:
        brca_path = Path(brca_top30_path)
    if brca_path is None or not brca_path.exists():
        return {"status": "skipped", "reason": "brca_top30_file_missing"}

    brca = pd.read_csv(brca_path)
    drug_col = "drug_name" if "drug_name" in brca.columns else brca.columns[0]
    brca_names = brca[drug_col].astype(str).tolist()
    thyroid_names = ranked["drug_name"].astype(str).tolist()

    def key(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", value.lower())

    out: dict[str, Any] = {"status": "completed", "brca_top30_path": str(brca_path)}
    for k in [10, 15, 20, 30]:
        b = {key(x): x for x in brca_names[:k]}
        t = {key(x): x for x in thyroid_names[:k]}
        common = [t[n] for n in t if n in b]
        out[f"top{k}_overlap_count"] = int(len(common))
        out[f"top{k}_overlap_ratio"] = float(len(common) / min(k, len(thyroid_names), len(brca_names)))
        out[f"top{k}_overlap_drugs"] = common
    return out


def build_thyroid_selective_ranking(
    cfg: dict[str, Any],
    input_name: str,
    brca_top30_path: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    paths = pipeline_paths(cfg)
    out_dir = paths.results_dir / "thyroid_selective_rerank"
    out_dir.mkdir(parents=True, exist_ok=True)

    known_upper = {str(x).upper() for x in cfg.get("known_thyroid_drugs", [])}
    model_scores = _model_drug_scores(cfg, input_name)
    thyroid = _thyroid_response_summary(cfg)
    pan = _pancancer_response_summary(cfg)
    ranked = model_scores.merge(thyroid, on="canonical_drug_id", how="left").merge(pan, on="canonical_drug_id", how="left")

    ranked["selectivity_delta_ln_ic50"] = ranked["pancancer_mean_ln_ic50"] - ranked["thyroid_mean_ln_ic50"]
    ranked["model_rank_score"] = _percentile_high_good(ranked["ensemble_score"])
    ranked["thyroid_response_rank_score"] = _percentile_high_good(ranked["thyroid_observed_score"])
    ranked["selectivity_rank_score"] = _percentile_high_good(ranked["selectivity_delta_ln_ic50"])
    ranked["pancancer_activity_rank_score"] = _percentile_high_good(ranked["pancancer_activity_score"])
    ranked["thyroid_context_score"] = ranked.apply(lambda row: _context_score(row, known_upper), axis=1)
    ranked["broad_cytotoxic_flag"] = ranked.apply(_broad_cytotoxic_flag, axis=1)
    ranked["broad_cytotoxic_penalty"] = 0.0
    ranked.loc[ranked["broad_cytotoxic_flag"], "broad_cytotoxic_penalty"] += 0.20
    ranked.loc[
        (ranked["pancancer_activity_rank_score"] >= 0.75) & (ranked["selectivity_delta_ln_ic50"].fillna(0) <= 0),
        "broad_cytotoxic_penalty",
    ] += 0.20

    ranked["thyroid_selective_score"] = (
        0.30 * ranked["model_rank_score"]
        + 0.30 * ranked["selectivity_rank_score"]
        + 0.20 * ranked["thyroid_response_rank_score"]
        + 0.20 * ranked["thyroid_context_score"]
        - ranked["broad_cytotoxic_penalty"]
    )
    ranked["thyroid_selective_score"] = ranked["thyroid_selective_score"].clip(lower=0.0)
    ranked["ranking_layer"] = np.select(
        [
            ranked["thyroid_context_score"].ge(0.45) & ranked["selectivity_delta_ln_ic50"].gt(0),
            ranked["broad_cytotoxic_flag"],
        ],
        ["thyroid_context_supported", "broad_cytotoxic_downweighted"],
        default="screened_selectivity_candidate",
    )

    output_cols = [
        "drug_name",
        "canonical_drug_id",
        "canonical_smiles",
        "target_genes",
        "PATHWAY_NAME_NORMALIZED",
        "classification",
        "ranking_layer",
        "thyroid_selective_score",
        "ensemble_score",
        "mean_pred_ln_ic50",
        "thyroid_mean_ln_ic50",
        "pancancer_mean_ln_ic50",
        "selectivity_delta_ln_ic50",
        "thyroid_sensitive_rate",
        "thyroid_context_score",
        "broad_cytotoxic_flag",
        "broad_cytotoxic_penalty",
        "thyroid_observed_cell_lines",
        "pancancer_cell_lines",
        "pancancer_tcga_types",
    ]
    ranked = ranked.sort_values("thyroid_selective_score", ascending=False).reset_index(drop=True)
    ranked["thyroid_selective_rank"] = np.arange(1, len(ranked) + 1)
    ranked = ranked[["thyroid_selective_rank"] + [c for c in output_cols if c in ranked.columns]]

    targeted = ranked.loc[
        ranked["ranking_layer"].isin(["thyroid_context_supported", "screened_selectivity_candidate"])
        & ~ranked["broad_cytotoxic_flag"].astype(bool)
    ].copy()
    targeted = targeted.sort_values("thyroid_selective_score", ascending=False).reset_index(drop=True)
    targeted["thyroid_targeted_rank"] = np.arange(1, len(targeted) + 1)

    write_table(ranked, out_dir / "thyroid_selective_reranked_all_drugs.csv")
    write_table(ranked.head(30), out_dir / "thyroid_selective_reranked_top30.csv")
    write_table(targeted, out_dir / "thyroid_targeted_non_broad_reranked_all_drugs.csv")
    write_table(targeted.head(30), out_dir / "thyroid_targeted_non_broad_reranked_top30.csv")

    overlap_all = _compare_brca_overlap(cfg, ranked, brca_top30_path)
    overlap_targeted = _compare_brca_overlap(cfg, targeted, brca_top30_path)
    qc = {
        "step": "thyroid_selective_rerank",
        "input_set": input_name,
        "drug_count": int(len(ranked)),
        "targeted_non_broad_count": int(len(targeted)),
        "score_formula": {
            "model_rank_score": 0.30,
            "selectivity_rank_score": 0.30,
            "thyroid_response_rank_score": 0.20,
            "thyroid_context_score": 0.20,
            "broad_cytotoxic_penalty": "subtract 0.20 for broad target/pathway and 0.20 for broad pan-cancer activity without thyroid selectivity",
        },
        "layer_counts": ranked["ranking_layer"].value_counts(dropna=False).to_dict(),
        "top30_broad_cytotoxic_count": int(ranked.head(30)["broad_cytotoxic_flag"].sum()),
        "targeted_top30_broad_cytotoxic_count": int(targeted.head(30)["broad_cytotoxic_flag"].sum()) if not targeted.empty else 0,
        "brca_overlap_all_rerank": overlap_all,
        "brca_overlap_targeted_non_broad": overlap_targeted,
        "outputs": {
            "all_ranked": str(out_dir / "thyroid_selective_reranked_all_drugs.csv"),
            "top30": str(out_dir / "thyroid_selective_reranked_top30.csv"),
            "targeted_all_ranked": str(out_dir / "thyroid_targeted_non_broad_reranked_all_drugs.csv"),
            "targeted_top30": str(out_dir / "thyroid_targeted_non_broad_reranked_top30.csv"),
        },
    }
    write_json(paths.reports_dir / "qc_thyroid_selective_rerank_20260421.json", qc)
    return ranked, targeted, qc


def run_downstream_validation(cfg: dict[str, Any], candidates: pd.DataFrame) -> dict[str, Any]:
    cfg2 = copy.deepcopy(cfg)
    cfg2["paths"] = copy.deepcopy(cfg["paths"])
    cfg2["paths"]["external_validation_dir"] = "external_validation/thyroid_selective"
    cfg2["paths"]["admet_output_dir"] = "admet/thyroid_selective"
    cfg2["paths"]["knowledge_validation_dir"] = "knowledge_validation/thyroid_selective"
    cfg2["paths"]["phase5_dir"] = "phase5_final_results/thyroid_selective"
    paths2 = pipeline_paths(cfg2)

    top_for_validation = candidates.head(30).copy()
    if "rank" not in top_for_validation.columns:
        top_for_validation = top_for_validation.rename(columns={"thyroid_selective_rank": "rank", "thyroid_targeted_rank": "rank"})
    if "ensemble_score" not in top_for_validation.columns and "thyroid_selective_score" in top_for_validation.columns:
        top_for_validation["ensemble_score"] = top_for_validation["thyroid_selective_score"]
    if "mean_pred_ln_ic50" not in top_for_validation.columns:
        top_for_validation["mean_pred_ln_ic50"] = -top_for_validation["ensemble_score"]

    validated = external_validation(cfg2, top_for_validation)
    if "canonical_drug_id" in validated.columns:
        validated = validated.copy()
        validated["canonical_drug_id"] = validated["canonical_drug_id"].astype(str)
    admet = admet_assessment(cfg2, validated)
    final = knowledge_validation(cfg2, admet)
    summary = {
        "external_validation_rows": int(len(validated)),
        "admet_rows": int(len(admet)),
        "knowledge_rows": int(len(final)),
        "tier_counts": final["tier"].value_counts(dropna=False).to_dict() if not final.empty else {},
        "outputs": {
            "external_top15": str(paths2.external_validation_dir / "top15_validated.csv"),
            "admet_final": str(paths2.admet_output_dir / "final_drug_candidates.csv"),
            "knowledge_final": str(paths2.phase5_dir / "final_comprehensive_candidates.csv"),
        },
    }
    write_json(paths2.reports_dir / "qc_thyroid_selective_downstream_20260421.json", summary)
    return summary


def _fmt(value: Any) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    try:
        return f"{float(value):.4f}"
    except Exception:
        return str(value)


def _md_table(df: pd.DataFrame, cols: list[str]) -> str:
    if df.empty:
        return "_No rows._"
    subset = df[[c for c in cols if c in df.columns]].copy()
    lines = ["| " + " | ".join(subset.columns) + " |", "| " + " | ".join(["---"] * len(subset.columns)) + " |"]
    for row in subset.itertuples(index=False):
        lines.append("| " + " | ".join(_fmt(v) for v in row) + " |")
    return "\n".join(lines)


def write_report(cfg: dict[str, Any], ranked: pd.DataFrame, targeted: pd.DataFrame, qc: dict[str, Any], downstream: dict[str, Any]) -> Path:
    paths = pipeline_paths(cfg)
    report_path = paths.root / "docs" / "THYROID_SELECTIVE_RERANK_20260421.md"
    original_overlap = {
        "top10": "9/10",
        "top15": "14/15",
        "top30": "25/30",
    }
    cols = [
        "thyroid_selective_rank",
        "thyroid_targeted_rank",
        "drug_name",
        "thyroid_selective_score",
        "ranking_layer",
        "ensemble_score",
        "thyroid_mean_ln_ic50",
        "pancancer_mean_ln_ic50",
        "selectivity_delta_ln_ic50",
        "thyroid_context_score",
        "broad_cytotoxic_flag",
    ]
    text = f"""# Thyroid-Selective Re-ranking - 2026-04-21

## 목적

pan-cancer LINCS 모델의 상위 후보가 BRCA 추천 후보와 과도하게 겹치는 문제를 줄이기 위해, 모델 예측 점수만이 아니라 thyroid cell line에서의 상대 선택성, thyroid 관련 target/pathway context, broad cytotoxic penalty를 함께 반영한 재랭킹 레이어를 추가했다.

## 왜 필요한가

- 기존 pan-cancer LINCS Top 후보는 BRCA random3 strong-context-smiles Top 후보와 많이 겹쳤다.
- 기존 겹침: Top10 {original_overlap["top10"]}, Top15 {original_overlap["top15"]}, Top30 {original_overlap["top30"]}.
- 이는 오류라기보다 GDSC screened pool에서 전반적으로 강한 항암제와 broad cytotoxic 후보가 공통으로 상위에 올라오는 구조 때문이다.

## Score 정의

`thyroid_selective_score = 0.30 * model_rank_score + 0.30 * selectivity_rank_score + 0.20 * thyroid_response_rank_score + 0.20 * thyroid_context_score - broad_cytotoxic_penalty`

- `model_rank_score`: pan-cancer LINCS ensemble 예측 민감도 rank.
- `selectivity_rank_score`: `pancancer_mean_ln_ic50 - thyroid_mean_ln_ic50` rank. 값이 클수록 thyroid가 pan-cancer 평균보다 더 민감하다.
- `thyroid_response_rank_score`: thyroid screened response의 observed sensitivity rank.
- `thyroid_context_score`: BRAF/RET/NTRK/RAS/MAPK/VEGFR/PI3K/MTOR/CDK 등 thyroid 관련 target/pathway 점수.
- `broad_cytotoxic_penalty`: mitosis, microtubule, DNA replication/damage, proteasome 등 broad cytotoxic mechanism 및 pan-cancer broad activity에 대한 감점.

## 전체 재랭킹 Top 20

{_md_table(ranked.head(20), cols)}

## Broad Cytotoxic 제외 Targeted/Selective Top 20

{_md_table(targeted.head(20), cols)}

## BRCA 후보와 겹침 변화

### 전체 재랭킹

- Top10 overlap: {qc["brca_overlap_all_rerank"].get("top10_overlap_count")}/10
- Top15 overlap: {qc["brca_overlap_all_rerank"].get("top15_overlap_count")}/15
- Top30 overlap: {qc["brca_overlap_all_rerank"].get("top30_overlap_count")}/30

### Broad Cytotoxic 제외 Targeted/Selective 리스트

- Top10 overlap: {qc["brca_overlap_targeted_non_broad"].get("top10_overlap_count")}/10
- Top15 overlap: {qc["brca_overlap_targeted_non_broad"].get("top15_overlap_count")}/15
- Top30 overlap: {qc["brca_overlap_targeted_non_broad"].get("top30_overlap_count")}/30

## Downstream validation

- External validation rows: {downstream.get("external_validation_rows")}
- ADMET rows: {downstream.get("admet_rows")}
- Knowledge validation rows: {downstream.get("knowledge_rows")}
- Tier counts: `{json.dumps(downstream.get("tier_counts", {}), ensure_ascii=False)}`

## 산출물

- `results/thyroid_selective_rerank/thyroid_selective_reranked_top30.csv`
- `results/thyroid_selective_rerank/thyroid_targeted_non_broad_reranked_top30.csv`
- `external_validation/thyroid_selective/top15_validated.csv`
- `admet/thyroid_selective/final_drug_candidates.csv`
- `phase5_final_results/thyroid_selective/final_comprehensive_candidates.csv`
- `reports/qc_thyroid_selective_rerank_20260421.json`
"""
    report_path.write_text(text, encoding="utf-8")
    return report_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Re-rank thyroid candidates with thyroid-selective evidence")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--input-set", default="numeric_strong_context_smiles_pan_lincs")
    parser.add_argument("--brca-top30", default=None)
    parser.add_argument("--skip-downstream", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    ranked, targeted, qc = build_thyroid_selective_ranking(cfg, args.input_set, args.brca_top30)
    downstream = {}
    if not args.skip_downstream:
        downstream = run_downstream_validation(cfg, targeted if not targeted.empty else ranked)
    report = write_report(cfg, ranked, targeted, qc, downstream)
    print(
        json.dumps(
            {
                "status": "completed",
                "top_drug": ranked.iloc[0]["drug_name"] if not ranked.empty else None,
                "targeted_top_drug": targeted.iloc[0]["drug_name"] if not targeted.empty else None,
                "report": str(report),
                "qc": qc,
                "downstream": downstream,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
