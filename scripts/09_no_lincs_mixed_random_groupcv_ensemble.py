from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import (  # noqa: E402
    NAME_NORMALIZER,
    load_config,
    make_model,
    pipeline_paths,
    regression_metrics,
    safe_pearson,
    safe_spearman,
    write_json,
    write_table,
)


INPUT_NAME = "numeric_strong_context_smiles_no_lincs"
ENSEMBLE_NAME = "crossattention_residualmlp_lightgbm_no_lincs"
MEMBERS = ("CrossAttention", "ResidualMLP", "LightGBM")


def _read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _fmt(value: Any) -> str:
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.4f}"
    except Exception:
        return str(value)


def _markdown_table(df: pd.DataFrame, cols: list[str]) -> str:
    if df.empty:
        return "_No rows._"
    sub = df[[c for c in cols if c in df.columns]].copy()
    lines = ["| " + " | ".join(sub.columns) + " |", "| " + " | ".join(["---"] * len(sub.columns)) + " |"]
    for row in sub.itertuples(index=False):
        lines.append("| " + " | ".join(_fmt(v) for v in row) + " |")
    return "\n".join(lines)


def _load_rows(path: Path) -> pd.DataFrame:
    usecols = [
        "sample_id",
        "canonical_drug_id",
        "drug_name",
        "canonical_smiles",
        "target_genes",
        "PATHWAY_NAME_NORMALIZED",
        "classification",
    ]
    return pd.read_csv(path, usecols=lambda c: c in usecols)


def train_lightgbm_groupcv(cfg: dict[str, Any], input_name: str, out_dir: Path, force: bool = False) -> tuple[np.ndarray, dict[str, Any]]:
    paths = pipeline_paths(cfg)
    oof_path = out_dir / "groupcv" / "oof" / "LightGBM.npy"
    result_path = out_dir / "groupcv" / "folds" / "LightGBM_groupcv.json"
    if oof_path.exists() and result_path.exists() and not force:
        return np.load(oof_path), _read_json(result_path)

    X = np.load(paths.processed_dir / f"X_{input_name}.npy").astype(np.float32)
    y = np.load(paths.processed_dir / "y_train.npy").astype(np.float32)
    rows = pd.read_csv(paths.processed_dir / "row_metadata.csv", usecols=["canonical_drug_id"])
    groups = rows["canonical_drug_id"].astype(str).to_numpy()
    seed = int(cfg["training"].get("random_state", cfg["project"].get("seed", 42)))
    n_splits = min(3, len(np.unique(groups)))
    if n_splits < 2:
        raise ValueError("At least two drug groups are required for GroupCV")

    oof = np.zeros(len(y), dtype=np.float32)
    fold_metrics = []
    start = time.time()
    for fold, (train_idx, valid_idx) in enumerate(GroupKFold(n_splits=n_splits).split(X, y, groups), start=1):
        model = make_model("LightGBM", seed + fold, len(y))
        if model is None:
            raise RuntimeError("LightGBM is unavailable")
        fold_start = time.time()
        model.fit(X[train_idx], y[train_idx])
        train_pred = np.asarray(model.predict(X[train_idx]), dtype=float)
        valid_pred = np.asarray(model.predict(X[valid_idx]), dtype=float)
        oof[valid_idx] = valid_pred.astype(np.float32)
        fold_metrics.append(
            {
                "fold": fold,
                "valid_drug_groups": int(pd.Series(groups[valid_idx]).nunique()),
                "train": regression_metrics(y[train_idx], train_pred),
                "valid": regression_metrics(y[valid_idx], valid_pred),
                "elapsed_sec": float(time.time() - fold_start),
            }
        )
        valid_metrics = fold_metrics[-1]["valid"]
        print(
            f"[no-lincs-mixed] LightGBM groupcv fold={fold} "
            f"groups={fold_metrics[-1]['valid_drug_groups']} "
            f"spearman={valid_metrics['spearman']:.4f} rmse={valid_metrics['rmse']:.4f}",
            flush=True,
        )

    metrics = regression_metrics(y, oof)
    result = {
        "input_set": input_name,
        "model": "LightGBM",
        "cv": "groupcv_by_drug",
        "fold_metrics": fold_metrics,
        "oof_metrics": metrics,
        "train_oof_spearman_gap": float(np.nanmean([m["train"]["spearman"] for m in fold_metrics]) - metrics["spearman"]),
        "prediction_variance": float(np.nanvar(oof)),
        "elapsed_sec": float(time.time() - start),
    }
    oof_path.parent.mkdir(parents=True, exist_ok=True)
    result_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(oof_path, oof)
    write_json(result_path, result)
    return oof, result


def load_random3_predictions(cfg: dict[str, Any], input_name: str) -> dict[str, np.ndarray]:
    paths = pipeline_paths(cfg)
    return {
        "CrossAttention": np.load(paths.results_dir / "additional_dl" / input_name / "oof" / "CrossAttention.npy"),
        "ResidualMLP": np.load(paths.results_dir / "additional_dl" / input_name / "oof" / "ResidualMLP.npy"),
        "LightGBM": np.load(paths.results_dir / "no_lincs" / "oof" / input_name / "LightGBM.npy"),
    }


def load_groupcv_predictions(cfg: dict[str, Any], input_name: str, out_dir: Path, force_lightgbm: bool = False) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    paths = pipeline_paths(cfg)
    preds = {
        "CrossAttention": np.load(paths.results_dir / "additional_dl" / input_name / "groupcv" / "oof" / "CrossAttention.npy"),
        "ResidualMLP": np.load(paths.results_dir / "additional_dl" / input_name / "groupcv" / "oof" / "ResidualMLP.npy"),
    }
    preds["LightGBM"], lightgbm_result = train_lightgbm_groupcv(cfg, input_name, out_dir, force=force_lightgbm)
    return preds, lightgbm_result


def summarize_ensemble(
    cfg: dict[str, Any],
    input_name: str,
    cv_name: str,
    preds: dict[str, np.ndarray],
    out_dir: Path,
) -> dict[str, Any]:
    paths = pipeline_paths(cfg)
    y = np.load(paths.processed_dir / "y_train.npy").astype(np.float32)
    rows = _load_rows(paths.processed_dir / "row_metadata.csv")

    individual = pd.DataFrame(
        [
            {
                "cv": cv_name,
                "model": name,
                **regression_metrics(y, pred),
                "prediction_variance": float(np.nanvar(pred)),
            }
            for name, pred in preds.items()
        ]
    ).sort_values("spearman", ascending=False)

    raw_weights = {row.model: max(0.0, float(row.spearman)) for row in individual.itertuples()}
    total = sum(raw_weights.values())
    weighted_weights = {k: (v / total if total > 0 else 1.0 / len(preds)) for k, v in raw_weights.items()}
    equal_weights = {k: 1.0 / len(preds) for k in preds}

    def combine(weights: dict[str, float]) -> np.ndarray:
        mixed = np.zeros(len(y), dtype=np.float32)
        for name, weight in weights.items():
            mixed += float(weight) * preds[name].astype(np.float32)
        return mixed

    weighted_pred = combine(weighted_weights)
    equal_pred = combine(equal_weights)
    ensemble = pd.DataFrame(
        [
            {"cv": cv_name, "ensemble": "spearman_weighted", **regression_metrics(y, weighted_pred), "prediction_variance": float(np.nanvar(weighted_pred))},
            {"cv": cv_name, "ensemble": "equal_weight", **regression_metrics(y, equal_pred), "prediction_variance": float(np.nanvar(equal_pred))},
        ]
    ).sort_values("spearman", ascending=False)

    diversity_rows = []
    names = list(preds)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            lp = preds[left]
            rp = preds[right]
            diversity_rows.append(
                {
                    "cv": cv_name,
                    "model_a": left,
                    "model_b": right,
                    "prediction_pearson_corr": safe_pearson(lp, rp),
                    "prediction_spearman_corr": safe_spearman(lp, rp),
                    "residual_pearson_corr": safe_pearson(y - lp, y - rp),
                    "mean_abs_prediction_gap": float(np.mean(np.abs(lp - rp))),
                }
            )
    diversity = pd.DataFrame(diversity_rows)

    scored_rows = rows.copy()
    scored_rows["mixed_pred_ln_ic50"] = weighted_pred
    scored_rows["mixed_score"] = -weighted_pred
    group_cols = ["canonical_drug_id", "drug_name", "canonical_smiles", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"]
    top = (
        scored_rows.groupby(group_cols, dropna=False)
        .agg(
            mean_pred_ln_ic50=("mixed_pred_ln_ic50", "mean"),
            mixed_score=("mixed_score", "mean"),
            screened_rows=("mixed_score", "size"),
        )
        .reset_index()
        .sort_values("mixed_score", ascending=False)
    )
    top["_drug_name_norm"] = top["drug_name"].map(lambda x: NAME_NORMALIZER.sub("", str(x).lower()))
    top = top.drop_duplicates("_drug_name_norm", keep="first").drop(columns=["_drug_name_norm"])
    top["rank"] = np.arange(1, len(top) + 1)
    top = top[["rank"] + [c for c in top.columns if c != "rank"]]

    cv_dir = out_dir / cv_name
    cv_dir.mkdir(parents=True, exist_ok=True)
    np.save(cv_dir / "spearman_weighted_ensemble_oof.npy", weighted_pred)
    np.save(cv_dir / "equal_weight_ensemble_oof.npy", equal_pred)
    write_table(individual, cv_dir / "individual_metrics.csv")
    write_table(ensemble, cv_dir / "ensemble_metrics.csv")
    write_table(diversity, cv_dir / "ensemble_diversity.csv")
    write_table(top.head(30), cv_dir / "mixed_ensemble_top30_drugs.csv")

    return {
        "cv": cv_name,
        "weights": {"spearman_weighted": weighted_weights, "equal_weight": equal_weights},
        "individual_metrics": individual.to_dict(orient="records"),
        "ensemble_metrics": ensemble.to_dict(orient="records"),
        "diversity_rows": int(len(diversity)),
        "top30_path": str(cv_dir / "mixed_ensemble_top30_drugs.csv"),
        "individual_table": individual,
        "ensemble_table": ensemble,
        "diversity_table": diversity,
        "top_table": top,
    }


def write_report(cfg: dict[str, Any], input_name: str, random_result: dict[str, Any], group_result: dict[str, Any], out_dir: Path) -> Path:
    paths = pipeline_paths(cfg)
    report = paths.root / "docs" / "THYROID_NO_LINCS_MIXED_RANDOM3_GROUPCV_ENSEMBLE_20260421.md"
    metric_cols = ["cv", "model", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20", "prediction_variance"]
    ensemble_cols = ["cv", "ensemble", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20", "prediction_variance"]
    diversity_cols = ["cv", "model_a", "model_b", "prediction_spearman_corr", "residual_pearson_corr", "mean_abs_prediction_gap"]
    top_cols = ["rank", "drug_name", "mixed_score", "mean_pred_ln_ic50", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"]
    text = f"""# Thyroid No-LINCS Mixed Ensemble - 2026-04-21

## 목적

LINCS feature를 제거한 입력셋에서 `CrossAttention + ResidualMLP + LightGBM` 조합을 random sample 3-fold와 drug-level GroupCV 기준으로 각각 평가했다.

## 설정

- Input set: `{input_name}`
- Removed feature family: LINCS
- Members: `CrossAttention`, `ResidualMLP`, `LightGBM`
- CV 1: random sample 3-fold OOF
- CV 2: GroupKFold 3-fold by `canonical_drug_id`
- Ensemble: Spearman-weighted average와 equal-weight average

## Random Sample 3-fold 개별 모델 성능

{_markdown_table(random_result["individual_table"], metric_cols)}

## Random Sample 3-fold 앙상블 성능

{_markdown_table(random_result["ensemble_table"], ensemble_cols)}

## Random Sample 3-fold Diversity

{_markdown_table(random_result["diversity_table"], diversity_cols)}

## GroupCV 개별 모델 성능

{_markdown_table(group_result["individual_table"], metric_cols)}

## GroupCV 앙상블 성능

{_markdown_table(group_result["ensemble_table"], ensemble_cols)}

## GroupCV Diversity

{_markdown_table(group_result["diversity_table"], diversity_cols)}

## GroupCV Spearman-weighted Top15 Drugs

{_markdown_table(group_result["top_table"].head(15), top_cols)}

## Weights

```json
{json.dumps({"random3": random_result["weights"], "groupcv": group_result["weights"]}, ensure_ascii=False, indent=2)}
```

## 산출물

- `{out_dir}/random3/individual_metrics.csv`
- `{out_dir}/random3/ensemble_metrics.csv`
- `{out_dir}/random3/ensemble_diversity.csv`
- `{out_dir}/random3/mixed_ensemble_top30_drugs.csv`
- `{out_dir}/groupcv/individual_metrics.csv`
- `{out_dir}/groupcv/ensemble_metrics.csv`
- `{out_dir}/groupcv/ensemble_diversity.csv`
- `{out_dir}/groupcv/mixed_ensemble_top30_drugs.csv`
- `reports/qc_no_lincs_mixed_random3_groupcv_ensemble_20260421.json`
"""
    report.write_text(text, encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Run no-LINCS CrossAttention + ResidualMLP + LightGBM random3/GroupCV ensemble")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--input-set", default=INPUT_NAME)
    parser.add_argument("--force-lightgbm-groupcv", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = pipeline_paths(cfg)
    out_dir = paths.results_dir / "mixed_no_lincs" / ENSEMBLE_NAME
    out_dir.mkdir(parents=True, exist_ok=True)

    random_preds = load_random3_predictions(cfg, args.input_set)
    random_result = summarize_ensemble(cfg, args.input_set, "random3", random_preds, out_dir)
    group_preds, lightgbm_group = load_groupcv_predictions(cfg, args.input_set, out_dir, force_lightgbm=args.force_lightgbm_groupcv)
    group_result = summarize_ensemble(cfg, args.input_set, "groupcv", group_preds, out_dir)

    result = {
        "ensemble_name": ENSEMBLE_NAME,
        "input_set": args.input_set,
        "members": list(MEMBERS),
        "removed_feature_family": "LINCS",
        "random3": {k: v for k, v in random_result.items() if not k.endswith("_table")},
        "groupcv": {k: v for k, v in group_result.items() if not k.endswith("_table")},
        "lightgbm_groupcv": lightgbm_group,
        "outputs": {
            "root": str(out_dir),
            "report": str(paths.root / "docs" / "THYROID_NO_LINCS_MIXED_RANDOM3_GROUPCV_ENSEMBLE_20260421.md"),
        },
    }
    write_json(out_dir / "mixed_no_lincs_random3_groupcv_results.json", result)
    write_json(paths.reports_dir / "qc_no_lincs_mixed_random3_groupcv_ensemble_20260421.json", result)
    report = write_report(cfg, args.input_set, random_result, group_result, out_dir)

    print(
        json.dumps(
            {
                "status": "completed",
                "input_set": args.input_set,
                "random3_best_ensemble": random_result["ensemble_table"].iloc[0].to_dict(),
                "groupcv_best_ensemble": group_result["ensemble_table"].iloc[0].to_dict(),
                "report": str(report),
                "output_root": str(out_dir),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
