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


ENSEMBLE_NAME = "crossattention_residualmlp_lightgbm_groupcv"
INPUT_NAME = "numeric_strong_context_smiles_pan_lincs"


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


def load_rows_minimal(path: Path) -> pd.DataFrame:
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
    oof_path = out_dir / "oof" / "LightGBM.npy"
    result_path = out_dir / "folds" / "LightGBM_groupcv.json"
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
        tr = regression_metrics(y[train_idx], train_pred)
        va = regression_metrics(y[valid_idx], valid_pred)
        fold_metrics.append(
            {
                "fold": fold,
                "valid_drug_groups": int(pd.Series(groups[valid_idx]).nunique()),
                "train": tr,
                "valid": va,
                "elapsed_sec": float(time.time() - fold_start),
            }
        )
        print(
            f"[mixed-groupcv] LightGBM fold={fold} groups={fold_metrics[-1]['valid_drug_groups']} "
            f"spearman={va['spearman']:.4f} rmse={va['rmse']:.4f}",
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
    (out_dir / "oof").mkdir(parents=True, exist_ok=True)
    (out_dir / "folds").mkdir(parents=True, exist_ok=True)
    np.save(oof_path, oof)
    write_json(result_path, result)
    return oof, result


def load_dl_groupcv_oof(cfg: dict[str, Any], input_name: str, model_name: str) -> np.ndarray:
    paths = pipeline_paths(cfg)
    path = paths.results_dir / "additional_dl" / input_name / "groupcv" / "oof" / f"{model_name}.npy"
    if not path.exists():
        raise FileNotFoundError(path)
    return np.load(path)


def build_mixed_ensemble(cfg: dict[str, Any], input_name: str, force_lightgbm: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    paths = pipeline_paths(cfg)
    y = np.load(paths.processed_dir / "y_train.npy").astype(np.float32)
    rows = load_rows_minimal(paths.processed_dir / "row_metadata.csv")
    out_dir = paths.results_dir / "mixed_groupcv" / ENSEMBLE_NAME
    out_dir.mkdir(parents=True, exist_ok=True)

    preds: dict[str, np.ndarray] = {
        "CrossAttention": load_dl_groupcv_oof(cfg, input_name, "CrossAttention"),
        "ResidualMLP": load_dl_groupcv_oof(cfg, input_name, "ResidualMLP"),
    }
    preds["LightGBM"], lightgbm_result = train_lightgbm_groupcv(cfg, input_name, out_dir, force=force_lightgbm)

    individual_rows = []
    for name, pred in preds.items():
        individual_rows.append(
            {
                "model": name,
                **regression_metrics(y, pred),
                "prediction_variance": float(np.nanvar(pred)),
            }
        )
    individual = pd.DataFrame(individual_rows).sort_values("spearman", ascending=False)

    weighted_weights = {
        row.model: max(0.0, float(row.spearman)) for row in individual.itertuples()
    }
    total = sum(weighted_weights.values())
    weighted_weights = {k: (v / total if total > 0 else 1.0 / len(preds)) for k, v in weighted_weights.items()}
    equal_weights = {k: 1.0 / len(preds) for k in preds}

    def combine(weights: dict[str, float]) -> np.ndarray:
        out = np.zeros(len(y), dtype=np.float32)
        for name, weight in weights.items():
            out += float(weight) * preds[name].astype(np.float32)
        return out

    weighted_pred = combine(weighted_weights)
    equal_pred = combine(equal_weights)
    ensemble_rows = [
        {"ensemble": "spearman_weighted", **regression_metrics(y, weighted_pred), "prediction_variance": float(np.nanvar(weighted_pred))},
        {"ensemble": "equal_weight", **regression_metrics(y, equal_pred), "prediction_variance": float(np.nanvar(equal_pred))},
    ]
    ensemble_metrics = pd.DataFrame(ensemble_rows).sort_values("spearman", ascending=False)

    diversity_rows = []
    names = list(preds)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            lp = preds[left]
            rp = preds[right]
            diversity_rows.append(
                {
                    "model_a": left,
                    "model_b": right,
                    "prediction_pearson_corr": safe_pearson(lp, rp),
                    "prediction_spearman_corr": safe_spearman(lp, rp),
                    "residual_pearson_corr": safe_pearson(y - lp, y - rp),
                    "mean_abs_prediction_gap": float(np.mean(np.abs(lp - rp))),
                }
            )
    diversity = pd.DataFrame(diversity_rows)

    candidate_rows = rows.copy()
    candidate_rows["mixed_groupcv_pred_ln_ic50"] = weighted_pred
    candidate_rows["mixed_groupcv_score"] = -weighted_pred
    group_cols = ["canonical_drug_id", "drug_name", "canonical_smiles", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"]
    top = (
        candidate_rows.groupby(group_cols, dropna=False)
        .agg(
            mean_pred_ln_ic50=("mixed_groupcv_pred_ln_ic50", "mean"),
            mixed_groupcv_score=("mixed_groupcv_score", "mean"),
            screened_rows=("mixed_groupcv_score", "size"),
        )
        .reset_index()
        .sort_values("mixed_groupcv_score", ascending=False)
    )
    top["_drug_name_norm"] = top["drug_name"].map(lambda x: NAME_NORMALIZER.sub("", str(x).lower()))
    top = top.drop_duplicates("_drug_name_norm", keep="first").drop(columns=["_drug_name_norm"])
    top["rank"] = np.arange(1, len(top) + 1)
    top = top[["rank"] + [c for c in top.columns if c != "rank"]]

    np.save(out_dir / "spearman_weighted_ensemble_oof.npy", weighted_pred)
    np.save(out_dir / "equal_weight_ensemble_oof.npy", equal_pred)
    write_table(individual, out_dir / "individual_groupcv_metrics.csv")
    write_table(ensemble_metrics, out_dir / "ensemble_groupcv_metrics.csv")
    write_table(diversity, out_dir / "ensemble_diversity.csv")
    write_table(top.head(30), out_dir / "mixed_groupcv_ensemble_top30_drugs.csv")

    result = {
        "ensemble_name": ENSEMBLE_NAME,
        "input_set": input_name,
        "cv": "groupcv_by_drug",
        "members": list(preds),
        "weights": {
            "spearman_weighted": weighted_weights,
            "equal_weight": equal_weights,
        },
        "individual_metrics": individual.to_dict(orient="records"),
        "ensemble_metrics": ensemble_metrics.to_dict(orient="records"),
        "lightgbm_groupcv": lightgbm_result,
        "diversity_rows": int(len(diversity)),
        "outputs": {
            "individual_metrics": str(out_dir / "individual_groupcv_metrics.csv"),
            "ensemble_metrics": str(out_dir / "ensemble_groupcv_metrics.csv"),
            "diversity": str(out_dir / "ensemble_diversity.csv"),
            "top30": str(out_dir / "mixed_groupcv_ensemble_top30_drugs.csv"),
        },
    }
    write_json(out_dir / "mixed_groupcv_ensemble_results.json", result)
    write_json(paths.reports_dir / "qc_mixed_groupcv_crossattention_residualmlp_lightgbm_20260421.json", result)
    write_report(cfg, input_name, individual, ensemble_metrics, diversity, top, result)
    return individual, ensemble_metrics, result


def write_report(
    cfg: dict[str, Any],
    input_name: str,
    individual: pd.DataFrame,
    ensemble_metrics: pd.DataFrame,
    diversity: pd.DataFrame,
    top: pd.DataFrame,
    result: dict[str, Any],
) -> Path:
    paths = pipeline_paths(cfg)
    report = paths.root / "docs" / "THYROID_MIXED_GROUPCV_CROSSATTENTION_RESIDUALMLP_LIGHTGBM_20260421.md"
    text = f"""# Thyroid Mixed GroupCV Ensemble - 2026-04-21

## 목적

`CrossAttention + ResidualMLP + LightGBM` 조합을 drug-level GroupCV 기준으로 평가했다. CrossAttention/ResidualMLP는 기존 추가 DL GroupCV OOF를 재사용했고, LightGBM은 같은 `canonical_drug_id` GroupKFold split으로 새로 학습했다.

## 설정

- Input set: `{input_name}`
- CV: GroupKFold 3-fold by `canonical_drug_id`
- Members: `CrossAttention`, `ResidualMLP`, `LightGBM`
- Ensemble: Spearman-weighted average와 equal-weight average를 함께 산출

## 개별 모델 GroupCV 성능

{_markdown_table(individual, ["model", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20", "prediction_variance"])}

## 앙상블 GroupCV 성능

{_markdown_table(ensemble_metrics, ["ensemble", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20", "prediction_variance"])}

## Diversity

{_markdown_table(diversity, ["model_a", "model_b", "prediction_spearman_corr", "residual_pearson_corr", "mean_abs_prediction_gap"])}

## Spearman-weighted Top15 Drugs

{_markdown_table(top.head(15), ["rank", "drug_name", "mixed_groupcv_score", "mean_pred_ln_ic50", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"])}

## Weights

```json
{json.dumps(result["weights"], ensure_ascii=False, indent=2)}
```

## 산출물

- `results/mixed_groupcv/{ENSEMBLE_NAME}/individual_groupcv_metrics.csv`
- `results/mixed_groupcv/{ENSEMBLE_NAME}/ensemble_groupcv_metrics.csv`
- `results/mixed_groupcv/{ENSEMBLE_NAME}/ensemble_diversity.csv`
- `results/mixed_groupcv/{ENSEMBLE_NAME}/mixed_groupcv_ensemble_top30_drugs.csv`
- `results/mixed_groupcv/{ENSEMBLE_NAME}/mixed_groupcv_ensemble_results.json`
- `reports/qc_mixed_groupcv_crossattention_residualmlp_lightgbm_20260421.json`
"""
    report.write_text(text, encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Run CrossAttention + ResidualMLP + LightGBM GroupCV ensemble")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--input-set", default=INPUT_NAME)
    parser.add_argument("--force-lightgbm", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    individual, ensemble_metrics, result = build_mixed_ensemble(cfg, args.input_set, force_lightgbm=args.force_lightgbm)
    print(
        json.dumps(
            {
                "status": "completed",
                "ensemble_name": ENSEMBLE_NAME,
                "best_individual": individual.iloc[0].to_dict() if not individual.empty else {},
                "best_ensemble": ensemble_metrics.iloc[0].to_dict() if not ensemble_metrics.empty else {},
                "outputs": result["outputs"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
