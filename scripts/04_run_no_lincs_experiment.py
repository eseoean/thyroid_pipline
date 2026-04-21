from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold, KFold

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import (
    NAME_NORMALIZER,
    load_config,
    make_model,
    pipeline_paths,
    read_table,
    regression_metrics,
    safe_pearson,
    safe_spearman,
    write_json,
    write_table,
)


INPUT_DEFS = {
    "numeric": ("X_numeric.npy", "numeric_feature_names.json"),
    "numeric_smiles": ("X_numeric_smiles.npy", "numeric_smiles_feature_names.json"),
    "numeric_strong_context_smiles": (
        "X_numeric_strong_context_smiles.npy",
        "numeric_strong_context_smiles_feature_names.json",
    ),
}


def _read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _write_feature_names(path: Path, names: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(names, f, ensure_ascii=False, indent=2)


def build_no_lincs_inputs(cfg: dict[str, Any]) -> tuple[dict[str, Path], pd.DataFrame]:
    paths = pipeline_paths(cfg)
    input_paths: dict[str, Path] = {}
    rows = []
    for base_name, (array_name, feature_name_file) in INPUT_DEFS.items():
        x_path = paths.processed_dir / array_name
        names_path = paths.processed_dir / feature_name_file
        X = np.load(x_path)
        names = [str(x) for x in _read_json(names_path)]
        if X.shape[1] != len(names):
            raise ValueError(f"{base_name}: feature count mismatch: X={X.shape[1]} names={len(names)}")

        keep_mask = np.array(["lincs" not in name.lower() for name in names], dtype=bool)
        removed_names = [name for name, keep in zip(names, keep_mask) if not keep]
        kept_names = [name for name, keep in zip(names, keep_mask) if keep]
        no_lincs_name = f"{base_name}_no_lincs"
        out_path = paths.processed_dir / f"X_{no_lincs_name}.npy"
        out_names_path = paths.processed_dir / f"{no_lincs_name}_feature_names.json"
        np.save(out_path, X[:, keep_mask].astype(np.float32))
        _write_feature_names(out_names_path, kept_names)
        input_paths[no_lincs_name] = out_path
        rows.append(
            {
                "input_set": no_lincs_name,
                "source_input_set": base_name,
                "rows": int(X.shape[0]),
                "source_features": int(X.shape[1]),
                "removed_lincs_features": int(len(removed_names)),
                "remaining_features": int(len(kept_names)),
                "array_path": str(out_path),
                "feature_names_path": str(out_names_path),
            }
        )

    summary = pd.DataFrame(rows)
    write_table(summary, paths.reports_dir / "no_lincs_input_summary_20260421.csv")
    write_json(paths.reports_dir / "no_lincs_input_summary_20260421.json", summary.to_dict(orient="records"))
    return input_paths, summary


def _fold_coverage(rows: pd.DataFrame, splits: list[tuple[np.ndarray, np.ndarray]], col: str) -> list[dict[str, int]]:
    return [{"fold": i, f"unique_{col}": int(rows.iloc[valid_idx][col].nunique())} for i, (_, valid_idx) in enumerate(splits, start=1)]


def train_no_lincs(cfg: dict[str, Any], input_paths: dict[str, Path]) -> dict[str, pd.DataFrame]:
    paths = pipeline_paths(cfg)
    y = np.load(paths.processed_dir / "y_train.npy")
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    seed = int(cfg["training"].get("random_state", cfg["project"].get("seed", 42)))
    n_splits = min(int(cfg["training"].get("n_splits", 3)), len(y))
    models = list(cfg["training"]["models"])
    out_root = paths.results_dir / "no_lincs"
    summaries: dict[str, pd.DataFrame] = {}
    qc: dict[str, Any] = {"step": "no_lincs_random3_training", "input_sets": {}, "skipped_models": []}

    for input_name, input_path in input_paths.items():
        X = np.load(input_path)
        out_dir = out_root / "random3" / input_name
        oof_dir = out_root / "oof" / input_name
        out_dir.mkdir(parents=True, exist_ok=True)
        oof_dir.mkdir(parents=True, exist_ok=True)
        kfold = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
        splits = list(kfold.split(X))
        fold_rows = []
        summary_rows = []
        for model_name in models:
            if make_model(model_name, seed, len(y)) is None:
                qc["skipped_models"].append(
                    {"input_set": input_name, "model": model_name, "reason": "dependency_missing_or_not_implemented"}
                )
                continue

            oof = np.zeros(len(y), dtype=np.float32)
            train_metrics = []
            valid_metrics = []
            for fold, (train_idx, valid_idx) in enumerate(splits, start=1):
                model = make_model(model_name, seed + fold, len(y))
                model.fit(X[train_idx], y[train_idx])
                train_pred = np.asarray(model.predict(X[train_idx]), dtype=float)
                valid_pred = np.asarray(model.predict(X[valid_idx]), dtype=float)
                oof[valid_idx] = valid_pred.astype(np.float32)
                tr = regression_metrics(y[train_idx], train_pred)
                va = regression_metrics(y[valid_idx], valid_pred)
                train_metrics.append(tr)
                valid_metrics.append(va)
                fold_rows.append(
                    {
                        "input_set": input_name,
                        "model": model_name,
                        "fold": fold,
                        "train_rows": int(len(train_idx)),
                        "valid_rows": int(len(valid_idx)),
                        "train_spearman": tr["spearman"],
                        "valid_spearman": va["spearman"],
                        "valid_rmse": va["rmse"],
                    }
                )

            oof_metrics = regression_metrics(y, oof)
            train_gap = float(np.nanmean([m["spearman"] for m in train_metrics]) - oof_metrics["spearman"])
            payload = {
                "input_set": input_name,
                "model": model_name,
                "removed_feature_family": "LINCS",
                "fold_metrics": {"train": train_metrics, "valid": valid_metrics},
                "oof_metrics": oof_metrics,
                "train_oof_spearman_gap": train_gap,
                "prediction_variance": float(np.nanvar(oof)),
            }
            np.save(oof_dir / f"{model_name}.npy", oof)
            write_json(out_dir / f"{input_name}_{model_name}_random3.json", payload)
            summary_rows.append(
                {
                    "input_set": input_name,
                    "model": model_name,
                    **oof_metrics,
                    "train_oof_spearman_gap": train_gap,
                    "prediction_variance": float(np.nanvar(oof)),
                }
            )

        summary = pd.DataFrame(summary_rows).sort_values("spearman", ascending=False, na_position="last")
        write_table(summary, out_root / f"{input_name}_metrics_summary.csv")
        write_table(pd.DataFrame(fold_rows), out_root / f"{input_name}_fold_metrics.csv")
        summaries[input_name] = summary
        qc["input_sets"][input_name] = {
            "shape": [int(X.shape[0]), int(X.shape[1])],
            "models_completed": summary["model"].tolist() if len(summary) else [],
            "fold_count": int(n_splits),
            "drug_coverage_per_fold": _fold_coverage(rows, splits, "canonical_drug_id") if rows is not None else [],
            "cell_line_coverage_per_fold": _fold_coverage(rows, splits, "sample_id") if rows is not None else [],
        }

    write_json(paths.reports_dir / "qc_no_lincs_training_20260421.json", qc)
    return summaries


def run_groupcv_no_lincs(cfg: dict[str, Any], input_name: str = "numeric_strong_context_smiles_no_lincs") -> dict[str, Any]:
    paths = pipeline_paths(cfg)
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    if rows is None:
        return {"status": "skipped", "reason": "row_metadata_missing"}
    X = np.load(paths.processed_dir / f"X_{input_name}.npy")
    y = np.load(paths.processed_dir / "y_train.npy")
    groups = rows["canonical_drug_id"].astype(str).to_numpy()
    n_splits = min(3, len(np.unique(groups)))
    if n_splits < 2:
        return {"status": "skipped", "reason": "not_enough_drug_groups"}

    model_name = "ExtraTrees"
    seed = int(cfg["project"].get("seed", 42))
    oof = np.zeros(len(y), dtype=np.float32)
    fold_metrics = []
    for fold, (train_idx, valid_idx) in enumerate(GroupKFold(n_splits=n_splits).split(X, y, groups), start=1):
        model = make_model(model_name, seed + fold, len(y))
        model.fit(X[train_idx], y[train_idx])
        pred = np.asarray(model.predict(X[valid_idx]), dtype=float)
        oof[valid_idx] = pred.astype(np.float32)
        fold_metrics.append({"fold": fold, **regression_metrics(y[valid_idx], pred)})

    out_dir = paths.results_dir / "no_lincs" / "groupcv_stress_test"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"{input_name}_{model_name}_groupcv_oof.npy", oof)
    result = {
        "status": "completed",
        "input_set": input_name,
        "model": model_name,
        "removed_feature_family": "LINCS",
        "fold_metrics": fold_metrics,
        "oof_metrics": regression_metrics(y, oof),
    }
    write_json(out_dir / f"{input_name}_{model_name}_groupcv.json", result)
    return result


def run_ensemble_no_lincs(cfg: dict[str, Any], input_name: str = "numeric_strong_context_smiles_no_lincs") -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    paths = pipeline_paths(cfg)
    out_root = paths.results_dir / "no_lincs"
    y = np.load(paths.processed_dir / "y_train.npy")
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    metrics = read_table(out_root / f"{input_name}_metrics_summary.csv")
    if metrics is None or metrics.empty:
        raise FileNotFoundError(f"No no-LINCS metrics found for {input_name}")

    preds = {}
    weights = {}
    for row in metrics.itertuples():
        model_name = row.model
        path = out_root / "oof" / input_name / f"{model_name}.npy"
        if path.exists():
            preds[model_name] = np.load(path)
            score = getattr(row, "spearman", np.nan)
            weights[model_name] = max(0.0, float(score)) if np.isfinite(score) else 0.0
    if not preds:
        raise ValueError("No no-LINCS OOF predictions found")
    weight_sum = sum(weights.values())
    weights = {k: (v / weight_sum if weight_sum > 0 else 1.0 / len(preds)) for k, v in weights.items()}

    ensemble_pred = np.zeros(len(y), dtype=np.float32)
    for model_name, pred in preds.items():
        ensemble_pred += weights[model_name] * pred

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

    if rows is None:
        rows = pd.DataFrame({"canonical_drug_id": np.arange(len(y)), "drug_name": np.arange(len(y))})
    candidate_rows = rows.copy()
    candidate_rows["ensemble_pred_ln_ic50"] = ensemble_pred
    candidate_rows["ensemble_score"] = -ensemble_pred
    group_cols = ["canonical_drug_id", "drug_name", "canonical_smiles", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"]
    available = [c for c in group_cols if c in candidate_rows.columns]
    top = (
        candidate_rows.groupby(available, dropna=False)
        .agg(
            mean_pred_ln_ic50=("ensemble_pred_ln_ic50", "mean"),
            ensemble_score=("ensemble_score", "mean"),
            screened_rows=("ensemble_score", "size"),
        )
        .reset_index()
        .sort_values("ensemble_score", ascending=False)
    )
    if "drug_name" in top.columns:
        top["_drug_name_norm"] = top["drug_name"].map(lambda x: NAME_NORMALIZER.sub("", str(x).lower()))
        top = top.drop_duplicates("_drug_name_norm", keep="first").drop(columns=["_drug_name_norm"])
    top["rank"] = np.arange(1, len(top) + 1)
    top = top[["rank"] + [c for c in top.columns if c != "rank"]]

    out_dir = out_root / "ensemble"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"{input_name}_weighted_ensemble_oof.npy", ensemble_pred)
    write_table(top.head(30), out_dir / "thyroid_no_lincs_ensemble_top30_drugs.csv")
    write_table(diversity, out_dir / "thyroid_no_lincs_ensemble_diversity.csv")
    result = {
        "input_set": input_name,
        "removed_feature_family": "LINCS",
        "weights": weights,
        "ensemble_metrics": regression_metrics(y, ensemble_pred),
        "best_single_model": metrics.iloc[0].to_dict(),
        "diversity_rows": int(len(diversity)),
    }
    write_json(out_dir / "thyroid_no_lincs_ensemble_results.json", result)
    write_json(paths.reports_dir / "qc_no_lincs_ensemble_20260421.json", result)
    return top, diversity, result


def compare_with_baseline(cfg: dict[str, Any], input_summary: pd.DataFrame, groupcv: dict[str, Any], ensemble: dict[str, Any]) -> pd.DataFrame:
    paths = pipeline_paths(cfg)
    rows = []
    for no_lincs_input in input_summary["input_set"].tolist():
        source_input = str(input_summary.loc[input_summary["input_set"].eq(no_lincs_input), "source_input_set"].iloc[0])
        baseline = read_table(paths.results_dir / f"{source_input}_metrics_summary.csv")
        no_lincs = read_table(paths.results_dir / "no_lincs" / f"{no_lincs_input}_metrics_summary.csv")
        if baseline is not None and not baseline.empty:
            best = baseline.sort_values("spearman", ascending=False).iloc[0]
            rows.append(
                {
                    "input_set": source_input,
                    "variant": "baseline_with_lincs",
                    "model": best["model"],
                    "spearman": best["spearman"],
                    "rmse": best["rmse"],
                    "r2": best["r2"],
                    "ndcg_at_20": best.get("ndcg_at_20", np.nan),
                }
            )
        if no_lincs is not None and not no_lincs.empty:
            best = no_lincs.sort_values("spearman", ascending=False).iloc[0]
            rows.append(
                {
                    "input_set": source_input,
                    "variant": "no_lincs",
                    "model": best["model"],
                    "spearman": best["spearman"],
                    "rmse": best["rmse"],
                    "r2": best["r2"],
                    "ndcg_at_20": best.get("ndcg_at_20", np.nan),
                }
            )

    comparison = pd.DataFrame(rows)
    write_table(comparison, paths.results_dir / "no_lincs" / "baseline_vs_no_lincs_best_model_comparison.csv")

    groupcv_path = paths.results_dir / "groupcv_stress_test" / "numeric_strong_context_smiles_ExtraTrees_groupcv.json"
    baseline_groupcv = _read_json(groupcv_path) if groupcv_path.exists() else {}
    group_rows = []
    if baseline_groupcv:
        group_rows.append(
            {
                "variant": "baseline_with_lincs",
                "input_set": baseline_groupcv.get("input_set"),
                "model": baseline_groupcv.get("model"),
                **baseline_groupcv.get("oof_metrics", {}),
            }
        )
    if groupcv:
        group_rows.append(
            {
                "variant": "no_lincs",
                "input_set": groupcv.get("input_set"),
                "model": groupcv.get("model"),
                **groupcv.get("oof_metrics", {}),
            }
        )
    write_table(pd.DataFrame(group_rows), paths.results_dir / "no_lincs" / "baseline_vs_no_lincs_groupcv_comparison.csv")

    _write_markdown_report(cfg, input_summary, comparison, pd.DataFrame(group_rows), ensemble)
    return comparison


def _fmt(value: Any) -> str:
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.4f}"
    except Exception:
        return str(value)


def _markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for row in df.itertuples(index=False):
        lines.append("| " + " | ".join(_fmt(v) for v in row) + " |")
    return "\n".join(lines)


def _write_markdown_report(
    cfg: dict[str, Any],
    input_summary: pd.DataFrame,
    comparison: pd.DataFrame,
    groupcv_comparison: pd.DataFrame,
    ensemble: dict[str, Any],
) -> None:
    paths = pipeline_paths(cfg)
    report_path = paths.root / "docs" / "THYROID_NO_LINCS_EXPERIMENT_20260421.md"
    metrics = ensemble.get("ensemble_metrics", {})
    best = ensemble.get("best_single_model", {})
    text = f"""# Thyroid No-LINCS Ablation Experiment - 2026-04-21

## 목적

갑상선암 파이프라인에서 BRCA/MCF7 기반 LINCS feature가 질환 특이적이지 않다는 리스크가 있어, LINCS 관련 컬럼을 모두 제거한 입력셋으로 random sample 3-fold 학습과 GroupCV stress test를 다시 수행했다.

## 제거 기준

feature name에 `lincs`가 포함된 모든 컬럼을 제거했다. 여기에는 `drug__lincs__*` signature feature와 `drug__has_lincs_signature`, `drug__lincs_signature_source_direct`, `drug__lincs_signature_source_recovered` availability/source flag가 포함된다.

## 입력셋 차원 변화

{_markdown_table(input_summary)}

## Random Sample 3-fold Best Model 비교

{_markdown_table(comparison)}

## GroupCV Stress Test 비교

{_markdown_table(groupcv_comparison)}

## No-LINCS Ensemble

- Primary input set: `numeric_strong_context_smiles_no_lincs`
- Best single model: `{best.get("model", "")}`
- Ensemble Spearman: `{_fmt(metrics.get("spearman"))}`
- Ensemble RMSE: `{_fmt(metrics.get("rmse"))}`
- Ensemble R2: `{_fmt(metrics.get("r2"))}`

## 산출물

- `data/X_numeric_no_lincs.npy`
- `data/X_numeric_smiles_no_lincs.npy`
- `data/X_numeric_strong_context_smiles_no_lincs.npy`
- `results/no_lincs/random3/`
- `results/no_lincs/groupcv_stress_test/`
- `results/no_lincs/ensemble/`
- `reports/qc_no_lincs_training_20260421.json`
- `reports/qc_no_lincs_ensemble_20260421.json`
"""
    report_path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run thyroid no-LINCS ablation experiment")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    args = parser.parse_args()
    cfg = load_config(args.config)
    input_paths, input_summary = build_no_lincs_inputs(cfg)
    train_no_lincs(cfg, input_paths)
    groupcv = run_groupcv_no_lincs(cfg)
    _, _, ensemble = run_ensemble_no_lincs(cfg)
    compare_with_baseline(cfg, input_summary, groupcv, ensemble)
    print(json.dumps({"status": "completed", "input_sets": list(input_paths)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
