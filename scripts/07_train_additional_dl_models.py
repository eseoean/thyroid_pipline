from __future__ import annotations

import argparse
import copy
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import KFold
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import (  # noqa: E402
    NAME_NORMALIZER,
    load_config,
    pipeline_paths,
    read_json,
    read_table,
    regression_metrics,
    safe_pearson,
    safe_spearman,
    write_json,
    write_table,
)


MODEL_DEFAULTS = {
    "ResidualMLP": {"epochs": 80, "batch_size": 256, "lr": 8e-4, "weight_decay": 1e-4, "patience": 12},
    "WideDeep": {"epochs": 80, "batch_size": 256, "lr": 8e-4, "weight_decay": 1e-4, "patience": 12},
    "CrossAttention": {"epochs": 70, "batch_size": 256, "lr": 6e-4, "weight_decay": 1e-4, "patience": 10},
}


def _json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _torch_device(name: str) -> torch.device:
    if name != "auto":
        return torch.device(name)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def set_seed(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


class ResidualBlock(nn.Module):
    def __init__(self, hidden: int, dropout: float) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.LayerNorm(hidden),
            nn.Linear(hidden, hidden * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden * 2, hidden),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class ResidualMLP(nn.Module):
    def __init__(self, in_dim: int, hidden: int = 512, n_blocks: int = 3, dropout: float = 0.25) -> None:
        super().__init__()
        self.in_proj = nn.Sequential(nn.Linear(in_dim, hidden), nn.LayerNorm(hidden), nn.GELU(), nn.Dropout(dropout))
        self.blocks = nn.Sequential(*[ResidualBlock(hidden, dropout) for _ in range(n_blocks)])
        self.head = nn.Sequential(nn.LayerNorm(hidden), nn.Linear(hidden, hidden // 2), nn.GELU(), nn.Dropout(dropout), nn.Linear(hidden // 2, 1))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = self.in_proj(x)
        h = self.blocks(h)
        return self.head(h).squeeze(-1)


class WideDeep(nn.Module):
    def __init__(self, in_dim: int, hidden: int = 512, dropout: float = 0.25) -> None:
        super().__init__()
        self.wide = nn.Linear(in_dim, 1)
        self.deep = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.LayerNorm(hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            ResidualBlock(hidden, dropout),
            ResidualBlock(hidden, dropout),
            nn.LayerNorm(hidden),
            nn.Linear(hidden, hidden // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden // 2, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return (self.wide(x) + self.deep(x)).squeeze(-1)


class CrossAttention(nn.Module):
    def __init__(self, in_dim: int, block_indices: dict[str, list[int]], d_model: int = 128, nhead: int = 4, dropout: float = 0.20) -> None:
        super().__init__()
        self.block_names = [name for name, idx in block_indices.items() if idx]
        if not self.block_names:
            raise ValueError("CrossAttention requires at least one non-empty feature block")
        self.projectors = nn.ModuleDict()
        for name in self.block_names:
            idx = torch.tensor(block_indices[name], dtype=torch.long)
            self.register_buffer(f"idx_{name}", idx, persistent=False)
            self.projectors[name] = nn.Sequential(
                nn.Linear(len(block_indices[name]), d_model),
                nn.LayerNorm(d_model),
                nn.GELU(),
                nn.Dropout(dropout),
            )
        self.query = nn.Parameter(torch.zeros(1, 1, d_model))
        nn.init.normal_(self.query, std=0.02)
        self.cross_attn = nn.MultiheadAttention(d_model, nhead, dropout=dropout, batch_first=True)
        self.head = nn.Sequential(
            nn.LayerNorm(d_model * 2),
            nn.Linear(d_model * 2, d_model),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        tokens = []
        for name in self.block_names:
            idx = getattr(self, f"idx_{name}")
            tokens.append(self.projectors[name](x.index_select(1, idx)))
        kv = torch.stack(tokens, dim=1)
        query = self.query.expand(x.shape[0], -1, -1)
        attended, _ = self.cross_attn(query, kv, kv, need_weights=False)
        pooled = kv.mean(dim=1)
        h = torch.cat([attended.squeeze(1), pooled], dim=1)
        return self.head(h).squeeze(-1)


def feature_blocks(feature_names: list[str]) -> dict[str, list[int]]:
    blocks = {
        "sample": [],
        "drug_lincs": [],
        "drug_morgan": [],
        "drug_other": [],
        "strong_context": [],
        "smiles": [],
        "other": [],
    }
    for i, name in enumerate(feature_names):
        text = str(name)
        lower = text.lower()
        if text.startswith("sample__"):
            blocks["sample"].append(i)
        elif "lincs" in lower:
            blocks["drug_lincs"].append(i)
        elif text.startswith("drug_morgan_"):
            blocks["drug_morgan"].append(i)
        elif text.startswith("drug_") or text.startswith("drug__"):
            blocks["drug_other"].append(i)
        elif text.startswith("smiles_svd_"):
            blocks["smiles"].append(i)
        elif text.startswith("strong_context") or "=" in text:
            blocks["strong_context"].append(i)
        else:
            blocks["other"].append(i)
    return {k: v for k, v in blocks.items() if v}


def make_model(name: str, in_dim: int, blocks: dict[str, list[int]]) -> nn.Module:
    if name == "ResidualMLP":
        return ResidualMLP(in_dim)
    if name == "WideDeep":
        return WideDeep(in_dim)
    if name == "CrossAttention":
        return CrossAttention(in_dim, blocks)
    raise ValueError(f"Unknown DL model: {name}")


@dataclass
class FoldResult:
    fold: int
    best_epoch: int
    best_valid_rmse: float
    train_metrics: dict[str, float]
    valid_metrics: dict[str, float]
    elapsed_sec: float


def standardize_fold(X: np.ndarray, train_idx: np.ndarray, valid_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    X_train = X[train_idx].astype(np.float32, copy=True)
    X_valid = X[valid_idx].astype(np.float32, copy=True)
    mean = np.nanmean(X_train, axis=0, dtype=np.float64).astype(np.float32)
    std = np.nanstd(X_train, axis=0, dtype=np.float64).astype(np.float32)
    std[~np.isfinite(std) | (std < 1e-6)] = 1.0
    mean[~np.isfinite(mean)] = 0.0
    X_train = (np.nan_to_num(X_train, nan=mean) - mean) / std
    X_valid = (np.nan_to_num(X_valid, nan=mean) - mean) / std
    return X_train.astype(np.float32), X_valid.astype(np.float32)


def train_one_fold(
    model_name: str,
    X: np.ndarray,
    y: np.ndarray,
    train_idx: np.ndarray,
    valid_idx: np.ndarray,
    blocks: dict[str, list[int]],
    device: torch.device,
    seed: int,
    params: dict[str, Any],
) -> tuple[np.ndarray, FoldResult]:
    set_seed(seed)
    start = time.time()
    X_train, X_valid = standardize_fold(X, train_idx, valid_idx)
    y_train_raw = y[train_idx].astype(np.float32)
    y_valid_raw = y[valid_idx].astype(np.float32)
    y_mean = float(np.nanmean(y_train_raw))
    y_std = float(np.nanstd(y_train_raw))
    if not np.isfinite(y_std) or y_std < 1e-6:
        y_std = 1.0
    y_train = ((y_train_raw - y_mean) / y_std).astype(np.float32)
    y_valid = ((y_valid_raw - y_mean) / y_std).astype(np.float32)

    train_ds = TensorDataset(torch.from_numpy(X_train), torch.from_numpy(y_train))
    valid_x = torch.from_numpy(X_valid).to(device)
    model = make_model(model_name, X.shape[1], blocks).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=float(params["lr"]), weight_decay=float(params["weight_decay"]))
    loss_fn = nn.SmoothL1Loss()
    loader = DataLoader(train_ds, batch_size=int(params["batch_size"]), shuffle=True, drop_last=False)

    best_state = copy.deepcopy(model.state_dict())
    best_rmse = math.inf
    best_epoch = 0
    stale = 0
    epochs = int(params["epochs"])
    patience = int(params["patience"])

    for epoch in range(1, epochs + 1):
        model.train()
        for xb, yb in loader:
            xb = xb.to(device)
            yb = yb.to(device)
            optimizer.zero_grad(set_to_none=True)
            pred = model(xb)
            loss = loss_fn(pred, yb)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()

        model.eval()
        with torch.no_grad():
            valid_scaled = model(valid_x).detach().cpu().numpy().astype(np.float32)
        valid_pred = valid_scaled * y_std + y_mean
        rmse = regression_metrics(y_valid_raw, valid_pred)["rmse"]
        if rmse + 1e-5 < best_rmse:
            best_rmse = float(rmse)
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                break

    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        train_pred = (model(torch.from_numpy(X_train).to(device)).detach().cpu().numpy().astype(np.float32) * y_std) + y_mean
        valid_pred = (model(valid_x).detach().cpu().numpy().astype(np.float32) * y_std) + y_mean

    fold_result = FoldResult(
        fold=0,
        best_epoch=int(best_epoch),
        best_valid_rmse=float(best_rmse),
        train_metrics=regression_metrics(y_train_raw, train_pred),
        valid_metrics=regression_metrics(y_valid_raw, valid_pred),
        elapsed_sec=float(time.time() - start),
    )
    return valid_pred.astype(np.float32), fold_result


def train_dl_models(
    cfg: dict[str, Any],
    input_name: str,
    model_names: list[str],
    device: torch.device,
    max_epochs: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = pipeline_paths(cfg)
    X = np.load(paths.processed_dir / f"X_{input_name}.npy").astype(np.float32)
    y = np.load(paths.processed_dir / "y_train.npy").astype(np.float32)
    feature_names = _json(paths.processed_dir / f"{input_name}_feature_names.json")
    blocks = feature_blocks(feature_names)
    seed = int(cfg["training"].get("random_state", cfg["project"].get("seed", 42)))
    n_splits = min(int(cfg["training"].get("n_splits", 3)), len(y))

    out_root = paths.results_dir / "additional_dl" / input_name
    oof_dir = out_root / "oof"
    fold_dir = out_root / "folds"
    oof_dir.mkdir(parents=True, exist_ok=True)
    fold_dir.mkdir(parents=True, exist_ok=True)

    splits = list(KFold(n_splits=n_splits, shuffle=True, random_state=seed).split(X))
    summary_rows = []
    fold_rows = []
    qc = {
        "step": "additional_dl_training",
        "input_set": input_name,
        "device": str(device),
        "shape": [int(X.shape[0]), int(X.shape[1])],
        "feature_blocks": {k: len(v) for k, v in blocks.items()},
        "models": {},
    }

    for model_name in model_names:
        params = dict(MODEL_DEFAULTS[model_name])
        if max_epochs is not None:
            params["epochs"] = int(max_epochs)
        print(f"[dl] {model_name} start params={params}", flush=True)
        model_start = time.time()
        oof = np.zeros(len(y), dtype=np.float32)
        model_fold_results = []
        for fold, (train_idx, valid_idx) in enumerate(splits, start=1):
            pred, fold_result = train_one_fold(
                model_name,
                X,
                y,
                train_idx,
                valid_idx,
                blocks,
                device,
                seed + fold * 100 + len(model_name),
                params,
            )
            fold_result.fold = fold
            oof[valid_idx] = pred
            row = {
                "input_set": input_name,
                "model": model_name,
                "fold": fold,
                "best_epoch": fold_result.best_epoch,
                "best_valid_rmse": fold_result.best_valid_rmse,
                "elapsed_sec": fold_result.elapsed_sec,
                "train_spearman": fold_result.train_metrics["spearman"],
                "valid_spearman": fold_result.valid_metrics["spearman"],
                "valid_rmse": fold_result.valid_metrics["rmse"],
                "valid_r2": fold_result.valid_metrics["r2"],
            }
            fold_rows.append(row)
            model_fold_results.append(
                {
                    "fold": fold,
                    "best_epoch": fold_result.best_epoch,
                    "elapsed_sec": fold_result.elapsed_sec,
                    "train": fold_result.train_metrics,
                    "valid": fold_result.valid_metrics,
                }
            )
            print(
                f"[dl] {model_name} fold={fold} epoch={fold_result.best_epoch} "
                f"spearman={fold_result.valid_metrics['spearman']:.4f} rmse={fold_result.valid_metrics['rmse']:.4f}",
                flush=True,
            )

        metrics = regression_metrics(y, oof)
        fold_train_spearman = [r["train"]["spearman"] for r in model_fold_results]
        train_oof_gap = float(np.nanmean(fold_train_spearman) - metrics["spearman"])
        elapsed = float(time.time() - model_start)
        np.save(oof_dir / f"{model_name}.npy", oof)
        payload = {
            "input_set": input_name,
            "model": model_name,
            "params": params,
            "fold_metrics": model_fold_results,
            "oof_metrics": metrics,
            "train_oof_spearman_gap": train_oof_gap,
            "prediction_variance": float(np.nanvar(oof)),
            "elapsed_sec": elapsed,
        }
        write_json(fold_dir / f"{model_name}_random3.json", payload)
        summary_rows.append(
            {
                "input_set": input_name,
                "model": model_name,
                **metrics,
                "train_oof_spearman_gap": train_oof_gap,
                "prediction_variance": float(np.nanvar(oof)),
                "elapsed_sec": elapsed,
            }
        )
        qc["models"][model_name] = {
            "elapsed_sec": elapsed,
            "oof_metrics": metrics,
            "best_epochs": [r["best_epoch"] for r in model_fold_results],
        }
        print(f"[dl] {model_name} done spearman={metrics['spearman']:.4f} rmse={metrics['rmse']:.4f}", flush=True)

    summary = pd.DataFrame(summary_rows).sort_values("spearman", ascending=False, na_position="last")
    fold_table = pd.DataFrame(fold_rows)
    write_table(summary, out_root / "additional_dl_metrics_summary.csv")
    write_table(fold_table, out_root / "additional_dl_fold_metrics.csv")
    write_json(paths.reports_dir / "qc_additional_dl_training_20260421.json", qc)
    return summary, fold_table


def train_dl_groupcv(
    cfg: dict[str, Any],
    input_name: str,
    model_names: list[str],
    device: torch.device,
    max_epochs: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = pipeline_paths(cfg)
    X = np.load(paths.processed_dir / f"X_{input_name}.npy").astype(np.float32)
    y = np.load(paths.processed_dir / "y_train.npy").astype(np.float32)
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    if rows is None or "canonical_drug_id" not in rows.columns:
        raise FileNotFoundError("row_metadata.csv with canonical_drug_id is required for GroupCV")
    groups = rows["canonical_drug_id"].astype(str).to_numpy()
    feature_names = _json(paths.processed_dir / f"{input_name}_feature_names.json")
    blocks = feature_blocks(feature_names)
    seed = int(cfg["training"].get("random_state", cfg["project"].get("seed", 42)))
    n_splits = min(3, len(np.unique(groups)))
    if n_splits < 2:
        raise ValueError("At least two drug groups are required for GroupCV")

    out_root = paths.results_dir / "additional_dl" / input_name / "groupcv"
    oof_dir = out_root / "oof"
    fold_dir = out_root / "folds"
    oof_dir.mkdir(parents=True, exist_ok=True)
    fold_dir.mkdir(parents=True, exist_ok=True)

    splits = list(GroupKFold(n_splits=n_splits).split(X, y, groups))
    summary_rows = []
    fold_rows = []
    qc = {
        "step": "additional_dl_groupcv",
        "input_set": input_name,
        "group_column": "canonical_drug_id",
        "device": str(device),
        "shape": [int(X.shape[0]), int(X.shape[1])],
        "unique_groups": int(len(np.unique(groups))),
        "feature_blocks": {k: len(v) for k, v in blocks.items()},
        "models": {},
    }

    for model_name in model_names:
        params = dict(MODEL_DEFAULTS[model_name])
        if max_epochs is not None:
            params["epochs"] = int(max_epochs)
        print(f"[dl-groupcv] {model_name} start params={params}", flush=True)
        model_start = time.time()
        oof = np.zeros(len(y), dtype=np.float32)
        model_fold_results = []
        for fold, (train_idx, valid_idx) in enumerate(splits, start=1):
            pred, fold_result = train_one_fold(
                model_name,
                X,
                y,
                train_idx,
                valid_idx,
                blocks,
                device,
                seed + fold * 1000 + len(model_name),
                params,
            )
            fold_result.fold = fold
            oof[valid_idx] = pred
            valid_groups = int(pd.Series(groups[valid_idx]).nunique())
            row = {
                "input_set": input_name,
                "model": model_name,
                "cv": "groupcv_by_drug",
                "fold": fold,
                "valid_drug_groups": valid_groups,
                "best_epoch": fold_result.best_epoch,
                "best_valid_rmse": fold_result.best_valid_rmse,
                "elapsed_sec": fold_result.elapsed_sec,
                "train_spearman": fold_result.train_metrics["spearman"],
                "valid_spearman": fold_result.valid_metrics["spearman"],
                "valid_rmse": fold_result.valid_metrics["rmse"],
                "valid_r2": fold_result.valid_metrics["r2"],
            }
            fold_rows.append(row)
            model_fold_results.append(
                {
                    "fold": fold,
                    "valid_drug_groups": valid_groups,
                    "best_epoch": fold_result.best_epoch,
                    "elapsed_sec": fold_result.elapsed_sec,
                    "train": fold_result.train_metrics,
                    "valid": fold_result.valid_metrics,
                }
            )
            print(
                f"[dl-groupcv] {model_name} fold={fold} groups={valid_groups} epoch={fold_result.best_epoch} "
                f"spearman={fold_result.valid_metrics['spearman']:.4f} rmse={fold_result.valid_metrics['rmse']:.4f}",
                flush=True,
            )

        metrics = regression_metrics(y, oof)
        fold_train_spearman = [r["train"]["spearman"] for r in model_fold_results]
        train_oof_gap = float(np.nanmean(fold_train_spearman) - metrics["spearman"])
        elapsed = float(time.time() - model_start)
        np.save(oof_dir / f"{model_name}.npy", oof)
        payload = {
            "input_set": input_name,
            "model": model_name,
            "cv": "groupcv_by_drug",
            "params": params,
            "fold_metrics": model_fold_results,
            "oof_metrics": metrics,
            "train_oof_spearman_gap": train_oof_gap,
            "prediction_variance": float(np.nanvar(oof)),
            "elapsed_sec": elapsed,
        }
        write_json(fold_dir / f"{model_name}_groupcv.json", payload)
        summary_rows.append(
            {
                "input_set": input_name,
                "model": model_name,
                "cv": "groupcv_by_drug",
                **metrics,
                "train_oof_spearman_gap": train_oof_gap,
                "prediction_variance": float(np.nanvar(oof)),
                "elapsed_sec": elapsed,
            }
        )
        qc["models"][model_name] = {
            "elapsed_sec": elapsed,
            "oof_metrics": metrics,
            "best_epochs": [r["best_epoch"] for r in model_fold_results],
        }
        print(f"[dl-groupcv] {model_name} done spearman={metrics['spearman']:.4f} rmse={metrics['rmse']:.4f}", flush=True)

    summary = pd.DataFrame(summary_rows).sort_values("spearman", ascending=False, na_position="last")
    fold_table = pd.DataFrame(fold_rows)
    write_table(summary, out_root / "additional_dl_groupcv_metrics_summary.csv")
    write_table(fold_table, out_root / "additional_dl_groupcv_fold_metrics.csv")
    write_json(paths.reports_dir / "qc_additional_dl_groupcv_20260421.json", qc)
    return summary, fold_table


def build_dl_ensemble(cfg: dict[str, Any], input_name: str, include_flatmlp: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    paths = pipeline_paths(cfg)
    y = np.load(paths.processed_dir / "y_train.npy").astype(np.float32)
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    out_root = paths.results_dir / "additional_dl" / input_name
    add_metrics = read_table(out_root / "additional_dl_metrics_summary.csv")
    if add_metrics is None or add_metrics.empty:
        raise FileNotFoundError(out_root / "additional_dl_metrics_summary.csv")

    metric_frames = [add_metrics.copy()]
    preds: dict[str, np.ndarray] = {}
    for row in add_metrics.itertuples():
        path = out_root / "oof" / f"{row.model}.npy"
        if path.exists():
            preds[str(row.model)] = np.load(path)

    if include_flatmlp:
        flat_path = paths.results_dir / "pancancer_lincs" / "oof" / input_name / "FlatMLP.npy"
        flat_metrics_path = paths.results_dir / "pancancer_lincs" / f"{input_name}_metrics_summary.csv"
        flat_metrics = read_table(flat_metrics_path)
        if flat_path.exists() and flat_metrics is not None and not flat_metrics.empty:
            flat_row = flat_metrics.loc[flat_metrics["model"].eq("FlatMLP")].copy()
            if not flat_row.empty:
                preds["FlatMLP_existing"] = np.load(flat_path)
                flat_row["model"] = "FlatMLP_existing"
                metric_frames.append(flat_row)

    metrics = pd.concat(metric_frames, ignore_index=True)
    weights = {}
    for row in metrics.itertuples():
        model_name = str(row.model)
        if model_name in preds:
            score = float(getattr(row, "spearman", np.nan))
            weights[model_name] = max(0.0, score) if np.isfinite(score) else 0.0
    total = sum(weights.values())
    if total <= 0:
        weights = {k: 1.0 / len(preds) for k in preds}
    else:
        weights = {k: v / total for k, v in weights.items()}

    ensemble_pred = np.zeros(len(y), dtype=np.float32)
    for model_name, pred in preds.items():
        ensemble_pred += float(weights[model_name]) * pred.astype(np.float32)

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
    cand = rows.copy()
    cand["dl_ensemble_pred_ln_ic50"] = ensemble_pred
    cand["dl_ensemble_score"] = -ensemble_pred
    group_cols = ["canonical_drug_id", "drug_name", "canonical_smiles", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"]
    available = [c for c in group_cols if c in cand.columns]
    top = (
        cand.groupby(available, dropna=False)
        .agg(
            mean_pred_ln_ic50=("dl_ensemble_pred_ln_ic50", "mean"),
            dl_ensemble_score=("dl_ensemble_score", "mean"),
            screened_rows=("dl_ensemble_score", "size"),
        )
        .reset_index()
        .sort_values("dl_ensemble_score", ascending=False)
    )
    if "drug_name" in top.columns:
        top["_drug_name_norm"] = top["drug_name"].map(lambda x: NAME_NORMALIZER.sub("", str(x).lower()))
        top = top.drop_duplicates("_drug_name_norm", keep="first").drop(columns=["_drug_name_norm"])
    top["rank"] = np.arange(1, len(top) + 1)
    top = top[["rank"] + [c for c in top.columns if c != "rank"]]

    ensemble_dir = out_root / "ensemble"
    ensemble_dir.mkdir(parents=True, exist_ok=True)
    np.save(ensemble_dir / "additional_dl_weighted_ensemble_oof.npy", ensemble_pred)
    write_table(top.head(30), ensemble_dir / "additional_dl_ensemble_top30_drugs.csv")
    write_table(diversity, ensemble_dir / "additional_dl_ensemble_diversity.csv")
    result = {
        "input_set": input_name,
        "weights": weights,
        "ensemble_metrics": regression_metrics(y, ensemble_pred),
        "member_metrics": metrics.to_dict(orient="records"),
        "diversity_rows": int(len(diversity)),
    }
    write_json(ensemble_dir / "additional_dl_ensemble_results.json", result)
    write_json(paths.reports_dir / "qc_additional_dl_ensemble_20260421.json", result)
    return top, diversity, result


def _fmt(value: Any) -> str:
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.4f}"
    except Exception:
        return str(value)


def _md_table(df: pd.DataFrame, cols: list[str]) -> str:
    subset = df[[c for c in cols if c in df.columns]].copy()
    if subset.empty:
        return "_No rows._"
    lines = ["| " + " | ".join(subset.columns) + " |", "| " + " | ".join(["---"] * len(subset.columns)) + " |"]
    for row in subset.itertuples(index=False):
        lines.append("| " + " | ".join(_fmt(v) for v in row) + " |")
    return "\n".join(lines)


def write_report(cfg: dict[str, Any], input_name: str, summary: pd.DataFrame, top: pd.DataFrame, diversity: pd.DataFrame, ensemble: dict[str, Any]) -> Path:
    paths = pipeline_paths(cfg)
    report = paths.root / "docs" / "THYROID_ADDITIONAL_DL_MODELS_20260421.md"
    metrics_cols = ["input_set", "model", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20", "train_oof_spearman_gap", "elapsed_sec"]
    top_cols = ["rank", "drug_name", "dl_ensemble_score", "mean_pred_ln_ic50", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"]
    diversity_cols = ["model_a", "model_b", "prediction_spearman_corr", "residual_pearson_corr", "mean_abs_prediction_gap"]
    text = f"""# Thyroid Additional DL Models - 2026-04-21

## 목적

기존 thyroid 파이프라인은 ML tree 계열과 sklearn `FlatMLP` 중심이었다. BRCA 파이프라인에서 사용했던 DL 계열을 thyroid primary input set에 추가 적용해, DL 단독 성능과 DL-only ensemble/diversity를 확인했다.

## 실행 설정

- Input set: `{input_name}`
- CV: random sample 3-fold OOF
- 추가 학습 모델: `ResidualMLP`, `WideDeep`, `CrossAttention`
- 비교용 포함: 기존 pan-cancer LINCS `FlatMLP_existing` OOF
- Device: PyTorch auto device

## 추가 DL 모델 성능

{_md_table(summary, metrics_cols)}

## DL-only Weighted Ensemble

- Ensemble Spearman: `{_fmt(ensemble.get("ensemble_metrics", {}).get("spearman"))}`
- Ensemble RMSE: `{_fmt(ensemble.get("ensemble_metrics", {}).get("rmse"))}`
- Ensemble R2: `{_fmt(ensemble.get("ensemble_metrics", {}).get("r2"))}`
- Weights: `{json.dumps(ensemble.get("weights", {}), ensure_ascii=False)}`

## DL Ensemble Top 30 중 상위 15

{_md_table(top.head(15), top_cols)}

## DL Diversity

{_md_table(diversity, diversity_cols)}

## 산출물

- `results/additional_dl/{input_name}/additional_dl_metrics_summary.csv`
- `results/additional_dl/{input_name}/additional_dl_fold_metrics.csv`
- `results/additional_dl/{input_name}/oof/`
- `results/additional_dl/{input_name}/ensemble/additional_dl_ensemble_results.json`
- `results/additional_dl/{input_name}/ensemble/additional_dl_ensemble_top30_drugs.csv`
- `results/additional_dl/{input_name}/ensemble/additional_dl_ensemble_diversity.csv`
"""
    report.write_text(text, encoding="utf-8")
    return report


def write_groupcv_report(cfg: dict[str, Any], input_name: str, summary: pd.DataFrame, fold_table: pd.DataFrame) -> Path:
    paths = pipeline_paths(cfg)
    report = paths.root / "docs" / "THYROID_ADDITIONAL_DL_GROUPCV_20260421.md"
    baseline_path = paths.results_dir / "pancancer_lincs" / "groupcv_stress_test" / f"{input_name}_ExtraTrees_groupcv.json"
    baseline = read_json(baseline_path) if baseline_path.exists() else {}
    metrics_cols = ["input_set", "model", "cv", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20", "train_oof_spearman_gap", "elapsed_sec"]
    fold_cols = ["model", "fold", "valid_drug_groups", "best_epoch", "valid_spearman", "valid_rmse", "valid_r2"]
    text = f"""# Thyroid Additional DL GroupCV - 2026-04-21

## 목적

random sample 3-fold에서 좋게 나온 DL 모델들이 약물 ID 기준 GroupCV에서도 유지되는지 확인했다. GroupCV는 `canonical_drug_id`를 group으로 두기 때문에, validation fold의 약물은 train fold에서 빠진다. 따라서 random split보다 훨씬 보수적인 unseen-drug stress test이다.

## 실행 설정

- Input set: `{input_name}`
- CV: GroupKFold 3-fold by `canonical_drug_id`
- Models: `ResidualMLP`, `WideDeep`, `CrossAttention`
- Device: PyTorch auto device

## DL GroupCV 성능

{_md_table(summary, metrics_cols)}

## Fold별 성능

{_md_table(fold_table, fold_cols)}

## 기존 ML GroupCV 기준

- Baseline model: `{baseline.get("model", "")}`
- Baseline Spearman: `{_fmt(baseline.get("oof_metrics", {}).get("spearman"))}`
- Baseline RMSE: `{_fmt(baseline.get("oof_metrics", {}).get("rmse"))}`
- Baseline R2: `{_fmt(baseline.get("oof_metrics", {}).get("r2"))}`

## 해석

GroupCV는 random sample 3-fold 성능보다 낮게 나오는 것이 정상이다. 여기서 성능이 유지되는 모델은 drug identity leakage에 덜 의존하고, unseen-drug 일반화에 상대적으로 강하다고 볼 수 있다.

## 산출물

- `results/additional_dl/{input_name}/groupcv/additional_dl_groupcv_metrics_summary.csv`
- `results/additional_dl/{input_name}/groupcv/additional_dl_groupcv_fold_metrics.csv`
- `results/additional_dl/{input_name}/groupcv/oof/`
- `reports/qc_additional_dl_groupcv_20260421.json`
"""
    report.write_text(text, encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Train additional PyTorch DL models for thyroid pipeline")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--input-set", default="numeric_strong_context_smiles_pan_lincs")
    parser.add_argument("--models", default="ResidualMLP,WideDeep,CrossAttention")
    parser.add_argument("--device", default="auto")
    parser.add_argument("--max-epochs", type=int, default=None)
    parser.add_argument("--no-flatmlp", action="store_true")
    parser.add_argument("--groupcv", action="store_true", help="Run GroupKFold by canonical_drug_id instead of random3 training")
    args = parser.parse_args()

    cfg = load_config(args.config)
    model_names = [m.strip() for m in args.models.split(",") if m.strip()]
    unknown = [m for m in model_names if m not in MODEL_DEFAULTS]
    if unknown:
        raise ValueError(f"Unknown model names: {unknown}")
    device = _torch_device(args.device)
    if args.groupcv:
        summary, fold_table = train_dl_groupcv(cfg, args.input_set, model_names, device, args.max_epochs)
        report = write_groupcv_report(cfg, args.input_set, summary, fold_table)
        print(
            json.dumps(
                {
                    "status": "completed",
                    "cv": "groupcv_by_drug",
                    "input_set": args.input_set,
                    "models": model_names,
                    "device": str(device),
                    "best_model": summary.iloc[0].to_dict() if not summary.empty else {},
                    "report": str(report),
                },
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0
    summary, _ = train_dl_models(cfg, args.input_set, model_names, device, args.max_epochs)
    top, diversity, ensemble = build_dl_ensemble(cfg, args.input_set, include_flatmlp=not args.no_flatmlp)
    report = write_report(cfg, args.input_set, summary, top, diversity, ensemble)
    print(
        json.dumps(
            {
                "status": "completed",
                "input_set": args.input_set,
                "models": model_names,
                "device": str(device),
                "best_model": summary.iloc[0].to_dict() if not summary.empty else {},
                "ensemble_metrics": ensemble.get("ensemble_metrics", {}),
                "report": str(report),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
