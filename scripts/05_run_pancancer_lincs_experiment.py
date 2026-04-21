from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold, KFold

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import (  # noqa: E402
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

try:
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
except Exception:  # pragma: no cover - optional dependency
    Chem = None


DEFAULT_GCTX = (
    "/Users/skku_aws2_18/drug/derived/official_dataset_sources_20260406/lincs/"
    "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz"
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


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)


def log(message: str) -> None:
    print(message, flush=True)


def norm_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def canonicalize_smiles(smiles: object) -> str:
    raw = str(smiles).strip()
    if not raw or raw.lower() in {"nan", "none", "<na>", "", "restricted"} or Chem is None:
        return ""
    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def decode_hdf5_strings(values: object) -> list[str]:
    arr = np.asarray(values).reshape(-1)
    out = []
    for item in arr:
        if isinstance(item, (bytes, np.bytes_)):
            out.append(item.decode("utf-8", errors="ignore"))
        elif pd.isna(item):
            out.append("")
        else:
            out.append(str(item))
    return out


def ensure_decompressed_gctx(gctx_gz: Path, cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / gctx_gz.name.removesuffix(".gz")
    if out_path.exists() and out_path.stat().st_size > 0:
        log(f"[pan-lincs] Reusing decompressed GCTX: {out_path}")
        return out_path
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    log(f"[pan-lincs] Decompressing GCTX to {out_path}")
    with gzip.open(gctx_gz, "rb") as src, tmp_path.open("wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024 * 16)
    tmp_path.replace(out_path)
    return out_path


def build_signature_mapping(paths_processed: Path, staging_lincs: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    drug = pd.read_csv(paths_processed / "raw" / "drug_features.csv")[
        ["canonical_drug_id", "drug_name", "canonical_smiles"]
    ].copy()
    drug["canonical_drug_id"] = drug["canonical_drug_id"].astype(str)

    name_to_ids: dict[str, set[str]] = defaultdict(set)
    smiles_to_ids: dict[str, set[str]] = defaultdict(set)
    for row in drug.itertuples(index=False):
        key = norm_name(row.drug_name)
        if key:
            name_to_ids[key].add(str(row.canonical_drug_id))
        smi = canonicalize_smiles(row.canonical_smiles)
        if smi:
            smiles_to_ids[smi].add(str(row.canonical_drug_id))

    cell = pd.read_parquet(staging_lincs / "lincs_cell_info_basic_20260406.parquet")
    sig = pd.read_parquet(
        staging_lincs / "lincs_sig_info_basic_20260406.parquet",
        columns=["sig_id", "pert_id", "pert_iname", "pert_type", "cell_id"],
    )
    pert = pd.read_parquet(
        staging_lincs / "lincs_pert_info_basic_20260406.parquet",
        columns=["pert_id", "pert_iname", "canonical_smiles"],
    )

    tumor_cells = set(
        cell.loc[
            cell["sample_type"].astype(str).str.contains("tumor", case=False, na=False),
            "cell_id",
        ].astype(str)
    )
    sig = sig[(sig["pert_type"].eq("trt_cp")) & (sig["cell_id"].astype(str).isin(tumor_cells))].copy()
    sig = sig.merge(pert.drop_duplicates(["pert_id", "pert_iname"]), on=["pert_id", "pert_iname"], how="left")

    mapped = []
    for row in sig.itertuples(index=False):
        matched_ids = name_to_ids.get(norm_name(row.pert_iname), set())
        match_source = "name"
        if not matched_ids:
            smi = canonicalize_smiles(getattr(row, "canonical_smiles", ""))
            matched_ids = smiles_to_ids.get(smi, set()) if smi else set()
            match_source = "smiles"
        for drug_id in matched_ids:
            mapped.append(
                {
                    "sig_id": str(row.sig_id),
                    "canonical_drug_id": str(drug_id),
                    "cell_id": str(row.cell_id),
                    "match_source": match_source,
                }
            )

    mapping = pd.DataFrame(mapped).drop_duplicates().reset_index(drop=True)
    if mapping.empty:
        raise RuntimeError("No pan-cancer LINCS signatures mapped to thyroid drug IDs.")

    report = {
        "pan_cancer_definition": "LINCS GSE92742 trt_cp signatures from tumor sample_type cell lines",
        "tumor_cell_ids_total": int(len(tumor_cells)),
        "tumor_trt_cp_signature_rows": int(len(sig)),
        "mapped_signature_drug_pairs": int(len(mapping)),
        "mapped_unique_signatures": int(mapping["sig_id"].nunique()),
        "mapped_drugs": int(mapping["canonical_drug_id"].nunique()),
        "mapped_cells": int(mapping["cell_id"].nunique()),
        "mapping_source_breakdown": mapping["match_source"].value_counts().to_dict(),
        "top_mapped_cell_counts": mapping["cell_id"].value_counts().head(20).to_dict(),
        "unmapped_drugs": int(len(set(drug["canonical_drug_id"]) - set(mapping["canonical_drug_id"]))),
    }
    return mapping, report


def build_pan_lincs_signature(cfg: dict[str, Any], gctx_gz: Path, cache_dir: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    paths = pipeline_paths(cfg)
    staging_lincs = paths.processed_dir / "source_staging" / "lincs"
    numeric_names = _read_json(paths.processed_dir / "numeric_feature_names.json")
    lincs_gene_features = [
        name
        for name in numeric_names
        if str(name).startswith("drug__lincs__")
    ]
    target_symbols = [name.replace("drug__lincs__", "").upper() for name in lincs_gene_features]
    target_symbol_set = set(target_symbols)

    mapping, mapping_report = build_signature_mapping(paths.processed_dir, staging_lincs)
    gene_info = pd.read_parquet(
        staging_lincs / "lincs_gene_info_basic_20260406.parquet",
        columns=["pr_gene_id", "pr_gene_symbol"],
    )
    gene_to_symbol = dict(zip(gene_info["pr_gene_id"].astype(str), gene_info["pr_gene_symbol"].astype(str).str.upper()))

    gctx_path = ensure_decompressed_gctx(gctx_gz, cache_dir)
    log("[pan-lincs] Reading selected GCTX rows/columns")
    with h5py.File(gctx_path, "r") as handle:
        row_ids = decode_hdf5_strings(handle["0"]["META"]["ROW"]["id"])
        col_ids = decode_hdf5_strings(handle["0"]["META"]["COL"]["id"])
        matrix = handle["0"]["DATA"]["0"]["matrix"]

        if matrix.shape == (len(row_ids), len(col_ids)):
            row_major = True
        elif matrix.shape == (len(col_ids), len(row_ids)):
            row_major = False
        else:
            raise RuntimeError(f"Unexpected GCTX matrix shape: {matrix.shape}")

        row_symbol_to_idx: dict[str, int] = {}
        for idx, row_id in enumerate(row_ids):
            symbol = gene_to_symbol.get(str(row_id).strip(), "")
            if symbol in target_symbol_set and symbol not in row_symbol_to_idx:
                row_symbol_to_idx[symbol] = idx

        missing_symbols = [symbol for symbol in target_symbols if symbol not in row_symbol_to_idx]
        kept_symbols = [symbol for symbol in target_symbols if symbol in row_symbol_to_idx]
        row_idx = np.asarray([row_symbol_to_idx[symbol] for symbol in kept_symbols], dtype=np.int64)

        col_lookup = {sig_id: idx for idx, sig_id in enumerate(col_ids)}
        selected = mapping[mapping["sig_id"].isin(col_lookup)].copy()
        if selected.empty:
            raise RuntimeError("Mapped pan-cancer signatures are absent from GCTX column metadata.")
        selected["col_idx"] = selected["sig_id"].map(col_lookup).astype(int)
        selected = selected.sort_values(["col_idx", "canonical_drug_id"]).reset_index(drop=True)
        unique_sig = selected[["sig_id", "col_idx"]].drop_duplicates().sort_values("col_idx").reset_index(drop=True)
        unique_col_idx = unique_sig["col_idx"].to_numpy(dtype=np.int64)

        if row_major:
            # Read all genes for only the mapped signatures, then retain the target LINCS feature genes.
            unique_subset = np.asarray(matrix[:, unique_col_idx], dtype=np.float32)[row_idx, :]
        else:
            unique_subset = np.asarray(matrix[unique_col_idx, :], dtype=np.float32)[:, row_idx].T

    sig_position = {sig_id: idx for idx, sig_id in enumerate(unique_sig["sig_id"].astype(str))}
    col_expand = np.asarray([sig_position[sig_id] for sig_id in selected["sig_id"].astype(str)], dtype=np.int64)
    selected_drug_ids = selected["canonical_drug_id"].astype(str).to_numpy()
    expanded = unique_subset[:, col_expand]

    drug_ids = sorted(
        selected["canonical_drug_id"].astype(str).unique(),
        key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value),
    )
    drug_rows = []
    for drug_id in drug_ids:
        drug_rows.append(expanded[:, selected_drug_ids == drug_id].mean(axis=1, dtype=np.float32))
    drug_mat = np.vstack(drug_rows).astype(np.float32)

    symbol_to_values = {symbol: drug_mat[:, i] for i, symbol in enumerate(kept_symbols)}
    full_columns: dict[str, Any] = {"canonical_drug_id": drug_ids}
    for feature_name, symbol in zip(lincs_gene_features, target_symbols):
        full_columns[feature_name] = symbol_to_values.get(symbol, np.zeros(len(drug_ids), dtype=np.float32))
    full = pd.DataFrame(full_columns)

    out_path = staging_lincs / "lincs_drug_signature_pancancer_20260421.parquet"
    full.to_parquet(out_path, index=False)
    report = {
        **mapping_report,
        "gctx_source": str(gctx_gz),
        "decompressed_gctx": str(gctx_path),
        "row_major": row_major,
        "target_lincs_gene_features": int(len(lincs_gene_features)),
        "matched_gctx_gene_features": int(len(kept_symbols)),
        "missing_gctx_gene_features": int(len(missing_symbols)),
        "mapped_signatures_in_gctx": int(selected["sig_id"].nunique()),
        "mapped_signature_drug_pairs_in_gctx": int(len(selected)),
        "pan_lincs_drug_signature_path": str(out_path),
        "pan_lincs_drug_signature_shape": [int(full.shape[0]), int(full.shape[1])],
    }
    write_json(paths.reports_dir / "qc_pancancer_lincs_signature_20260421.json", report)
    return full, report


def build_pan_lincs_inputs(cfg: dict[str, Any], signature: pd.DataFrame) -> tuple[dict[str, Path], pd.DataFrame]:
    paths = pipeline_paths(cfg)
    rows_meta = pd.read_csv(paths.processed_dir / "row_metadata.csv", usecols=["canonical_drug_id"])
    row_drug_ids = rows_meta["canonical_drug_id"].astype(str).to_numpy()
    signature = signature.copy()
    signature["canonical_drug_id"] = signature["canonical_drug_id"].astype(str)
    sig_by_drug = signature.set_index("canonical_drug_id")

    input_paths: dict[str, Path] = {}
    summary_rows = []
    matched_drugs = set(sig_by_drug.index)

    for base_name, (array_name, feature_name_file) in INPUT_DEFS.items():
        X = np.load(paths.processed_dir / array_name)
        names = [str(x) for x in _read_json(paths.processed_dir / feature_name_file)]
        X_pan = X.copy()
        lincs_positions = [i for i, name in enumerate(names) if "lincs" in name.lower()]
        gene_positions = [(i, name) for i, name in enumerate(names) if name.startswith("drug__lincs__")]
        has_idx = names.index("drug__has_lincs_signature")
        direct_idx = names.index("drug__lincs_signature_source_direct")
        recovered_idx = names.index("drug__lincs_signature_source_recovered")

        X_pan[:, lincs_positions] = 0.0
        has_values = np.isin(row_drug_ids, list(matched_drugs)).astype(np.float32)
        X_pan[:, has_idx] = has_values
        X_pan[:, direct_idx] = has_values
        X_pan[:, recovered_idx] = 0.0

        for pos, feature_name in gene_positions:
            if feature_name not in sig_by_drug.columns:
                continue
            value_map = sig_by_drug[feature_name]
            X_pan[:, pos] = pd.Series(row_drug_ids).map(value_map).fillna(0.0).to_numpy(dtype=np.float32)

        out_name = f"{base_name}_pan_lincs"
        out_path = paths.processed_dir / f"X_{out_name}.npy"
        names_path = paths.processed_dir / f"{out_name}_feature_names.json"
        np.save(out_path, X_pan.astype(np.float32))
        _write_json(names_path, names)
        input_paths[out_name] = out_path
        summary_rows.append(
            {
                "input_set": out_name,
                "source_input_set": base_name,
                "rows": int(X.shape[0]),
                "features": int(X.shape[1]),
                "replaced_lincs_features": int(len(lincs_positions)),
                "pan_lincs_drugs": int(len(matched_drugs)),
                "row_level_has_pan_lincs_ratio": float(has_values.mean()),
                "array_path": str(out_path),
                "feature_names_path": str(names_path),
            }
        )

    summary = pd.DataFrame(summary_rows)
    write_table(summary, paths.reports_dir / "pancancer_lincs_input_summary_20260421.csv")
    write_json(paths.reports_dir / "pancancer_lincs_input_summary_20260421.json", summary.to_dict(orient="records"))
    return input_paths, summary


def _fold_coverage(rows: pd.DataFrame, splits: list[tuple[np.ndarray, np.ndarray]], col: str) -> list[dict[str, int]]:
    return [{"fold": i, f"unique_{col}": int(rows.iloc[valid_idx][col].nunique())} for i, (_, valid_idx) in enumerate(splits, start=1)]


def train_pan_lincs(cfg: dict[str, Any], input_paths: dict[str, Path]) -> dict[str, pd.DataFrame]:
    paths = pipeline_paths(cfg)
    y = np.load(paths.processed_dir / "y_train.npy")
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    seed = int(cfg["training"].get("random_state", cfg["project"].get("seed", 42)))
    n_splits = min(int(cfg["training"].get("n_splits", 3)), len(y))
    models = list(cfg["training"]["models"])
    out_root = paths.results_dir / "pancancer_lincs"
    summaries: dict[str, pd.DataFrame] = {}
    qc: dict[str, Any] = {"step": "pancancer_lincs_random3_training", "input_sets": {}, "skipped_models": []}

    for input_name, input_path in input_paths.items():
        X = np.load(input_path)
        out_dir = out_root / "random3" / input_name
        oof_dir = out_root / "oof" / input_name
        out_dir.mkdir(parents=True, exist_ok=True)
        oof_dir.mkdir(parents=True, exist_ok=True)
        splits = list(KFold(n_splits=n_splits, shuffle=True, random_state=seed).split(X))
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
                "lincs_variant": "pan_cancer_tumor_cells",
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

    write_json(paths.reports_dir / "qc_pancancer_lincs_training_20260421.json", qc)
    return summaries


def run_groupcv_pan_lincs(cfg: dict[str, Any], input_name: str = "numeric_strong_context_smiles_pan_lincs") -> dict[str, Any]:
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

    out_dir = paths.results_dir / "pancancer_lincs" / "groupcv_stress_test"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"{input_name}_{model_name}_groupcv_oof.npy", oof)
    result = {
        "status": "completed",
        "input_set": input_name,
        "model": model_name,
        "lincs_variant": "pan_cancer_tumor_cells",
        "fold_metrics": fold_metrics,
        "oof_metrics": regression_metrics(y, oof),
    }
    write_json(out_dir / f"{input_name}_{model_name}_groupcv.json", result)
    return result


def run_ensemble_pan_lincs(cfg: dict[str, Any], input_name: str = "numeric_strong_context_smiles_pan_lincs") -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    paths = pipeline_paths(cfg)
    out_root = paths.results_dir / "pancancer_lincs"
    y = np.load(paths.processed_dir / "y_train.npy")
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    metrics = read_table(out_root / f"{input_name}_metrics_summary.csv")
    if metrics is None or metrics.empty:
        raise FileNotFoundError(f"No pan-LINCS metrics found for {input_name}")

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
        raise ValueError("No pan-LINCS OOF predictions found")
    weight_sum = sum(weights.values())
    weights = {k: (v / weight_sum if weight_sum > 0 else 1.0 / len(preds)) for k, v in weights.items()}

    ensemble_pred = np.zeros(len(y), dtype=np.float32)
    for model_name, pred in preds.items():
        ensemble_pred += weights[model_name] * pred

    diversity_rows = []
    model_names = list(preds)
    for i, left in enumerate(model_names):
        for right in model_names[i + 1 :]:
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
    write_table(top.head(30), out_dir / "thyroid_pancancer_lincs_ensemble_top30_drugs.csv")
    write_table(diversity, out_dir / "thyroid_pancancer_lincs_ensemble_diversity.csv")
    result = {
        "input_set": input_name,
        "lincs_variant": "pan_cancer_tumor_cells",
        "weights": weights,
        "ensemble_metrics": regression_metrics(y, ensemble_pred),
        "best_single_model": metrics.iloc[0].to_dict(),
        "diversity_rows": int(len(diversity)),
    }
    write_json(out_dir / "thyroid_pancancer_lincs_ensemble_results.json", result)
    write_json(paths.reports_dir / "qc_pancancer_lincs_ensemble_20260421.json", result)
    return top, diversity, result


def _read_result_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def compare_variants(cfg: dict[str, Any], input_summary: pd.DataFrame, groupcv: dict[str, Any], ensemble: dict[str, Any]) -> pd.DataFrame:
    paths = pipeline_paths(cfg)
    rows = []
    for pan_input in input_summary["input_set"].tolist():
        source_input = str(input_summary.loc[input_summary["input_set"].eq(pan_input), "source_input_set"].iloc[0])
        variants = [
            ("baseline_mcf7_lincs", paths.results_dir / f"{source_input}_metrics_summary.csv"),
            ("no_lincs", paths.results_dir / "no_lincs" / f"{source_input}_no_lincs_metrics_summary.csv"),
            ("pancancer_lincs", paths.results_dir / "pancancer_lincs" / f"{pan_input}_metrics_summary.csv"),
        ]
        for variant, path in variants:
            table = read_table(path)
            if table is None or table.empty:
                continue
            best = table.sort_values("spearman", ascending=False).iloc[0]
            rows.append(
                {
                    "source_input_set": source_input,
                    "variant": variant,
                    "model": best["model"],
                    "spearman": best["spearman"],
                    "rmse": best["rmse"],
                    "r2": best["r2"],
                    "ndcg_at_20": best.get("ndcg_at_20", np.nan),
                }
            )

    comparison = pd.DataFrame(rows)
    write_table(comparison, paths.results_dir / "pancancer_lincs" / "lincs_variant_best_model_comparison.csv")

    baseline_groupcv = _read_result_json(paths.results_dir / "groupcv_stress_test" / "numeric_strong_context_smiles_ExtraTrees_groupcv.json")
    no_groupcv = _read_result_json(paths.results_dir / "no_lincs" / "groupcv_stress_test" / "numeric_strong_context_smiles_no_lincs_ExtraTrees_groupcv.json")
    group_rows = []
    for variant, payload in [("baseline_mcf7_lincs", baseline_groupcv), ("no_lincs", no_groupcv), ("pancancer_lincs", groupcv)]:
        if not payload:
            continue
        group_rows.append(
            {
                "variant": variant,
                "input_set": payload.get("input_set"),
                "model": payload.get("model"),
                **payload.get("oof_metrics", {}),
            }
        )
    group_comparison = pd.DataFrame(group_rows)
    write_table(group_comparison, paths.results_dir / "pancancer_lincs" / "lincs_variant_groupcv_comparison.csv")
    _write_markdown_report(cfg, input_summary, comparison, group_comparison, ensemble)
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
    signature_qc = _read_result_json(paths.reports_dir / "qc_pancancer_lincs_signature_20260421.json")
    metrics = ensemble.get("ensemble_metrics", {})
    best = ensemble.get("best_single_model", {})
    text = f"""# Thyroid Pan-Cancer LINCS Experiment - 2026-04-21

## 목적

MCF7-only LINCS feature의 유방암 세포주 편향을 줄이기 위해, LINCS GSE92742의 tumor sample-type cancer cell line 전체를 사용한 pan-cancer drug perturbation signature로 LINCS block을 교체했다.

## Pan-Cancer LINCS 정의

- Source: `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz`
- Cell filter: `sample_type == tumor`
- Perturbation filter: `pert_type == trt_cp`
- Drug matching: thyroid screened drug의 normalized name 우선, canonical SMILES 보조
- Feature replacement: 기존 `drug__lincs__*` 1024개 gene feature와 LINCS availability flag 3개를 pan-cancer 값으로 교체

## Signature QC

- Mapped drugs: `{signature_qc.get("mapped_drugs", "")}`
- Mapped cells: `{signature_qc.get("mapped_cells", "")}`
- Mapped unique signatures: `{signature_qc.get("mapped_unique_signatures", "")}`
- Matched GCTX gene features: `{signature_qc.get("matched_gctx_gene_features", "")}`

## 입력셋

{_markdown_table(input_summary)}

## Random Sample 3-fold Best Model 비교

{_markdown_table(comparison)}

## GroupCV Stress Test 비교

{_markdown_table(groupcv_comparison)}

## Pan-Cancer LINCS Ensemble

- Primary input set: `numeric_strong_context_smiles_pan_lincs`
- Best single model: `{best.get("model", "")}`
- Ensemble Spearman: `{_fmt(metrics.get("spearman"))}`
- Ensemble RMSE: `{_fmt(metrics.get("rmse"))}`
- Ensemble R2: `{_fmt(metrics.get("r2"))}`

## 산출물

- `data/source_staging/lincs/lincs_drug_signature_pancancer_20260421.parquet`
- `data/X_numeric_pan_lincs.npy`
- `data/X_numeric_smiles_pan_lincs.npy`
- `data/X_numeric_strong_context_smiles_pan_lincs.npy`
- `results/pancancer_lincs/random3/`
- `results/pancancer_lincs/groupcv_stress_test/`
- `results/pancancer_lincs/ensemble/`
- `reports/qc_pancancer_lincs_signature_20260421.json`
- `reports/qc_pancancer_lincs_training_20260421.json`
"""
    (paths.root / "docs" / "THYROID_PANCANCER_LINCS_EXPERIMENT_20260421.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build pan-cancer LINCS thyroid inputs and run model ablation")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--gctx", default=DEFAULT_GCTX)
    parser.add_argument("--cache-dir", default="/tmp/thyroid_lincs_pancancer_cache")
    args = parser.parse_args()

    cfg = load_config(args.config)
    signature, signature_report = build_pan_lincs_signature(cfg, Path(args.gctx), Path(args.cache_dir))
    log(f"[pan-lincs] Signature completed: {signature_report['pan_lincs_drug_signature_shape']}")
    input_paths, input_summary = build_pan_lincs_inputs(cfg, signature)
    log("[pan-lincs] Pan-LINCS input arrays built")
    train_pan_lincs(cfg, input_paths)
    log("[pan-lincs] Random3 training completed")
    groupcv = run_groupcv_pan_lincs(cfg)
    log("[pan-lincs] GroupCV completed")
    _, _, ensemble = run_ensemble_pan_lincs(cfg)
    compare_variants(cfg, input_summary, groupcv, ensemble)
    print(json.dumps({"status": "completed", "input_sets": list(input_paths)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
