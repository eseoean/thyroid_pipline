#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

PAAD_SCRIPT_PATH = REPO_ROOT / "scripts" / "11_run_paad_full_pipeline.py"
PAN_INPUT = "numeric_strong_context_smiles_pan_lincs"
PAAD_FOCUSED_INPUT = "numeric_strong_context_smiles_paad_lincs"
NO_LINCS_INPUT = "numeric_strong_context_smiles_no_lincs"
YAPC_INPUT = "numeric_strong_context_smiles_yapc_lincs"


def _load_paad_module() -> Any:
    spec = importlib.util.spec_from_file_location("paad_full_pipeline", PAAD_SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load PAAD helper script: {PAAD_SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["paad_full_pipeline"] = module
    spec.loader.exec_module(module)
    return module


PIPE = _load_paad_module()


def log(message: str) -> None:
    print(message, flush=True)


def read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def norm_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def decode_hdf5_strings(values: object) -> list[str]:
    arr = np.asarray(values).reshape(-1)
    out: list[str] = []
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
        log(f"[yapc-lincs] Reusing decompressed GCTX: {out_path}")
        return out_path
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    log(f"[yapc-lincs] Decompressing GCTX to {out_path}")
    with gzip.open(gctx_gz, "rb") as src, tmp_path.open("wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024 * 16)
    tmp_path.replace(out_path)
    return out_path


def _fmt(value: Any) -> str:
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.4f}"
    except Exception:
        return str(value)


def _md_table(df: pd.DataFrame, columns: list[str], n: int = 20) -> str:
    if df is None or df.empty:
        return "_No rows._"
    sub = df[[c for c in columns if c in df.columns]].head(n).copy()
    if sub.empty:
        return "_No matching columns._"
    lines = ["| " + " | ".join(sub.columns) + " |", "| " + " | ".join(["---"] * len(sub.columns)) + " |"]
    for row in sub.itertuples(index=False):
        lines.append("| " + " | ".join(_fmt(v) for v in row) + " |")
    return "\n".join(lines)


def lincs_phase2_paths(root: Path, lincs_dir: str) -> dict[str, Path]:
    base = root / lincs_dir
    return {
        "base": base,
        "cell": base / "GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz",
        "sig": base / "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
        "inst": base / "GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz",
        "pert": base / "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
        "gene": base / "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
        "gctx_gz": base / "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz",
    }


def build_yapc_mapping(cfg: dict[str, Any], lincs_paths: dict[str, Path]) -> tuple[pd.DataFrame, dict[str, Any]]:
    paths = PIPE.pipeline_paths(cfg)
    drugs = pd.read_csv(paths.raw_dir / "drug_features.csv")[
        ["canonical_drug_id", "drug_name", "canonical_smiles"]
    ].drop_duplicates()
    drugs["canonical_drug_id"] = drugs["canonical_drug_id"].astype(str)

    name_to_ids: dict[str, set[str]] = {}
    smiles_to_ids: dict[str, set[str]] = {}
    for row in drugs.itertuples(index=False):
        key = norm_name(row.drug_name)
        if key:
            name_to_ids.setdefault(key, set()).add(str(row.canonical_drug_id))
        can, ok = PIPE.SRC.canonicalize_smiles(row.canonical_smiles)
        if ok and can:
            smiles_to_ids.setdefault(can, set()).add(str(row.canonical_drug_id))

    cell = pd.read_csv(lincs_paths["cell"], sep="\t")
    pancreas_mask = pd.Series(False, index=cell.index)
    for col in cell.columns:
        pancreas_mask |= cell[col].astype(str).str.contains(
            "pancre|paad|panc|yapc",
            case=False,
            regex=True,
            na=False,
        )
    pancreas_cells = cell.loc[pancreas_mask].copy()
    pancreas_ids = set(pancreas_cells["cell_id"].astype(str))
    target_cells = sorted(c for c in pancreas_ids if c.upper().startswith("YAPC"))

    sig = pd.read_csv(
        lincs_paths["sig"],
        sep="\t",
        usecols=["sig_id", "pert_id", "pert_iname", "pert_type", "cell_id", "pert_idose", "pert_itime"],
    )
    sig = sig[sig["cell_id"].astype(str).isin(target_cells)].copy()
    trt = sig[sig["pert_type"].eq("trt_cp")].copy()
    pert = pd.read_csv(lincs_paths["pert"], sep="\t")
    trt = trt.merge(
        pert[["pert_id", "pert_iname", "canonical_smiles", "inchi_key"]].drop_duplicates(["pert_id", "pert_iname"]),
        on=["pert_id", "pert_iname"],
        how="left",
    )

    records = []
    for row in trt.itertuples(index=False):
        matched_ids = name_to_ids.get(norm_name(row.pert_iname), set())
        match_source = "name"
        if not matched_ids:
            can, ok = PIPE.SRC.canonicalize_smiles(getattr(row, "canonical_smiles", ""))
            matched_ids = smiles_to_ids.get(can, set()) if ok and can else set()
            match_source = "smiles"
        for drug_id in matched_ids:
            records.append(
                {
                    "sig_id": str(row.sig_id),
                    "canonical_drug_id": str(drug_id),
                    "pert_id": str(row.pert_id),
                    "pert_iname": str(row.pert_iname),
                    "cell_id": str(row.cell_id),
                    "pert_idose": str(row.pert_idose),
                    "pert_itime": str(row.pert_itime),
                    "match_source": match_source,
                }
            )
    mapping = pd.DataFrame(records).drop_duplicates().reset_index(drop=True)
    if mapping.empty:
        raise RuntimeError("No GSE70138 YAPC LINCS signatures mapped to PAAD drug IDs.")

    mapping_out = paths.reports_dir / "paad_yapc_lincs_phase2_mapping_20260422.csv"
    mapping_out.parent.mkdir(parents=True, exist_ok=True)
    mapping.merge(drugs[["canonical_drug_id", "drug_name"]], on="canonical_drug_id", how="left").to_csv(mapping_out, index=False)

    report = {
        "phase2_lincs_dir": str(lincs_paths["base"]),
        "pancreas_metadata_cell_ids": sorted(pancreas_ids),
        "target_cell_ids": target_cells,
        "target_cell_metadata_rows": pancreas_cells[pancreas_cells["cell_id"].astype(str).isin(target_cells)].to_dict(orient="records"),
        "target_cell_signature_rows": int(len(sig)),
        "target_cell_trt_cp_signature_rows": int(len(trt)),
        "target_cell_unique_trt_cp_perturbagens": int(trt["pert_iname"].nunique()),
        "mapped_signature_drug_pairs": int(len(mapping)),
        "mapped_unique_signatures": int(mapping["sig_id"].nunique()),
        "mapped_drugs": int(mapping["canonical_drug_id"].nunique()),
        "mapped_cells": mapping["cell_id"].value_counts().to_dict(),
        "mapping_source_breakdown": mapping["match_source"].value_counts().to_dict(),
        "mapping_path": str(mapping_out),
    }
    return mapping, report


def build_yapc_signature(
    cfg: dict[str, Any],
    lincs_paths: dict[str, Path],
    cache_dir: Path,
    feature_limit: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    paths = PIPE.pipeline_paths(cfg)
    mapping, mapping_report = build_yapc_mapping(cfg, lincs_paths)
    gene_info = pd.read_csv(lincs_paths["gene"], sep="\t", usecols=["pr_gene_id", "pr_gene_symbol"])
    gene_to_symbol = dict(zip(gene_info["pr_gene_id"].astype(str), gene_info["pr_gene_symbol"].astype(str).str.upper()))

    gctx_path = ensure_decompressed_gctx(lincs_paths["gctx_gz"], cache_dir)
    log("[yapc-lincs] Reading selected GCTX signatures")
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

        col_lookup = {sig_id: idx for idx, sig_id in enumerate(col_ids)}
        selected = mapping[mapping["sig_id"].isin(col_lookup)].copy()
        if selected.empty:
            raise RuntimeError("Mapped YAPC signatures are absent from GSE70138 Level5 column metadata.")
        selected["col_idx"] = selected["sig_id"].map(col_lookup).astype(int)
        selected = selected.sort_values(["col_idx", "canonical_drug_id"]).reset_index(drop=True)
        unique_sig = selected[["sig_id", "col_idx"]].drop_duplicates().sort_values("col_idx").reset_index(drop=True)
        unique_col_idx = unique_sig["col_idx"].to_numpy(dtype=np.int64)

        if row_major:
            gene_by_sig = np.asarray(matrix[:, unique_col_idx], dtype=np.float32)
        else:
            gene_by_sig = np.asarray(matrix[unique_col_idx, :], dtype=np.float32).T

    sig_position = {sig_id: idx for idx, sig_id in enumerate(unique_sig["sig_id"].astype(str))}
    col_expand = np.asarray([sig_position[sig_id] for sig_id in selected["sig_id"].astype(str)], dtype=np.int64)
    selected_drug_ids = selected["canonical_drug_id"].astype(str).to_numpy()
    expanded = gene_by_sig[:, col_expand]

    drug_ids = sorted(
        selected["canonical_drug_id"].astype(str).unique(),
        key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value),
    )
    drug_rows = []
    for drug_id in drug_ids:
        drug_rows.append(expanded[:, selected_drug_ids == drug_id].mean(axis=1, dtype=np.float32))
    drug_mat = np.vstack(drug_rows).astype(np.float32)

    gene_symbols = []
    seen: dict[str, int] = {}
    for row_id in row_ids:
        symbol = PIPE.SRC.clean_gene_symbol(gene_to_symbol.get(str(row_id).strip(), str(row_id).strip()))
        if not symbol:
            symbol = f"GENE_{len(gene_symbols)}"
        seen[symbol] = seen.get(symbol, 0) + 1
        gene_symbols.append(symbol if seen[symbol] == 1 else f"{symbol}_{seen[symbol]}")

    full = pd.DataFrame(drug_mat, columns=gene_symbols)
    full.insert(0, "canonical_drug_id", drug_ids)
    selected_features = PIPE.SRC.top_variance_columns(full, "canonical_drug_id", int(feature_limit), set())
    out = full[["canonical_drug_id"] + selected_features].copy()
    out = out.rename(columns={c: f"drug__yapc_lincs__{c}" for c in selected_features})

    signature_out = lincs_paths["base"] / "lincs_drug_signature_yapc_phase2_20260422.parquet"
    out.to_parquet(signature_out, index=False)
    feature_out = paths.reports_dir / "paad_yapc_lincs_selected_features_20260422.csv"
    pd.DataFrame({"feature_name": selected_features}).to_csv(feature_out, index=False)

    report = {
        **mapping_report,
        "gctx_source": str(lincs_paths["gctx_gz"]),
        "decompressed_gctx": str(gctx_path),
        "row_major": bool(row_major),
        "mapped_signatures_in_gctx": int(selected["sig_id"].nunique()),
        "mapped_signature_drug_pairs_in_gctx": int(len(selected)),
        "all_gctx_gene_features": int(len(gene_symbols)),
        "selected_yapc_lincs_features": int(len(selected_features)),
        "signature_path": str(signature_out),
        "selected_feature_path": str(feature_out),
        "signature_shape": [int(out.shape[0]), int(out.shape[1])],
    }
    PIPE.write_json(paths.reports_dir / "qc_paad_yapc_lincs_signature_20260422.json", report)
    return out, report


def build_yapc_input(cfg: dict[str, Any], signature: pd.DataFrame, output_input: str) -> dict[str, Any]:
    paths = PIPE.pipeline_paths(cfg)
    X_base = np.load(paths.processed_dir / f"X_{NO_LINCS_INPUT}.npy").astype(np.float32)
    base_names = [str(x) for x in read_json(paths.processed_dir / f"{NO_LINCS_INPUT}_feature_names.json")]
    rows = pd.read_csv(paths.processed_dir / "row_metadata.csv", usecols=["canonical_drug_id"])
    row_drug_ids = rows["canonical_drug_id"].astype(str)
    signature = signature.copy()
    signature["canonical_drug_id"] = signature["canonical_drug_id"].astype(str)
    sig_by_drug = signature.set_index("canonical_drug_id")

    feature_cols = [c for c in signature.columns if c != "canonical_drug_id"]
    mapped = row_drug_ids.isin(sig_by_drug.index).astype(np.float32).to_numpy()
    lincs_values = np.zeros((len(row_drug_ids), len(feature_cols)), dtype=np.float32)
    for i, col in enumerate(feature_cols):
        lincs_values[:, i] = row_drug_ids.map(sig_by_drug[col]).fillna(0.0).to_numpy(dtype=np.float32)

    X_out = np.concatenate([X_base, lincs_values, mapped.reshape(-1, 1)], axis=1).astype(np.float32)
    names_out = base_names + feature_cols + ["drug__has_yapc_lincs_signature"]
    np.save(paths.processed_dir / f"X_{output_input}.npy", X_out)
    write_json(paths.processed_dir / f"{output_input}_feature_names.json", names_out)

    report = {
        "step": "paad_yapc_lincs_input",
        "base_input": NO_LINCS_INPUT,
        "output_input": output_input,
        "base_shape": [int(X_base.shape[0]), int(X_base.shape[1])],
        "output_shape": [int(X_out.shape[0]), int(X_out.shape[1])],
        "yapc_lincs_feature_count": int(len(feature_cols)),
        "yapc_lincs_drugs": int(len(sig_by_drug.index)),
        "row_level_yapc_lincs_coverage": float(mapped.mean()),
        "matrix_path": str(paths.processed_dir / f"X_{output_input}.npy"),
        "feature_names_path": str(paths.processed_dir / f"{output_input}_feature_names.json"),
    }
    PIPE.write_json(paths.reports_dir / "qc_paad_yapc_lincs_input_20260422.json", report)
    return report


def run_training(cfg: dict[str, Any], input_name: str, force: bool) -> dict[str, Any]:
    cv_results: dict[str, Any] = {}
    for cv in ["random4", "groupcv4_drug"]:
        log(f"[yapc-lincs] train ML {cv} input={input_name}")
        ml_summary, _ = PIPE.train_ml_cv(cfg, input_name, cv, force=force)
        log(f"[yapc-lincs] train DL {cv} input={input_name}")
        dl_summary, _ = PIPE.train_dl_cv(cfg, input_name, cv, force=force)
        log(f"[yapc-lincs] build ensemble {cv} input={input_name}")
        ensemble = PIPE.build_cv_ensemble(cfg, input_name, cv)
        cv_results[cv] = {
            "ml_best": ml_summary.iloc[0].to_dict() if not ml_summary.empty else {},
            "dl_best": dl_summary.iloc[0].to_dict() if not dl_summary.empty else {},
            "ensemble_best": ensemble.get("ensemble_metrics", [{}])[0],
            "top50_path": ensemble.get("top50_path"),
        }
    return cv_results


def write_variant_report(cfg: dict[str, Any], input_name: str, run_summary: dict[str, Any]) -> Path:
    paths = PIPE.pipeline_paths(cfg)
    random_ens = pd.read_csv(paths.results_dir / "ensemble" / "random4" / input_name / "ensemble_metrics.csv")
    group_ens = pd.read_csv(paths.results_dir / "ensemble" / "groupcv4_drug" / input_name / "ensemble_metrics.csv")
    random_ind = pd.read_csv(paths.results_dir / "ensemble" / "random4" / input_name / "individual_metrics.csv")
    group_ind = pd.read_csv(paths.results_dir / "ensemble" / "groupcv4_drug" / input_name / "individual_metrics.csv")
    top = pd.read_csv(paths.results_dir / "ensemble" / "groupcv4_drug" / input_name / "ensemble_top50_drugs.csv")

    comparison_rows = []
    for variant, cv, candidate_input in [
        ("pan_cancer_lincs", "random4", PAN_INPUT),
        ("paad_focused_lincs", "random4", PAAD_FOCUSED_INPUT),
        ("yapc_cell_lincs", "random4", input_name),
        ("pan_cancer_lincs", "groupcv4_drug", PAN_INPUT),
        ("paad_focused_lincs", "groupcv4_drug", PAAD_FOCUSED_INPUT),
        ("yapc_cell_lincs", "groupcv4_drug", input_name),
    ]:
        path = paths.results_dir / "ensemble" / cv / candidate_input / "ensemble_metrics.csv"
        if path.exists():
            row = pd.read_csv(path).sort_values("spearman", ascending=False).iloc[0].to_dict()
            comparison_rows.append({"variant": variant, "cv": cv, **row})
    comparison = pd.DataFrame(comparison_rows)
    PIPE.write_table(comparison, paths.results_dir / "ensemble" / "paad_yapc_lincs_variant_comparison.csv")

    sig_qc = run_summary.get("signature_qc", {})
    input_qc = run_summary.get("input_qc", {})
    text = f"""# PAAD YAPC LINCS Variant - 2026-04-22

## Scope

- Cancer type: pancreatic cancer / TCGA-PAAD
- LINCS source: GSE70138 Phase II Level5
- Cell filter: `YAPC`, `YAPC.311`
- Perturbation filter: `pert_type == trt_cp`
- Input set: `{input_name}`
- Baseline comparator: `{PAN_INPUT}`
- CV checks: random sample 4-fold and drug GroupCV 4-fold
- Model families: ML, DL, and ML+DL ensembles

## Signature QC

- Target cells: `{sig_qc.get("target_cell_ids")}`
- YAPC/YAPC.311 signature rows: `{sig_qc.get("target_cell_signature_rows")}`
- YAPC/YAPC.311 trt_cp signature rows: `{sig_qc.get("target_cell_trt_cp_signature_rows")}`
- Unique YAPC perturbagens: `{sig_qc.get("target_cell_unique_trt_cp_perturbagens")}`
- Mapped PAAD GDSC drugs: `{sig_qc.get("mapped_drugs")}`
- Mapped Level5 signatures: `{sig_qc.get("mapped_signatures_in_gctx")}`
- Selected YAPC LINCS features: `{sig_qc.get("selected_yapc_lincs_features")}`
- Signature table: `{sig_qc.get("signature_path")}`

## Input QC

- Base no-LINCS shape: `{input_qc.get("base_shape")}`
- Output shape: `{input_qc.get("output_shape")}`
- Row-level YAPC LINCS coverage: `{_fmt(input_qc.get("row_level_yapc_lincs_coverage"))}`
- YAPC LINCS drugs: `{input_qc.get("yapc_lincs_drugs")}`

## Variant Comparison

{_md_table(comparison, ["variant", "cv", "ensemble", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20"], 20)}

## Random sample 4-fold ensemble

{_md_table(random_ens, ["ensemble", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20"], 10)}

## Drug GroupCV 4-fold ensemble

{_md_table(group_ens, ["ensemble", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20"], 10)}

## Random sample member models

{_md_table(random_ind, ["member", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20"], 20)}

## Drug GroupCV member models

{_md_table(group_ind, ["member", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20"], 20)}

## GroupCV ensemble Top 15

{_md_table(top, ["rank", "drug_name", "ensemble_score", "mean_pred_ln_ic50", "screened_rows", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"], 15)}

## Key outputs

- `data/source_staging/lincs_gse70138/lincs_drug_signature_yapc_phase2_20260422.parquet`
- `data/paad/X_{input_name}.npy`
- `data/paad/{input_name}_feature_names.json`
- `reports/paad/qc_paad_yapc_lincs_signature_20260422.json`
- `reports/paad/qc_paad_yapc_lincs_input_20260422.json`
- `results/paad/ensemble/paad_yapc_lincs_variant_comparison.csv`
- `results/paad/ensemble/groupcv4_drug/{input_name}/ensemble_metrics.csv`
"""
    out_path = paths.root / "docs" / "PAAD_YAPC_LINCS_VARIANT_20260422.md"
    out_path.write_text(text, encoding="utf-8")
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Build and train a GSE70138 YAPC LINCS PAAD feature variant")
    parser.add_argument("--config", default="config/paad_pipeline_config.json")
    parser.add_argument("--lincs-dir", default="data/source_staging/lincs_gse70138")
    parser.add_argument("--cache-dir", default="/tmp/paad_yapc_lincs_gse70138_cache")
    parser.add_argument("--output-input", default=YAPC_INPUT)
    parser.add_argument("--feature-limit", type=int, default=512)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-training", action="store_true")
    args = parser.parse_args()

    cfg = PIPE.load_config(args.config)
    paths = PIPE.pipeline_paths(cfg)
    lincs_paths = lincs_phase2_paths(paths.root, args.lincs_dir)
    missing = [str(path) for key, path in lincs_paths.items() if key != "base" and not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required GSE70138 LINCS files:\n" + "\n".join(missing))

    run_summary: dict[str, Any] = {
        "status": "started",
        "output_input": args.output_input,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    signature, signature_qc = build_yapc_signature(cfg, lincs_paths, Path(args.cache_dir), args.feature_limit)
    run_summary["signature_qc"] = signature_qc
    run_summary["input_qc"] = build_yapc_input(cfg, signature, args.output_input)
    if args.skip_training:
        run_summary["status"] = "completed_input_only"
    else:
        run_summary["cv_results"] = run_training(cfg, args.output_input, force=args.force)
        report = write_variant_report(cfg, args.output_input, run_summary)
        run_summary["report"] = str(report)
        run_summary["status"] = "completed"
    run_summary["completed_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    write_json(paths.reports_dir / "paad_yapc_lincs_variant_run_summary_20260422.json", run_summary)
    print(json.dumps(run_summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
