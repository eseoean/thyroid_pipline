#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

PAAD_SCRIPT_PATH = REPO_ROOT / "scripts" / "11_run_paad_full_pipeline.py"
PAN_INPUT = "numeric_strong_context_smiles_pan_lincs"
PAAD_INPUT = "numeric_strong_context_smiles_paad_lincs"


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


def _split_gene_list(value: object) -> set[str]:
    genes: set[str] = set()
    for part in re.split(r"[;,/| ]+", str(value)):
        part = part.strip().upper()
        if part and part not in {"NAN", "NONE", "<NA>"}:
            genes.add(part)
    return genes


def exact_paad_lincs_qc(paths: Any) -> dict[str, Any]:
    lincs_dir = paths.root / "data" / "source_staging" / "lincs"
    cell_path = lincs_dir / "lincs_cell_info_basic_20260406.parquet"
    sig_path = lincs_dir / "lincs_sig_info_basic_20260406.parquet"
    inst_path = lincs_dir / "lincs_inst_info_basic_20260406.parquet"
    if not cell_path.exists():
        return {"status": "missing_cell_metadata", "cell_path": str(cell_path)}

    cell = pd.read_parquet(cell_path)
    mask = pd.Series(False, index=cell.index)
    for col in cell.columns:
        mask |= cell[col].astype(str).str.contains("pancre|paad|panc", case=False, regex=True, na=False)
    pancreas_cells = cell.loc[mask].copy()
    pancreas_cell_ids = set(pancreas_cells["cell_id"].astype(str))

    sig_rows = 0
    sig_trt_cp_rows = 0
    inst_rows = 0
    inst_trt_cp_rows = 0
    if sig_path.exists() and pancreas_cell_ids:
        sig = pd.read_parquet(sig_path, columns=["sig_id", "pert_id", "pert_iname", "pert_type", "cell_id"])
        sig_sub = sig[sig["cell_id"].astype(str).isin(pancreas_cell_ids)]
        sig_rows = int(len(sig_sub))
        sig_trt_cp_rows = int(sig_sub["pert_type"].astype(str).eq("trt_cp").sum())
    if inst_path.exists() and pancreas_cell_ids:
        inst = pd.read_parquet(inst_path, columns=["inst_id", "pert_id", "pert_iname", "pert_type", "cell_id"])
        inst_sub = inst[inst["cell_id"].astype(str).isin(pancreas_cell_ids)]
        inst_rows = int(len(inst_sub))
        inst_trt_cp_rows = int(inst_sub["pert_type"].astype(str).eq("trt_cp").sum())

    return {
        "status": "checked",
        "definition": "Exact PAAD LINCS means LINCS signatures whose cell metadata primary_site/subtype/cell_id contains pancreas/PAAD/PANC.",
        "pancreas_lincs_cell_count": int(len(pancreas_cells)),
        "pancreas_lincs_cell_ids": sorted(pancreas_cell_ids),
        "pancreas_lincs_cells": pancreas_cells.to_dict(orient="records"),
        "sig_info_rows_for_pancreas_cells": sig_rows,
        "sig_info_trt_cp_rows_for_pancreas_cells": sig_trt_cp_rows,
        "inst_info_rows_for_pancreas_cells": inst_rows,
        "inst_info_trt_cp_rows_for_pancreas_cells": inst_trt_cp_rows,
        "exact_paad_lincs_available": bool(sig_trt_cp_rows > 0),
    }


def collect_paad_evidence_genes(cfg: dict[str, Any], paths: Any) -> dict[str, set[str]]:
    rows = pd.read_csv(paths.processed_dir / "row_metadata.csv")
    expression_path = paths.external_dir / "paad_expression.csv"

    priority = PIPE.priority_genes(cfg)
    target_genes: set[str] = set()
    if "target_genes" in rows.columns:
        for value in rows["target_genes"].dropna().astype(str).unique():
            target_genes.update(_split_gene_list(value))

    tcga_expression: set[str] = set()
    if expression_path.exists():
        expr = pd.read_csv(expression_path, usecols=["Hugo_Symbol"])
        tcga_expression = {str(x).upper() for x in expr["Hugo_Symbol"].dropna()}

    depmap_paad: set[str] = set()
    for col in rows.columns:
        for prefix in (
            "sample__depmap_dependency__",
            "sample__depmap_expression__",
            "sample__depmap_cnv__",
            "sample__depmap_mutation__",
        ):
            if col.startswith(prefix):
                depmap_paad.add(col.replace(prefix, "").upper())

    return {
        "priority": priority,
        "drug_targets": target_genes,
        "tcga_expression": tcga_expression,
        "depmap_paad_features": depmap_paad,
    }


def build_paad_lincs_input(cfg: dict[str, Any], source_input: str, output_input: str) -> dict[str, Any]:
    paths = PIPE.pipeline_paths(cfg)
    X = np.load(paths.processed_dir / f"X_{source_input}.npy").astype(np.float32)
    names = [str(x) for x in read_json(paths.processed_dir / f"{source_input}_feature_names.json")]
    evidence = collect_paad_evidence_genes(cfg, paths)

    lincs_gene_positions: list[tuple[int, str, str]] = []
    lincs_flag_positions: list[int] = []
    for idx, name in enumerate(names):
        if name.startswith("drug__lincs__"):
            gene = name.replace("drug__lincs__", "").upper()
            lincs_gene_positions.append((idx, name, gene))
        elif "lincs" in name.lower():
            lincs_flag_positions.append(idx)

    selected_rows = []
    for idx, name, gene in lincs_gene_positions:
        score = 0
        reasons = []
        if gene in evidence["priority"]:
            score += 3
            reasons.append("paad_priority")
        if gene in evidence["drug_targets"]:
            score += 2
            reasons.append("candidate_target")
        if gene in evidence["tcga_expression"]:
            score += 2
            reasons.append("tcga_paad_expression")
        if gene in evidence["depmap_paad_features"]:
            score += 1
            reasons.append("paad_depmap_feature")
        if score > 0:
            selected_rows.append(
                {
                    "position": idx,
                    "feature_name": name,
                    "gene": gene,
                    "paad_evidence_score": score,
                    "evidence": ";".join(reasons),
                }
            )

    selected = pd.DataFrame(selected_rows).sort_values(
        ["paad_evidence_score", "gene"],
        ascending=[False, True],
    )
    selected_positions = set(selected["position"].astype(int).tolist()) if not selected.empty else set()
    non_lincs_positions = [idx for idx, name in enumerate(names) if "lincs" not in name.lower()]
    keep_positions = non_lincs_positions + [idx for idx, _, _ in lincs_gene_positions if idx in selected_positions] + lincs_flag_positions
    keep_names = [names[idx] for idx in keep_positions]
    X_out = X[:, keep_positions].astype(np.float32)

    np.save(paths.processed_dir / f"X_{output_input}.npy", X_out)
    write_json(paths.processed_dir / f"{output_input}_feature_names.json", keep_names)

    selected_path = paths.reports_dir / "paad_lincs_selected_features_20260422.csv"
    selected_path.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(selected_path, index=False)

    exact_qc = exact_paad_lincs_qc(paths)
    qc = {
        "step": "paad_lincs_variant_input",
        "source_input": source_input,
        "output_input": output_input,
        "definition": "PAAD-focused LINCS uses the pan-cancer drug perturbation values but keeps only LINCS gene features with PAAD evidence: curated PAAD priority genes, candidate drug targets, TCGA-PAAD expression coverage, or PAAD DepMap feature availability.",
        "exact_paad_lincs_qc": exact_qc,
        "source_shape": [int(X.shape[0]), int(X.shape[1])],
        "output_shape": [int(X_out.shape[0]), int(X_out.shape[1])],
        "source_lincs_gene_features": int(len(lincs_gene_positions)),
        "selected_paad_lincs_gene_features": int(len(selected_positions)),
        "kept_lincs_flag_features": int(len(lincs_flag_positions)),
        "non_lincs_features": int(len(non_lincs_positions)),
        "evidence_gene_counts": {key: int(len(value)) for key, value in evidence.items()},
        "selected_features_path": str(selected_path),
        "matrix_path": str(paths.processed_dir / f"X_{output_input}.npy"),
        "feature_names_path": str(paths.processed_dir / f"{output_input}_feature_names.json"),
    }
    write_json(paths.reports_dir / "qc_paad_lincs_variant_20260422.json", qc)
    return qc


def run_training(cfg: dict[str, Any], input_name: str, force: bool) -> dict[str, Any]:
    cv_results: dict[str, Any] = {}
    for cv in ["random4", "groupcv4_drug"]:
        log(f"[paad-lincs] train ML {cv} input={input_name}")
        ml_summary, _ = PIPE.train_ml_cv(cfg, input_name, cv, force=force)
        log(f"[paad-lincs] train DL {cv} input={input_name}")
        dl_summary, _ = PIPE.train_dl_cv(cfg, input_name, cv, force=force)
        log(f"[paad-lincs] build ensemble {cv} input={input_name}")
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
    input_qc = run_summary.get("input_qc", {})
    exact = input_qc.get("exact_paad_lincs_qc", {})

    random_ens = pd.read_csv(paths.results_dir / "ensemble" / "random4" / input_name / "ensemble_metrics.csv")
    group_ens = pd.read_csv(paths.results_dir / "ensemble" / "groupcv4_drug" / input_name / "ensemble_metrics.csv")
    random_ind = pd.read_csv(paths.results_dir / "ensemble" / "random4" / input_name / "individual_metrics.csv")
    group_ind = pd.read_csv(paths.results_dir / "ensemble" / "groupcv4_drug" / input_name / "individual_metrics.csv")
    top = pd.read_csv(paths.results_dir / "ensemble" / "groupcv4_drug" / input_name / "ensemble_top50_drugs.csv")

    pan_group_path = paths.results_dir / "ensemble" / "groupcv4_drug" / PAN_INPUT / "ensemble_metrics.csv"
    pan_random_path = paths.results_dir / "ensemble" / "random4" / PAN_INPUT / "ensemble_metrics.csv"
    comparison_rows = []
    for variant, cv, path in [
        ("pan_cancer_lincs", "random4", pan_random_path),
        ("paad_focused_lincs", "random4", paths.results_dir / "ensemble" / "random4" / input_name / "ensemble_metrics.csv"),
        ("pan_cancer_lincs", "groupcv4_drug", pan_group_path),
        ("paad_focused_lincs", "groupcv4_drug", paths.results_dir / "ensemble" / "groupcv4_drug" / input_name / "ensemble_metrics.csv"),
    ]:
        if path.exists():
            row = pd.read_csv(path).sort_values("spearman", ascending=False).iloc[0].to_dict()
            comparison_rows.append({"variant": variant, "cv": cv, **row})
    comparison = pd.DataFrame(comparison_rows)
    PIPE.write_table(comparison, paths.results_dir / "ensemble" / "paad_lincs_variant_comparison.csv")

    text = f"""# PAAD-Focused LINCS Variant - 2026-04-22

## Scope

- Cancer type: pancreatic cancer / TCGA-PAAD
- Input set: `{input_name}`
- Baseline comparator: `{PAN_INPUT}`
- CV checks: random sample 4-fold and drug GroupCV 4-fold
- Model families: ML, DL, and ML+DL ensembles

## Exact PAAD LINCS Availability

- Pancreas-like LINCS cells in metadata: `{exact.get("pancreas_lincs_cell_ids", [])}`
- `sig_info` rows for those cells: `{exact.get("sig_info_rows_for_pancreas_cells")}`
- `sig_info` trt_cp rows for those cells: `{exact.get("sig_info_trt_cp_rows_for_pancreas_cells")}`
- `inst_info` rows for those cells: `{exact.get("inst_info_rows_for_pancreas_cells")}`
- Exact PAAD-cell LINCS available: `{exact.get("exact_paad_lincs_available")}`

Because exact pancreatic-cell LINCS signatures were absent in the available source metadata, this run uses a PAAD-focused LINCS feature selection:
pan-cancer drug perturbation values are retained only for LINCS genes supported by PAAD priority biology, candidate drug targets, TCGA-PAAD expression coverage, or PAAD DepMap feature availability.

## Input QC

- Source shape: `{input_qc.get("source_shape")}`
- Output shape: `{input_qc.get("output_shape")}`
- Source LINCS gene features: `{input_qc.get("source_lincs_gene_features")}`
- Selected PAAD-focused LINCS gene features: `{input_qc.get("selected_paad_lincs_gene_features")}`
- Kept LINCS flag features: `{input_qc.get("kept_lincs_flag_features")}`
- Selected feature table: `{input_qc.get("selected_features_path")}`

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

- `data/paad/X_{input_name}.npy`
- `data/paad/{input_name}_feature_names.json`
- `reports/paad/qc_paad_lincs_variant_20260422.json`
- `reports/paad/paad_lincs_selected_features_20260422.csv`
- `results/paad/ml/random4/{input_name}/`
- `results/paad/ml/groupcv4_drug/{input_name}/`
- `results/paad/dl/random4/{input_name}/`
- `results/paad/dl/groupcv4_drug/{input_name}/`
- `results/paad/ensemble/random4/{input_name}/`
- `results/paad/ensemble/groupcv4_drug/{input_name}/`
"""
    out_path = paths.root / "docs" / "PAAD_LINCS_VARIANT_20260422.md"
    out_path.write_text(text, encoding="utf-8")
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Build and train a PAAD-focused LINCS feature variant")
    parser.add_argument("--config", default="config/paad_pipeline_config.json")
    parser.add_argument("--source-input", default=PAN_INPUT)
    parser.add_argument("--output-input", default=PAAD_INPUT)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-training", action="store_true")
    args = parser.parse_args()

    cfg = PIPE.load_config(args.config)
    run_summary: dict[str, Any] = {
        "status": "started",
        "source_input": args.source_input,
        "output_input": args.output_input,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    log("[paad-lincs] build PAAD-focused LINCS input")
    run_summary["input_qc"] = build_paad_lincs_input(cfg, args.source_input, args.output_input)
    if args.skip_training:
        run_summary["status"] = "completed_input_only"
    else:
        run_summary["cv_results"] = run_training(cfg, args.output_input, force=args.force)
        report = write_variant_report(cfg, args.output_input, run_summary)
        run_summary["report"] = str(report)
        run_summary["status"] = "completed"
    run_summary["completed_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    paths = PIPE.pipeline_paths(cfg)
    write_json(paths.reports_dir / "paad_lincs_variant_run_summary_20260422.json", run_summary)
    print(json.dumps(run_summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
