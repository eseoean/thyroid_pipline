from __future__ import annotations

import argparse
import html
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import TruncatedSVD
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import mean_absolute_error, mean_squared_error, ndcg_score, r2_score
from sklearn.model_selection import KFold, GroupKFold
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

try:
    from lightgbm import LGBMRegressor
except Exception:  # pragma: no cover - optional dependency
    LGBMRegressor = None

try:
    from xgboost import XGBRegressor
except Exception:  # pragma: no cover - optional dependency
    XGBRegressor = None

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
except Exception:  # pragma: no cover - optional dependency
    Chem = None
    AllChem = None
    DataStructs = None


ID_COLUMNS = {
    "sample_id",
    "cell_line_name",
    "canonical_drug_id",
    "drug_id",
    "drug_name",
    "canonical_smiles",
    "smiles",
    "disease_label",
    "disease_label_normalized",
    "thyroid_subtype",
    "target_genes",
    "PATHWAY_NAME_NORMALIZED",
    "TCGA_DESC",
    "classification",
    "drug_bridge_strength",
    "stage3_resolution_status",
}
LABEL_COLUMNS = {"IC50", "LN_IC50", "ln_ic50", "auc", "AUC", "response", "response_label", "y"}
TEXT_SEPARATORS = re.compile(r"[;,|/]+")


@dataclass
class PipelinePaths:
    root: Path
    raw_dir: Path
    processed_dir: Path
    results_dir: Path
    reports_dir: Path
    external_validation_dir: Path
    admet_output_dir: Path
    knowledge_validation_dir: Path
    phase5_dir: Path
    external_dir: Path
    admet_source_dir: Path


def load_config(config_path: str | Path) -> dict[str, Any]:
    path = Path(config_path).resolve()
    with path.open("r", encoding="utf-8") as f:
        cfg = json.load(f)
    cfg["_config_path"] = str(path)
    cfg["_root"] = str(path.parent.parent)
    return cfg


def _root(cfg: dict[str, Any]) -> Path:
    return Path(cfg.get("_root", ".")).resolve()


def _path(cfg: dict[str, Any], value: str | Path | None) -> Path:
    if value is None:
        return _root(cfg)
    p = Path(value)
    return p if p.is_absolute() else _root(cfg) / p


def pipeline_paths(cfg: dict[str, Any]) -> PipelinePaths:
    paths = cfg["paths"]
    out = PipelinePaths(
        root=_root(cfg),
        raw_dir=_path(cfg, paths["raw_dir"]),
        processed_dir=_path(cfg, paths["processed_dir"]),
        results_dir=_path(cfg, paths["results_dir"]),
        reports_dir=_path(cfg, paths["reports_dir"]),
        external_validation_dir=_path(cfg, paths["external_validation_dir"]),
        admet_output_dir=_path(cfg, paths["admet_output_dir"]),
        knowledge_validation_dir=_path(cfg, paths["knowledge_validation_dir"]),
        phase5_dir=_path(cfg, paths["phase5_dir"]),
        external_dir=_path(cfg, paths["external_dir"]),
        admet_source_dir=_path(cfg, paths["admet_source_dir"]),
    )
    for directory in out.__dict__.values():
        if isinstance(directory, Path):
            directory.mkdir(parents=True, exist_ok=True)
    return out


def write_json(path: Path, payload: dict[str, Any] | list[Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(_json_safe(payload), f, ensure_ascii=False, indent=2)


def read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_safe(v) for v in value]
    if isinstance(value, tuple):
        return [_json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        if not np.isfinite(value):
            return None
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def read_table(path: Path | str | None) -> pd.DataFrame | None:
    if path is None:
        return None
    p = Path(path)
    if not p.exists():
        return None
    suffix = p.suffix.lower()
    if suffix in {".parquet", ".pq"}:
        return pd.read_parquet(p)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(p, sep="\t")
    return pd.read_csv(p)


def write_table(df: pd.DataFrame, path: Path, also_parquet: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    if also_parquet:
        try:
            df.to_parquet(path.with_suffix(".parquet"), index=False)
        except Exception:
            pass


def summarize_series(s: pd.Series) -> dict[str, Any]:
    numeric = pd.to_numeric(s, errors="coerce")
    return {
        "count": int(numeric.notna().sum()),
        "missing": int(numeric.isna().sum()),
        "min": float(numeric.min()) if numeric.notna().any() else None,
        "q1": float(numeric.quantile(0.25)) if numeric.notna().any() else None,
        "median": float(numeric.quantile(0.50)) if numeric.notna().any() else None,
        "q3": float(numeric.quantile(0.75)) if numeric.notna().any() else None,
        "max": float(numeric.max()) if numeric.notna().any() else None,
    }


def safe_spearman(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 3 or np.nanstd(y_true) == 0 or np.nanstd(y_pred) == 0:
        return float("nan")
    return float(stats.spearmanr(y_true, y_pred, nan_policy="omit").correlation)


def safe_pearson(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 3 or np.nanstd(y_true) == 0 or np.nanstd(y_pred) == 0:
        return float("nan")
    return float(stats.pearsonr(y_true, y_pred)[0])


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true = np.asarray(y_true)[mask]
    y_pred = np.asarray(y_pred)[mask]
    if len(y_true) == 0:
        return {"n": 0, "spearman": np.nan, "pearson": np.nan, "rmse": np.nan, "mae": np.nan, "r2": np.nan, "ndcg_at_20": np.nan}
    relevance_values = -y_true
    relevance_values = relevance_values - np.nanmin(relevance_values)
    if np.nanmax(relevance_values) == 0:
        relevance_values = np.ones_like(relevance_values)
    relevance = relevance_values.reshape(1, -1)
    scores = -y_pred.reshape(1, -1)
    try:
        ndcg = float(ndcg_score(relevance, scores, k=min(20, len(y_true))))
    except Exception:
        ndcg = float("nan")
    return {
        "n": int(len(y_true)),
        "spearman": safe_spearman(y_true, y_pred),
        "pearson": safe_pearson(y_true, y_pred),
        "rmse": float(mean_squared_error(y_true, y_pred, squared=False)),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "r2": float(r2_score(y_true, y_pred)) if len(y_true) >= 2 else float("nan"),
        "ndcg_at_20": ndcg,
    }


def canonicalize_smiles(smiles: Any) -> tuple[str, bool]:
    if pd.isna(smiles) or str(smiles).strip() == "":
        return "", False
    text = str(smiles).strip()
    if Chem is None:
        return text, True
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return "", False
    return Chem.MolToSmiles(mol, canonical=True), True


def split_genes(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    parts = [p.strip().upper() for p in TEXT_SEPARATORS.split(str(value)) if p.strip()]
    return sorted(set(parts))


def normalize_disease_label(value: Any, aliases: Iterable[str]) -> tuple[str, bool]:
    text = "" if pd.isna(value) else str(value).strip()
    low = text.lower()
    alias_lows = [a.lower() for a in aliases]
    hit = any(alias in low for alias in alias_lows)
    if "thca" in low:
        hit = True
    return ("THCA" if hit else text.upper(), hit)


def infer_subtype(label: Any, subtype_aliases: dict[str, list[str]]) -> str:
    low = "" if pd.isna(label) else str(label).lower()
    for subtype, aliases in subtype_aliases.items():
        if any(alias.lower() in low for alias in aliases):
            return subtype
    return "Unknown"


def seed_demo_data(config_path: str | Path) -> dict[str, str]:
    cfg = load_config(config_path)
    paths = pipeline_paths(cfg)
    rng = np.random.default_rng(int(cfg["project"].get("seed", 42)))

    samples = pd.DataFrame(
        {
            "sample_id": [f"THCA_CL_{i:02d}" for i in range(1, 13)],
            "cell_line_name": [f"thyroid_cell_{i:02d}" for i in range(1, 13)],
            "disease_label": [
                "papillary thyroid carcinoma",
                "papillary thyroid carcinoma",
                "follicular thyroid carcinoma",
                "thyroid cancer",
                "anaplastic thyroid carcinoma",
                "medullary thyroid carcinoma",
                "thyroid carcinoma",
                "THCA",
                "papillary thyroid carcinoma",
                "thyroid cancer",
                "follicular thyroid carcinoma",
                "anaplastic thyroid carcinoma",
            ],
            "sample__crispr__BRAF": rng.normal(0, 1, 12),
            "sample__crispr__RET": rng.normal(0, 1, 12),
            "sample__crispr__NTRK1": rng.normal(0, 1, 12),
            "sample__expr__MAPK_SCORE": rng.normal(0.5, 0.2, 12),
            "sample__expr__VEGFR_SCORE": rng.normal(0.2, 0.2, 12),
        }
    )

    drugs = pd.DataFrame(
        [
            ("D001", "Lenvatinib", "CC1=C(C(=O)NC2=CC=C(C=C2)OC3=CC=CC=C3)N=C(N1)N", "VEGFR;FGFR;RET", "VEGFR/angiogenesis", "approved"),
            ("D002", "Sorafenib", "CNC(=O)C1=NC=CC(=C1)OC2=CC=C(C=C2)NC(=O)NC3=CC(=C(C=C3)Cl)C(F)(F)F", "BRAF;VEGFR", "MAPK signaling", "approved"),
            ("D003", "Cabozantinib", "COC1=CC2=C(C=C1)N=CN=C2NC3=CC(=C(C=C3)F)OC4=CC=CC=C4", "MET;VEGFR;RET", "VEGFR/angiogenesis", "approved"),
            ("D004", "Vandetanib", "COC1=CC2=C(C=C1OCCCN3CCOCC3)N=CN=C2NC4=CC(=C(C=C4)Br)F", "EGFR;RET;VEGFR", "RET signaling", "approved"),
            ("D005", "Selpercatinib", "CC1=CC(=NC=C1)NC2=NC=NC3=C2C=C(C=C3)OC4=CC=CC=C4", "RET", "RET signaling", "approved"),
            ("D006", "Pralsetinib", "CC1=CC(=NC=C1)NC2=NC=NC3=C2C=C(C=C3)N4CCN(CC4)C", "RET", "RET signaling", "approved"),
            ("D007", "Dabrafenib", "CC(C)(C)NC(=O)C1=CC(=C(C=C1)F)NS(=O)(=O)C2=C(C=CC=C2F)F", "BRAF", "MAPK signaling", "approved"),
            ("D008", "Trametinib", "CC1=C(C(=O)N(N1C2=CC=CC=C2)C3=CC=C(C=C3)I)C(=O)NC4=CC=CC=C4", "MEK1;MEK2", "MAPK signaling", "approved"),
            ("D009", "Larotrectinib", "CC1=C(C=CC(=N1)NC2=NC=NC3=C2C=CN3)C4=CN=CC=C4", "NTRK1;NTRK2;NTRK3", "NTRK fusion", "approved"),
            ("D010", "Entrectinib", "CN1CCN(CC1)C2=NC3=C(C=NN3C4=CC=CC=C4)C(=N2)C5=CC=CC=C5", "NTRK1;ROS1;ALK", "NTRK fusion", "approved"),
            ("D011", "Everolimus", "CC1CC(C(C(C1)OC2CC(C(C(O2)C)O)OC)C)OC", "mTOR", "PI3K/AKT/mTOR", "indication_expansion"),
            ("D012", "Palbociclib", "CC1=C(C(=O)N=C(N1)N)C2=CC=C(C=C2)N3CCN(CC3)C", "CDK4;CDK6", "cell-cycle", "indication_expansion"),
            ("D013", "Olaparib", "C1CC1C(=O)N2CCC(CC2)NC(=O)C3=CC=CC=C3", "PARP1;PARP2", "DNA damage", "indication_expansion"),
            ("D014", "Alpelisib", "CC1=NC(=NC=C1)N2C=NC3=C2C=CC(=C3)S(=O)(=O)N", "PIK3CA", "PI3K/AKT/mTOR", "indication_expansion"),
            ("D015", "Gemcitabine", "C1=NC(=O)N(C=C1F)C2C(C(C(O2)CO)O)F", "RRM1;RRM2", "DNA damage", "exploratory"),
            ("D016", "Docetaxel", "CC(C)C1=CC=C(C=C1)C(C)C(=O)O", "TUBB", "cell-cycle", "exploratory"),
            ("D017", "Dasatinib", "CC1=NC(=CC=N1)NC2=NC=C(S2)C3=CC=CC=C3", "SRC;ABL1", "kinase signaling", "exploratory"),
            ("D018", "Temsirolimus", "CC1CC(C(C(C1)OC2CC(C(C(O2)C)O)OC)C)O", "mTOR", "PI3K/AKT/mTOR", "indication_expansion"),
        ],
        columns=["canonical_drug_id", "drug_name", "canonical_smiles", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"],
    )

    for i in range(32):
        drugs[f"drug_morgan_{i:03d}"] = rng.integers(0, 2, len(drugs))
    for col, loc in [
        ("drug_desc_hba", 4),
        ("drug_desc_hbd", 2),
        ("drug_desc_heavy_atoms", 28),
        ("drug_desc_ring_count", 3),
        ("drug_desc_rot_bonds", 5),
    ]:
        drugs[col] = np.maximum(0, rng.normal(loc, 1, len(drugs))).round(2)
    drugs["drug__lincs_score"] = rng.normal(0, 1, len(drugs))
    drugs["drug_bridge_strength"] = np.where(drugs["classification"].eq("approved"), "multi_source", "single_source")
    drugs["stage3_resolution_status"] = np.where(drugs["target_genes"].str.len() > 0, "resolved", "target_partial")
    drugs["TCGA_DESC"] = "THCA"

    pairs = samples[["sample_id", "cell_line_name", "disease_label"]].merge(
        drugs[["canonical_drug_id", "drug_name"]], how="cross"
    )
    drug_effect = dict(zip(drugs["canonical_drug_id"], rng.normal(0, 0.7, len(drugs))))
    for did in ["D001", "D002", "D003", "D004", "D005", "D007", "D008", "D009", "D010"]:
        drug_effect[did] -= 1.2
    sample_effect = dict(zip(samples["sample_id"], rng.normal(0, 0.3, len(samples))))
    pairs["LN_IC50"] = [
        3.5 + drug_effect[row.canonical_drug_id] + sample_effect[row.sample_id] + rng.normal(0, 0.18)
        for row in pairs.itertuples()
    ]

    expression_genes = ["BRAF", "RET", "NTRK1", "NTRK2", "NTRK3", "VEGFR", "FGFR", "MTOR", "PIK3CA", "CDK4", "CDK6", "PARP1"]
    expression = pd.DataFrame({"Hugo_Symbol": expression_genes})
    for i in range(1, 41):
        expression[f"THCA_PAT_{i:03d}"] = rng.gamma(2.0, 1.0, len(expression))
    clinical = pd.DataFrame(
        {
            "PATIENT_ID": [f"THCA_PAT_{i:03d}" for i in range(1, 41)],
            "OS_MONTHS": rng.gamma(4.0, 12.0, 40).round(2),
            "OS_STATUS": rng.integers(0, 2, 40),
            "RFS_MONTHS": rng.gamma(3.0, 10.0, 40).round(2),
            "RFS_STATUS": rng.integers(0, 2, 40),
        }
    )

    admet_smiles = drugs[["canonical_smiles"]].rename(columns={"canonical_smiles": "smiles"}).copy()
    admet_smiles["label"] = rng.integers(0, 2, len(admet_smiles))
    for assay in ["ames", "dili", "herg", "solubility", "cyp3a4"]:
        table = admet_smiles.copy()
        if assay in {"ames", "dili", "herg"}:
            table["label"] = rng.binomial(1, 0.25, len(table))
        else:
            table["label"] = rng.normal(0, 1, len(table)).round(3)
        write_table(table, paths.admet_source_dir / f"{assay}.csv")

    write_table(samples, paths.raw_dir / "sample_features.csv")
    write_table(drugs, paths.raw_dir / "drug_features.csv")
    write_table(drugs[["canonical_drug_id", "drug_name", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification", "drug_bridge_strength", "stage3_resolution_status", "TCGA_DESC"]], paths.raw_dir / "drug_annotations.csv")
    write_table(pairs, paths.raw_dir / "thyroid_response_pairs.csv")
    write_table(expression, paths.external_dir / "thyroid_expression.csv")
    write_table(clinical, paths.external_dir / "thyroid_clinical.csv")

    return {
        "response_table": str(paths.raw_dir / "thyroid_response_pairs.csv"),
        "sample_features": str(paths.raw_dir / "sample_features.csv"),
        "drug_features": str(paths.raw_dir / "drug_features.csv"),
        "drug_annotations": str(paths.raw_dir / "drug_annotations.csv"),
        "external_expression": str(paths.external_dir / "thyroid_expression.csv"),
        "external_clinical": str(paths.external_dir / "thyroid_clinical.csv"),
        "admet_dir": str(paths.admet_source_dir),
    }


def prepare_response_subset(cfg: dict[str, Any]) -> pd.DataFrame:
    paths = pipeline_paths(cfg)
    inputs = cfg["inputs"]
    response = read_table(_path(cfg, inputs["response_table"]))
    if response is None:
        raise FileNotFoundError(f"Missing response table: {inputs['response_table']}")

    sample_features = read_table(_path(cfg, inputs.get("sample_features")))
    drug_features = read_table(_path(cfg, inputs.get("drug_features")))
    drug_annotations = read_table(_path(cfg, inputs.get("drug_annotations")))

    response = response.copy()
    if "sample_id" not in response.columns or "canonical_drug_id" not in response.columns:
        raise ValueError("response_table must contain sample_id and canonical_drug_id")

    label_col = "LN_IC50" if "LN_IC50" in response.columns else None
    if label_col is None and "IC50" in response.columns:
        label_col = "IC50"
        response["LN_IC50"] = np.log1p(pd.to_numeric(response["IC50"], errors="coerce"))
        label_col = "LN_IC50"
    if label_col is None:
        raise ValueError("response_table must contain LN_IC50 or IC50")

    if "disease_label" not in response.columns and sample_features is not None and "disease_label" in sample_features.columns:
        response = response.merge(sample_features[["sample_id", "disease_label"]].drop_duplicates(), on="sample_id", how="left")
    if "disease_label" not in response.columns:
        response["disease_label"] = cfg["project"].get("disease", "thyroid cancer")

    normalized = response["disease_label"].map(lambda x: normalize_disease_label(x, cfg["disease_aliases"]))
    response["disease_label_normalized"] = [x[0] for x in normalized]
    response["_thyroid_alias_hit"] = [x[1] for x in normalized]
    response["thyroid_subtype"] = response["disease_label"].map(lambda x: infer_subtype(x, cfg["thyroid_subtype_aliases"]))

    original_rows = len(response)
    subset = response.loc[response["_thyroid_alias_hit"]].copy()
    subset[label_col] = pd.to_numeric(subset[label_col], errors="coerce")
    label_missing = int(subset[label_col].isna().sum())
    subset = subset.loc[subset[label_col].notna()].copy()

    duplicate_rows = int(subset.duplicated(["sample_id", "canonical_drug_id"]).sum())
    grouped = (
        subset.groupby(["sample_id", "canonical_drug_id"], as_index=False)
        .agg(
            {
                label_col: "mean",
                "cell_line_name": "first" if "cell_line_name" in subset.columns else "size",
                "drug_name": "first" if "drug_name" in subset.columns else "size",
                "disease_label": "first",
                "disease_label_normalized": "first",
                "thyroid_subtype": "first",
            }
        )
        .rename(columns={label_col: "LN_IC50"})
    )
    if "cell_line_name" not in grouped.columns or pd.api.types.is_numeric_dtype(grouped["cell_line_name"]):
        grouped["cell_line_name"] = grouped["sample_id"]
    if "drug_name" not in grouped.columns or pd.api.types.is_numeric_dtype(grouped["drug_name"]):
        grouped["drug_name"] = grouped["canonical_drug_id"]

    merged = grouped
    if sample_features is not None:
        sample_cols = [c for c in sample_features.columns if c not in {"cell_line_name", "disease_label"}]
        merged = merged.merge(sample_features[sample_cols].drop_duplicates("sample_id"), on="sample_id", how="left")
    if drug_features is not None:
        drug_cols = [c for c in drug_features.columns if c not in {"drug_name"}]
        merged = merged.merge(drug_features[drug_cols].drop_duplicates("canonical_drug_id"), on="canonical_drug_id", how="left")
    if drug_annotations is not None:
        ann = drug_annotations.drop_duplicates("canonical_drug_id")
        duplicate_ann_cols = [c for c in ann.columns if c in merged.columns and c != "canonical_drug_id"]
        ann = ann.drop(columns=duplicate_ann_cols)
        merged = merged.merge(ann, on="canonical_drug_id", how="left")

    for context_col in cfg["features"]["context_columns"]:
        if context_col not in merged.columns:
            merged[context_col] = "THCA" if context_col == "TCGA_DESC" else "UNCLASSIFIED"
    if "target_genes" not in merged.columns:
        merged["target_genes"] = ""
    if "canonical_smiles" not in merged.columns and "smiles" in merged.columns:
        merged["canonical_smiles"] = merged["smiles"]
    elif "canonical_smiles" not in merged.columns:
        merged["canonical_smiles"] = ""

    qc = {
        "step": "step1_response_subset",
        "original_rows": int(original_rows),
        "thyroid_alias_rows": int(len(subset) + label_missing),
        "rows_after_label_filter": int(len(subset)),
        "rows_after_dedup": int(len(merged)),
        "unique_cell_lines": int(merged["sample_id"].nunique()),
        "unique_drugs": int(merged["canonical_drug_id"].nunique()),
        "duplicate_cell_line_drug_rows_before_dedup": duplicate_rows,
        "label_missing": label_missing,
        "label_summary": summarize_series(merged["LN_IC50"]),
        "subtype_counts": merged["thyroid_subtype"].value_counts(dropna=False).to_dict(),
        "drug_count_summary": summarize_series(merged.groupby("canonical_drug_id").size()),
        "cell_line_count_summary": summarize_series(merged.groupby("sample_id").size()),
    }
    write_table(merged, paths.processed_dir / "thyroid_response_pairs.csv", also_parquet=True)
    write_table(merged, paths.processed_dir / "row_metadata.csv", also_parquet=True)
    write_json(paths.reports_dir / "qc_step1_response_subset.json", qc)
    write_simple_html_qc(paths.reports_dir / "qc_step1_response_subset.html", "Step 1 Response Subset QC", qc)
    return merged


def build_numeric_base(cfg: dict[str, Any], rows: pd.DataFrame | None = None) -> tuple[np.ndarray, np.ndarray, list[str]]:
    paths = pipeline_paths(cfg)
    if rows is None:
        rows = read_table(paths.processed_dir / "row_metadata.csv")
    if rows is None:
        raise FileNotFoundError("row_metadata.csv is missing; run prepare_response_subset first")

    physchem = set(cfg["features"]["physchem_context_columns"])
    excluded = ID_COLUMNS | LABEL_COLUMNS | physchem | {"_thyroid_alias_hit"}
    numeric_cols = []
    for col in rows.columns:
        if col in excluded:
            continue
        if pd.api.types.is_numeric_dtype(rows[col]):
            numeric_cols.append(col)

    numeric = rows[numeric_cols].apply(pd.to_numeric, errors="coerce").astype(float)
    nan_before = int(numeric.isna().sum().sum())
    inf_before = int(np.isinf(numeric.to_numpy()).sum())
    numeric = numeric.replace([np.inf, -np.inf], np.nan)
    medians = numeric.median(axis=0, numeric_only=True).fillna(0)
    numeric = numeric.fillna(medians).fillna(0)
    values = numeric.to_numpy(dtype=np.float32)

    constant_mask = np.nanstd(values, axis=0) == 0
    all_zero_mask = np.all(values == 0, axis=0)
    removed_constant = []
    if cfg["features"].get("remove_constant_numeric_features", True) and values.shape[1] > 0:
        keep = ~constant_mask
        removed_constant = [c for c, drop in zip(numeric_cols, ~keep) if drop]
        values = values[:, keep]
        numeric_cols = [c for c, k in zip(numeric_cols, keep) if k]

    y = rows["LN_IC50"].to_numpy(dtype=np.float32)
    np.save(paths.processed_dir / "X_numeric.npy", values)
    np.save(paths.processed_dir / "y_train.npy", y)
    write_json(paths.processed_dir / "numeric_feature_names.json", numeric_cols)

    qc = {
        "step": "step2_numeric_base",
        "policy": cfg["features"]["numeric_base_policy"],
        "row_count": int(values.shape[0]),
        "feature_count": int(values.shape[1]),
        "excluded_physchem_columns": sorted([c for c in physchem if c in rows.columns]),
        "nan_before_imputation": nan_before,
        "inf_before_imputation": inf_before,
        "all_zero_columns_before_constant_filter": int(all_zero_mask.sum()) if len(all_zero_mask) else 0,
        "constant_columns_removed": len(removed_constant),
        "removed_constant_feature_names": removed_constant[:200],
        "expected_direction": "Choi-style numeric base; BRCA reference was 5524 after excluding five physchem columns",
    }
    write_json(paths.reports_dir / "qc_step2_numeric_base.json", qc)
    return values, y, numeric_cols


def build_smiles_svd(cfg: dict[str, Any], rows: pd.DataFrame | None = None) -> tuple[np.ndarray, list[str]]:
    paths = pipeline_paths(cfg)
    if rows is None:
        rows = read_table(paths.processed_dir / "row_metadata.csv")
    if rows is None:
        raise FileNotFoundError("row_metadata.csv is missing; run prepare_response_subset first")

    dim = int(cfg["features"].get("smiles_svd_dim", 64))
    drug_smiles = rows[["canonical_drug_id", "drug_name", "canonical_smiles"]].drop_duplicates("canonical_drug_id").copy()
    canonical = drug_smiles["canonical_smiles"].map(canonicalize_smiles)
    drug_smiles["smiles_canonicalized"] = [x[0] for x in canonical]
    drug_smiles["rdkit_parse_ok"] = [x[1] for x in canonical]

    valid = drug_smiles.loc[drug_smiles["smiles_canonicalized"].str.len() > 0].copy()
    if len(valid) >= 2:
        vectorizer = TfidfVectorizer(analyzer="char", ngram_range=(2, 4), lowercase=False)
        tfidf = vectorizer.fit_transform(valid["smiles_canonicalized"])
        n_components = min(dim, max(1, min(tfidf.shape[0] - 1, tfidf.shape[1] - 1)))
        if n_components >= 1:
            svd = TruncatedSVD(n_components=n_components, random_state=int(cfg["project"].get("seed", 42)))
            emb_valid = svd.fit_transform(tfidf).astype(np.float32)
            explained = float(np.sum(svd.explained_variance_ratio_))
        else:
            emb_valid = np.zeros((len(valid), 0), dtype=np.float32)
            explained = 0.0
    else:
        emb_valid = np.zeros((len(valid), 0), dtype=np.float32)
        explained = 0.0

    if emb_valid.shape[1] < dim:
        emb_valid = np.pad(emb_valid, ((0, 0), (0, dim - emb_valid.shape[1])), constant_values=0)
    elif emb_valid.shape[1] > dim:
        emb_valid = emb_valid[:, :dim]

    emb = pd.DataFrame(0.0, index=drug_smiles["canonical_drug_id"], columns=[f"smiles_svd_{i:02d}" for i in range(dim)])
    if len(valid) > 0:
        emb.loc[valid["canonical_drug_id"], :] = emb_valid
    mapped = rows[["canonical_drug_id"]].join(emb, on="canonical_drug_id")
    X = mapped[[f"smiles_svd_{i:02d}" for i in range(dim)]].to_numpy(dtype=np.float32)
    feature_names = [f"smiles_svd_{i:02d}" for i in range(dim)]

    np.save(paths.processed_dir / "X_smiles_svd64.npy", X)
    write_json(paths.processed_dir / "smiles_feature_names.json", feature_names)
    write_table(drug_smiles, paths.processed_dir / "drug_smiles_qc.csv")
    qc = {
        "step": "step3_smiles_svd",
        "drug_count": int(len(drug_smiles)),
        "smiles_present_drugs": int((drug_smiles["canonical_smiles"].fillna("").astype(str).str.len() > 0).sum()),
        "rdkit_parse_success_drugs": int(drug_smiles["rdkit_parse_ok"].sum()),
        "row_count": int(X.shape[0]),
        "feature_count": int(X.shape[1]),
        "svd_explained_variance_sum": explained,
        "all_zero_row_count": int(np.all(X == 0, axis=1).sum()),
        "missing_smiles_drugs": drug_smiles.loc[~drug_smiles["rdkit_parse_ok"], "drug_name"].tolist(),
        "duplicate_smiles_count": int(drug_smiles["smiles_canonicalized"].duplicated().sum()),
    }
    write_json(paths.reports_dir / "qc_step3_smiles.json", qc)
    return X, feature_names


def build_strong_context(cfg: dict[str, Any], rows: pd.DataFrame | None = None) -> tuple[np.ndarray, list[str]]:
    paths = pipeline_paths(cfg)
    if rows is None:
        rows = read_table(paths.processed_dir / "row_metadata.csv")
    if rows is None:
        raise FileNotFoundError("row_metadata.csv is missing; run prepare_response_subset first")

    context_cols = cfg["features"]["context_columns"]
    target_dim = int(cfg["features"].get("strong_context_dim", 32))
    context = pd.DataFrame(index=rows.index)
    for col in context_cols:
        if col == "TCGA_DESC":
            default = cfg["project"].get("tcga_code", "THCA")
        else:
            default = "UNCLASSIFIED"
        values = rows[col] if col in rows.columns else default
        context[col] = pd.Series(values).fillna(default).astype(str).str.strip().replace("", default)

    try:
        encoder = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    except TypeError:
        encoder = OneHotEncoder(sparse=False, handle_unknown="ignore")
    onehot = encoder.fit_transform(context).astype(np.float32)
    raw_names = []
    for col, cats in zip(context_cols, encoder.categories_):
        raw_names.extend([f"{col}={cat}" for cat in cats])

    hashed_due_to_dim = False
    if onehot.shape[1] > target_dim:
        hashed_due_to_dim = True
        compressed = np.zeros((onehot.shape[0], target_dim), dtype=np.float32)
        for j in range(onehot.shape[1]):
            compressed[:, j % target_dim] += onehot[:, j]
        X = compressed
        feature_names = [f"strong_context_hash_{i:02d}" for i in range(target_dim)]
    elif onehot.shape[1] < target_dim:
        X = np.pad(onehot, ((0, 0), (0, target_dim - onehot.shape[1])), constant_values=0)
        feature_names = raw_names + [f"strong_context_pad_{i:02d}" for i in range(target_dim - onehot.shape[1])]
    else:
        X = onehot
        feature_names = raw_names

    mapping = context.copy()
    mapping["sample_id"] = rows["sample_id"].values
    mapping["canonical_drug_id"] = rows["canonical_drug_id"].values
    mapping["drug_name"] = rows["drug_name"].values
    write_table(mapping, paths.processed_dir / "context_mapping_table.csv")
    np.save(paths.processed_dir / "X_strong_context.npy", X.astype(np.float32))
    write_json(paths.processed_dir / "strong_context_feature_names.json", feature_names)
    qc = {
        "step": "step4_strong_context",
        "row_count": int(X.shape[0]),
        "feature_count": int(X.shape[1]),
        "raw_onehot_feature_count": int(onehot.shape[1]),
        "target_dim": target_dim,
        "hashed_due_to_dim": hashed_due_to_dim,
        "unique_class_counts": {col: int(context[col].nunique()) for col in context_cols},
        "unknown_or_unclassified_ratio": {
            col: float(context[col].str.upper().isin(["UNKNOWN", "UNCLASSIFIED", ""]).mean())
            for col in context_cols
        },
    }
    write_json(paths.reports_dir / "qc_step4_strong_context.json", qc)
    return X.astype(np.float32), feature_names


def assemble_model_inputs(cfg: dict[str, Any]) -> dict[str, Path]:
    paths = pipeline_paths(cfg)
    X_num = np.load(paths.processed_dir / "X_numeric.npy")
    X_smiles = np.load(paths.processed_dir / "X_smiles_svd64.npy")
    X_ctx = np.load(paths.processed_dir / "X_strong_context.npy")
    y = np.load(paths.processed_dir / "y_train.npy")
    numeric_names = read_json(paths.processed_dir / "numeric_feature_names.json")
    smiles_names = read_json(paths.processed_dir / "smiles_feature_names.json")
    context_names = read_json(paths.processed_dir / "strong_context_feature_names.json")

    if not (X_num.shape[0] == X_smiles.shape[0] == X_ctx.shape[0] == len(y)):
        raise ValueError("Input row counts do not match")

    X_num_smiles = np.concatenate([X_num, X_smiles], axis=1).astype(np.float32)
    X_full = np.concatenate([X_num, X_ctx, X_smiles], axis=1).astype(np.float32)
    np.save(paths.processed_dir / "X_numeric_smiles.npy", X_num_smiles)
    np.save(paths.processed_dir / "X_numeric_strong_context_smiles.npy", X_full)
    write_json(paths.processed_dir / "numeric_smiles_feature_names.json", numeric_names + smiles_names)
    write_json(paths.processed_dir / "numeric_strong_context_smiles_feature_names.json", numeric_names + context_names + smiles_names)
    qc = {
        "step": "step5_model_inputs",
        "row_count": int(len(y)),
        "numeric_dim": int(X_num.shape[1]),
        "smiles_dim": int(X_smiles.shape[1]),
        "strong_context_dim": int(X_ctx.shape[1]),
        "numeric_smiles_dim": int(X_num_smiles.shape[1]),
        "numeric_strong_context_smiles_dim": int(X_full.shape[1]),
        "guidance_expected_full_dim_when_numeric_5524": 5620,
        "dtype": str(X_full.dtype),
    }
    write_json(paths.reports_dir / "qc_step5_model_inputs.json", qc)
    return {
        "numeric": paths.processed_dir / "X_numeric.npy",
        "numeric_smiles": paths.processed_dir / "X_numeric_smiles.npy",
        "numeric_strong_context_smiles": paths.processed_dir / "X_numeric_strong_context_smiles.npy",
    }


def make_model(name: str, seed: int, n_rows: int) -> Any | None:
    if name == "LightGBM":
        if LGBMRegressor is None:
            return None
        return LGBMRegressor(
            n_estimators=180,
            learning_rate=0.04,
            max_depth=-1,
            num_leaves=31,
            subsample=0.85,
            colsample_bytree=0.85,
            random_state=seed,
            verbose=-1,
        )
    if name == "LightGBM_DART":
        if LGBMRegressor is None:
            return None
        return LGBMRegressor(
            boosting_type="dart",
            n_estimators=160,
            learning_rate=0.04,
            num_leaves=31,
            subsample=0.85,
            colsample_bytree=0.85,
            random_state=seed,
            verbose=-1,
        )
    if name == "XGBoost":
        if XGBRegressor is None:
            return None
        return XGBRegressor(
            n_estimators=180,
            max_depth=4,
            learning_rate=0.04,
            subsample=0.85,
            colsample_bytree=0.85,
            objective="reg:squarederror",
            random_state=seed,
            n_jobs=2,
            verbosity=0,
        )
    if name == "ExtraTrees":
        return ExtraTreesRegressor(n_estimators=180, min_samples_leaf=2, random_state=seed, n_jobs=-1)
    if name == "RandomForest":
        return RandomForestRegressor(n_estimators=180, min_samples_leaf=2, random_state=seed, n_jobs=-1)
    if name == "FlatMLP":
        hidden = (96, 48) if n_rows >= 120 else (48, 24)
        return make_pipeline(
            StandardScaler(),
            MLPRegressor(hidden_layer_sizes=hidden, alpha=1e-4, learning_rate_init=1e-3, max_iter=350, random_state=seed, early_stopping=True),
        )
    return None


def train_models(cfg: dict[str, Any], input_paths: dict[str, Path] | None = None) -> dict[str, pd.DataFrame]:
    paths = pipeline_paths(cfg)
    if input_paths is None:
        input_paths = {
            "numeric": paths.processed_dir / "X_numeric.npy",
            "numeric_smiles": paths.processed_dir / "X_numeric_smiles.npy",
            "numeric_strong_context_smiles": paths.processed_dir / "X_numeric_strong_context_smiles.npy",
        }
    y = np.load(paths.processed_dir / "y_train.npy")
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    seed = int(cfg["training"].get("random_state", cfg["project"].get("seed", 42)))
    n_splits = min(int(cfg["training"].get("n_splits", 3)), len(y))
    if n_splits < 2:
        raise ValueError("At least two rows are required for cross-validation")
    models = cfg["training"]["models"]
    summaries: dict[str, pd.DataFrame] = {}
    qc: dict[str, Any] = {"step": "step6_training", "input_sets": {}, "skipped_models": []}

    for input_name, input_path in input_paths.items():
        X = np.load(input_path)
        out_dir = paths.results_dir / "random3" / input_name
        oof_dir = paths.results_dir / "oof" / input_name
        out_dir.mkdir(parents=True, exist_ok=True)
        oof_dir.mkdir(parents=True, exist_ok=True)
        fold_rows = []
        summary_rows = []
        kfold = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
        for model_name in models:
            model_template = make_model(model_name, seed, len(y))
            if model_template is None:
                qc["skipped_models"].append({"input_set": input_name, "model": model_name, "reason": "dependency_missing_or_not_implemented"})
                continue
            oof = np.zeros(len(y), dtype=np.float32)
            train_metrics = []
            val_metrics = []
            for fold, (train_idx, valid_idx) in enumerate(kfold.split(X), start=1):
                model = make_model(model_name, seed + fold, len(y))
                model.fit(X[train_idx], y[train_idx])
                train_pred = np.asarray(model.predict(X[train_idx]), dtype=float)
                valid_pred = np.asarray(model.predict(X[valid_idx]), dtype=float)
                oof[valid_idx] = valid_pred.astype(np.float32)
                tr = regression_metrics(y[train_idx], train_pred)
                va = regression_metrics(y[valid_idx], valid_pred)
                train_metrics.append(tr)
                val_metrics.append(va)
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
                "fold_metrics": {"train": train_metrics, "valid": val_metrics},
                "oof_metrics": oof_metrics,
                "train_oof_spearman_gap": train_gap,
                "prediction_variance": float(np.nanvar(oof)),
            }
            np.save(oof_dir / f"{model_name}.npy", oof)
            write_json(out_dir / f"{input_name}_{model_name}_random3.json", payload)
            row = {"input_set": input_name, "model": model_name, **oof_metrics, "train_oof_spearman_gap": train_gap, "prediction_variance": float(np.nanvar(oof))}
            summary_rows.append(row)
        summary = pd.DataFrame(summary_rows).sort_values("spearman", ascending=False, na_position="last")
        write_table(summary, paths.results_dir / f"{input_name}_metrics_summary.csv")
        write_table(pd.DataFrame(fold_rows), paths.results_dir / f"{input_name}_fold_metrics.csv")
        summaries[input_name] = summary
        qc["input_sets"][input_name] = {
            "shape": [int(X.shape[0]), int(X.shape[1])],
            "models_completed": summary["model"].tolist() if len(summary) else [],
            "fold_count": int(n_splits),
            "label_summary": summarize_series(pd.Series(y)),
            "drug_coverage_per_fold": _fold_coverage(rows, kfold.split(X), "canonical_drug_id") if rows is not None else [],
            "cell_line_coverage_per_fold": _fold_coverage(rows, kfold.split(X), "sample_id") if rows is not None else [],
        }
    write_json(paths.reports_dir / "qc_step6_training.json", qc)
    return summaries


def _fold_coverage(rows: pd.DataFrame, splits: Iterable[tuple[np.ndarray, np.ndarray]], col: str) -> list[dict[str, int]]:
    coverage = []
    for fold, (_, valid_idx) in enumerate(splits, start=1):
        coverage.append({"fold": fold, f"unique_{col}": int(rows.iloc[valid_idx][col].nunique())})
    return coverage


def groupcv_stress_test(cfg: dict[str, Any], input_name: str = "numeric_strong_context_smiles") -> dict[str, Any]:
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
    oof = np.zeros(len(y), dtype=np.float32)
    fold_metrics = []
    for fold, (train_idx, valid_idx) in enumerate(GroupKFold(n_splits=n_splits).split(X, y, groups), start=1):
        model = make_model(model_name, int(cfg["project"].get("seed", 42)) + fold, len(y))
        model.fit(X[train_idx], y[train_idx])
        pred = model.predict(X[valid_idx])
        oof[valid_idx] = pred
        fold_metrics.append({"fold": fold, **regression_metrics(y[valid_idx], pred)})
    out_dir = paths.results_dir / "groupcv_stress_test"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"{input_name}_{model_name}_groupcv_oof.npy", oof)
    result = {"status": "completed", "input_set": input_name, "model": model_name, "fold_metrics": fold_metrics, "oof_metrics": regression_metrics(y, oof)}
    write_json(out_dir / f"{input_name}_{model_name}_groupcv.json", result)
    return result


def ensemble_and_diversity(cfg: dict[str, Any], input_name: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = pipeline_paths(cfg)
    input_name = input_name or cfg["training"].get("primary_input_set", "numeric_strong_context_smiles")
    y = np.load(paths.processed_dir / "y_train.npy")
    rows = read_table(paths.processed_dir / "row_metadata.csv")
    metrics_path = paths.results_dir / f"{input_name}_metrics_summary.csv"
    metrics = read_table(metrics_path)
    if metrics is None or metrics.empty:
        raise FileNotFoundError(f"No metrics found for input set {input_name}; run train_models first")
    oof_dir = paths.results_dir / "oof" / input_name
    preds = {}
    weights = {}
    for row in metrics.itertuples():
        model_name = row.model
        path = oof_dir / f"{model_name}.npy"
        if not path.exists():
            continue
        pred = np.load(path)
        preds[model_name] = pred
        score = getattr(row, "spearman", np.nan)
        weights[model_name] = max(0.0, float(score)) if np.isfinite(score) else 0.0
    if not preds:
        raise ValueError("No OOF predictions found for ensemble")
    weight_sum = sum(weights.values())
    if weight_sum <= 0:
        weights = {k: 1.0 / len(preds) for k in preds}
    else:
        weights = {k: v / weight_sum for k, v in weights.items()}
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
    available_group_cols = [c for c in group_cols if c in candidate_rows.columns]
    top = (
        candidate_rows.groupby(available_group_cols, dropna=False)
        .agg(mean_pred_ln_ic50=("ensemble_pred_ln_ic50", "mean"), ensemble_score=("ensemble_score", "mean"), screened_rows=("ensemble_score", "size"))
        .reset_index()
        .sort_values("ensemble_score", ascending=False)
    )
    top["rank"] = np.arange(1, len(top) + 1)
    top = top[["rank"] + [c for c in top.columns if c != "rank"]]

    out_dir = paths.results_dir / "ensemble"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"{input_name}_weighted_ensemble_oof.npy", ensemble_pred)
    write_table(top.head(30), out_dir / "thyroid_ensemble_top30_drugs.csv")
    write_table(diversity, out_dir / "thyroid_ensemble_diversity.csv")
    result = {
        "input_set": input_name,
        "weights": weights,
        "ensemble_metrics": regression_metrics(y, ensemble_pred),
        "best_single_model": metrics.iloc[0].to_dict(),
        "diversity_summary": diversity.describe(include="all").to_dict() if not diversity.empty else {},
    }
    write_json(out_dir / "thyroid_ensemble_results.json", result)
    write_json(paths.reports_dir / "qc_step7_ensemble_diversity.json", result)
    return top, diversity


def external_validation(cfg: dict[str, Any], top_candidates: pd.DataFrame | None = None) -> pd.DataFrame:
    paths = pipeline_paths(cfg)
    if top_candidates is None:
        top_candidates = read_table(paths.results_dir / "ensemble" / "thyroid_ensemble_top30_drugs.csv")
    if top_candidates is None:
        raise FileNotFoundError("Top candidate file missing; run ensemble first")

    expr = read_table(_path(cfg, cfg["inputs"].get("external_expression")))
    clin = read_table(_path(cfg, cfg["inputs"].get("external_clinical")))
    out_dir = paths.external_validation_dir
    known = {x.upper() for x in cfg["known_thyroid_drugs"]}

    if expr is None or expr.empty:
        validation = top_candidates.head(15).copy()
        validation["target_match_genes"] = ""
        validation["target_expression_pct"] = np.nan
        validation["target_expressed"] = False
        validation["survival_p_value"] = np.nan
        validation["known_thyroid_control"] = validation["drug_name"].astype(str).str.upper().isin(known)
        write_table(validation, out_dir / "top15_validated.csv")
        qc = {"step": "step8_external_validation", "status": "no_external_expression", "topk_precision": _known_precision(top_candidates, known)}
        write_json(paths.reports_dir / "qc_step8_external_validation.json", qc)
        return validation

    gene_col = "Hugo_Symbol" if "Hugo_Symbol" in expr.columns else expr.columns[0]
    expr = expr.copy()
    expr[gene_col] = expr[gene_col].astype(str).str.upper()
    patient_cols = [c for c in expr.columns if c != gene_col]
    expr_values = expr.set_index(gene_col)[patient_cols].apply(pd.to_numeric, errors="coerce")
    global_median = float(np.nanmedian(expr_values.to_numpy()))
    expr_index = set(expr_values.index)

    records = []
    for row in top_candidates.head(30).itertuples():
        genes = split_genes(getattr(row, "target_genes", ""))
        matches = [g for g in genes if g in expr_index]
        if matches:
            values = expr_values.loc[matches].mean(axis=0)
            pct = float((values >= global_median).mean())
        else:
            pct = np.nan
        surv = _survival_signal(values if matches else None, clin)
        records.append(
            {
                "rank": int(row.rank),
                "canonical_drug_id": row.canonical_drug_id,
                "drug_name": row.drug_name,
                "target_genes": getattr(row, "target_genes", ""),
                "target_match_genes": ";".join(matches),
                "target_gene_match_rate": float(len(matches) / len(genes)) if genes else np.nan,
                "target_expression_pct": pct,
                "target_expressed": bool(np.isfinite(pct) and pct >= 0.30),
                "survival_p_value": surv["p_value"],
                "survival_direction": surv["direction"],
                "known_thyroid_control": str(row.drug_name).upper() in known,
                "ensemble_score": row.ensemble_score,
                "mean_pred_ln_ic50": row.mean_pred_ln_ic50,
            }
        )
    validation = pd.DataFrame(records)
    write_table(validation[["canonical_drug_id", "drug_name", "target_genes", "target_match_genes", "target_expression_pct", "target_expressed"]], out_dir / "thyroid_target_expression.csv")
    write_table(validation[["canonical_drug_id", "drug_name", "survival_p_value", "survival_direction"]], out_dir / "thyroid_survival_validation.csv")
    precision = _known_precision(top_candidates, known)
    write_table(pd.DataFrame([precision]), out_dir / "thyroid_known_drug_precision.csv")
    write_table(validation.head(15), out_dir / "top15_validated.csv")
    qc = {
        "step": "step8_external_validation",
        "status": "completed",
        "expression_gene_count": int(expr_values.shape[0]),
        "patient_count": int(expr_values.shape[1]),
        "target_gene_match_rate_mean": float(validation["target_gene_match_rate"].mean(skipna=True)),
        "target_expressed_count_top15": int(validation.head(15)["target_expressed"].sum()),
        "clinical_available": clin is not None,
        "topk_precision": precision,
    }
    write_json(paths.reports_dir / "qc_step8_external_validation.json", qc)
    return validation.head(15)


def _survival_signal(values: pd.Series | None, clinical: pd.DataFrame | None) -> dict[str, Any]:
    if values is None or clinical is None or clinical.empty:
        return {"p_value": np.nan, "direction": "no_clinical_data"}
    endpoint = "OS_MONTHS" if "OS_MONTHS" in clinical.columns else None
    patient_col = "PATIENT_ID" if "PATIENT_ID" in clinical.columns else clinical.columns[0]
    if endpoint is None:
        return {"p_value": np.nan, "direction": "no_endpoint"}
    merged = clinical[[patient_col, endpoint]].copy()
    merged[endpoint] = pd.to_numeric(merged[endpoint], errors="coerce")
    expr = values.rename("target_expression").reset_index().rename(columns={"index": patient_col})
    merged = merged.merge(expr, on=patient_col, how="inner").dropna()
    if len(merged) < 8:
        return {"p_value": np.nan, "direction": "too_few_patients"}
    high = merged.loc[merged["target_expression"] >= merged["target_expression"].median(), endpoint]
    low = merged.loc[merged["target_expression"] < merged["target_expression"].median(), endpoint]
    if len(high) < 3 or len(low) < 3:
        return {"p_value": np.nan, "direction": "too_few_patients"}
    p = float(stats.mannwhitneyu(high, low, alternative="two-sided").pvalue)
    direction = "high_expression_longer_survival" if high.median() > low.median() else "high_expression_shorter_survival"
    return {"p_value": p, "direction": direction}


def _known_precision(top_candidates: pd.DataFrame, known_upper: set[str]) -> dict[str, float]:
    out = {}
    names = top_candidates["drug_name"].astype(str).str.upper().tolist()
    for k in [5, 10, 20]:
        topk = names[: min(k, len(names))]
        out[f"p_at_{k}"] = float(sum(name in known_upper for name in topk) / len(topk)) if topk else np.nan
    return out


def admet_assessment(cfg: dict[str, Any], candidates: pd.DataFrame | None = None) -> pd.DataFrame:
    paths = pipeline_paths(cfg)
    if candidates is None:
        candidates = read_table(paths.external_validation_dir / "top15_validated.csv")
    if candidates is None:
        candidates = read_table(paths.results_dir / "ensemble" / "thyroid_ensemble_top30_drugs.csv")
    if candidates is None:
        raise FileNotFoundError("Candidate file missing; run external_validation first")

    rows = read_table(paths.processed_dir / "row_metadata.csv")
    drug_meta = rows[["canonical_drug_id", "drug_name", "canonical_smiles", "target_genes", "PATHWAY_NAME_NORMALIZED", "classification"]].drop_duplicates("canonical_drug_id") if rows is not None else candidates
    candidates = candidates.merge(drug_meta, on=["canonical_drug_id", "drug_name"], how="left", suffixes=("", "_meta"))
    if "canonical_smiles" not in candidates.columns and "canonical_smiles_meta" in candidates.columns:
        candidates["canonical_smiles"] = candidates["canonical_smiles_meta"]

    assay_files = sorted(paths.admet_source_dir.glob(cfg["admet"].get("assay_glob", "*.csv")))
    assay_cache = [_load_assay(path, cfg) for path in assay_files]
    details = []
    summary_rows = []
    for row in candidates.head(15).itertuples():
        smiles = getattr(row, "canonical_smiles", "")
        can, ok = canonicalize_smiles(smiles)
        candidate_fp = _fingerprint(can, cfg) if ok else None
        no_match = 0
        toxic_flags = []
        assay_hits = 0
        for assay_name, assay_df in assay_cache:
            match = _nearest_admet_match(candidate_fp, assay_df, cfg)
            if match["match_type"] == "no_match":
                no_match += 1
            else:
                assay_hits += 1
            if assay_name.lower() in {"ames", "dili", "herg"} and match.get("predicted_label") == 1:
                toxic_flags.append(assay_name)
            details.append(
                {
                    "canonical_drug_id": row.canonical_drug_id,
                    "drug_name": row.drug_name,
                    "assay": assay_name,
                    **match,
                }
            )
        coverage = float(assay_hits / len(assay_cache)) if assay_cache else 0.0
        if not ok:
            admet_category = "NO_SMILES"
        elif toxic_flags:
            admet_category = "Caution"
        elif coverage >= 0.5:
            admet_category = "Approved"
        else:
            admet_category = "Candidate"
        summary_rows.append(
            {
                "canonical_drug_id": row.canonical_drug_id,
                "drug_name": row.drug_name,
                "canonical_smiles": smiles,
                "target_genes": getattr(row, "target_genes", getattr(row, "target_genes_meta", "")),
                "PATHWAY_NAME_NORMALIZED": getattr(row, "PATHWAY_NAME_NORMALIZED", getattr(row, "PATHWAY_NAME_NORMALIZED_meta", "")),
                "classification": getattr(row, "classification", getattr(row, "classification_meta", "")),
                "ensemble_score": getattr(row, "ensemble_score", np.nan),
                "target_expressed": getattr(row, "target_expressed", False),
                "admet_coverage": coverage,
                "admet_no_match_assays": no_match,
                "toxicity_flags": ";".join(toxic_flags),
                "admet_category": admet_category,
            }
        )
    detail_df = pd.DataFrame(details)
    summary = pd.DataFrame(summary_rows)
    write_table(detail_df, paths.admet_output_dir / "admet_detailed_candidates.csv")
    write_table(summary, paths.admet_output_dir / "final_drug_candidates.csv")
    qc = {
        "step": "step9_admet",
        "candidate_count": int(len(summary)),
        "assay_count": int(len(assay_cache)),
        "category_counts": summary["admet_category"].value_counts(dropna=False).to_dict() if not summary.empty else {},
        "toxicity_flag_count": int(summary["toxicity_flags"].astype(str).str.len().gt(0).sum()) if not summary.empty else 0,
    }
    write_json(paths.admet_output_dir / "admet_summary.json", qc)
    write_json(paths.reports_dir / "qc_step9_admet.json", qc)
    return summary


def _fingerprint(smiles: str, cfg: dict[str, Any]) -> Any | None:
    if Chem is None or AllChem is None or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(
        mol,
        int(cfg["admet"].get("fingerprint_radius", 2)),
        nBits=int(cfg["admet"].get("fingerprint_bits", 2048)),
    )


def _load_assay(path: Path, cfg: dict[str, Any]) -> tuple[str, pd.DataFrame]:
    df = read_table(path)
    if df is None:
        return path.stem, pd.DataFrame()
    smiles_col = "smiles" if "smiles" in df.columns else "SMILES" if "SMILES" in df.columns else df.columns[0]
    label_col = "label" if "label" in df.columns else "Y" if "Y" in df.columns else df.columns[-1]
    rows = []
    for r in df[[smiles_col, label_col]].dropna().itertuples(index=False):
        can, ok = canonicalize_smiles(r[0])
        fp = _fingerprint(can, cfg) if ok else None
        if fp is not None:
            rows.append({"smiles": can, "label": r[1], "fp": fp})
    return path.stem, pd.DataFrame(rows)


def _nearest_admet_match(candidate_fp: Any | None, assay_df: pd.DataFrame, cfg: dict[str, Any]) -> dict[str, Any]:
    if candidate_fp is None:
        return {"match_type": "NO_SMILES", "similarity": np.nan, "predicted_label": np.nan}
    if assay_df.empty or "fp" not in assay_df.columns:
        return {"match_type": "no_match", "similarity": np.nan, "predicted_label": np.nan}
    sims = [float(DataStructs.TanimotoSimilarity(candidate_fp, fp)) for fp in assay_df["fp"]]
    best_idx = int(np.argmax(sims))
    best = sims[best_idx]
    if best >= float(cfg["admet"]["exact_threshold"]):
        match_type = "exact"
    elif best >= float(cfg["admet"]["close_threshold"]):
        match_type = "close_analog"
    elif best >= float(cfg["admet"]["analog_threshold"]):
        match_type = "analog"
    else:
        match_type = "no_match"
    label = assay_df.iloc[best_idx]["label"]
    try:
        label = int(float(label))
    except Exception:
        pass
    return {"match_type": match_type, "similarity": best, "predicted_label": label}


def knowledge_validation(cfg: dict[str, Any], candidates: pd.DataFrame | None = None) -> pd.DataFrame:
    paths = pipeline_paths(cfg)
    if candidates is None:
        candidates = read_table(paths.admet_output_dir / "final_drug_candidates.csv")
    if candidates is None:
        raise FileNotFoundError("ADMET final candidates missing")
    known = {x.upper() for x in cfg["known_thyroid_drugs"]}
    terms = [t.upper() for t in cfg["thyroid_biology_terms"]]
    records = []
    for row in candidates.itertuples():
        target_genes = str(getattr(row, "target_genes", ""))
        pathway = str(getattr(row, "PATHWAY_NAME_NORMALIZED", ""))
        text = f"{target_genes} {pathway}".upper()
        drug_name = str(row.drug_name)
        drug_target = 1.0 if split_genes(target_genes) else 0.0
        target_relevance = 1.0 if any(term in text for term in terms) else 0.0
        clinical = 1.0 if drug_name.upper() in known or str(getattr(row, "classification", "")).lower() == "approved" else 0.4
        mechanism = 1.0 if target_relevance else 0.3
        safety = 1.0 if getattr(row, "admet_category", "") == "Approved" else 0.6 if getattr(row, "admet_category", "") == "Candidate" else 0.2
        target_expr_bonus = 0.5 if bool(getattr(row, "target_expressed", False)) else 0.0
        score = drug_target + target_relevance + clinical + mechanism + safety + target_expr_bonus
        if score >= 4.2 and getattr(row, "admet_category", "") != "Caution":
            tier = "Tier 1"
        elif score >= 3.0:
            tier = "Tier 2"
        elif getattr(row, "admet_category", "") in {"NO_SMILES", "Caution"}:
            tier = "Excluded"
        else:
            tier = "Tier 3"
        if drug_name.upper() in known:
            final_category = "Known thyroid cancer positive control"
        elif str(getattr(row, "classification", "")).lower() in {"approved", "indication_expansion"}:
            final_category = "Thyroid indication expansion candidate"
        else:
            final_category = "True repurposing or exploratory candidate"
        records.append(
            {
                "drug_name": drug_name,
                "canonical_drug_id": row.canonical_drug_id,
                "target_genes": target_genes,
                "pathway": pathway,
                "ensemble_score": getattr(row, "ensemble_score", np.nan),
                "admet_category": getattr(row, "admet_category", ""),
                "drug_target_evidence": drug_target,
                "target_thyroid_relevance": target_relevance,
                "clinical_evidence": clinical,
                "mechanism_rationale": mechanism,
                "safety_profile": safety,
                "target_expression_bonus": target_expr_bonus,
                "knowledge_score": score,
                "tier": tier,
                "final_category": final_category,
            }
        )
    result = pd.DataFrame(records)
    tier_order = {"Tier 1": 1, "Tier 2": 2, "Tier 3": 3, "Excluded": 4}
    result["_tier_order"] = result["tier"].map(tier_order).fillna(99)
    result = result.sort_values(["_tier_order", "knowledge_score", "ensemble_score"], ascending=[True, False, False]).drop(columns=["_tier_order"])
    write_table(result, paths.knowledge_validation_dir / "validation_summary.csv")
    write_json(paths.knowledge_validation_dir / "thyroid_knowledge_validation_results.json", result.to_dict(orient="records"))
    write_table(result, paths.phase5_dir / "final_comprehensive_candidates.csv")
    write_table(result.loc[result["tier"].eq("Tier 1")], paths.phase5_dir / "tier1_high_confidence.csv")
    qc = {
        "step": "step10_knowledge_validation",
        "candidate_count": int(len(result)),
        "tier_counts": result["tier"].value_counts(dropna=False).to_dict(),
        "category_counts": result["final_category"].value_counts(dropna=False).to_dict(),
    }
    write_json(paths.reports_dir / "qc_step10_knowledge_validation.json", qc)
    return result


def build_final_report(cfg: dict[str, Any]) -> Path:
    paths = pipeline_paths(cfg)
    metric_tables = []
    for input_name in ["numeric", "numeric_smiles", "numeric_strong_context_smiles"]:
        table = read_table(paths.results_dir / f"{input_name}_metrics_summary.csv")
        if table is not None and not table.empty:
            table = table.copy()
            table["input_set"] = input_name
            metric_tables.append(table)
    metrics = pd.concat(metric_tables, ignore_index=True) if metric_tables else pd.DataFrame()
    ensemble = read_json(paths.results_dir / "ensemble" / "thyroid_ensemble_results.json") if (paths.results_dir / "ensemble" / "thyroid_ensemble_results.json").exists() else {}
    final = read_table(paths.phase5_dir / "final_comprehensive_candidates.csv")
    qcs = {}
    for qc_path in sorted(paths.reports_dir.glob("qc_step*.json")):
        qcs[qc_path.stem] = read_json(qc_path)

    md_lines = [
        "# 갑상선암 약물 재창출 Hybrid Pipeline 최종 리포트",
        "",
        f"- 프로젝트: {cfg['project']['name']}",
        f"- 질환: {cfg['project']['disease']} ({cfg['project']['tcga_code']})",
        "- 기준: Choi numeric base + say2 SMILES/strong context/random3 ensemble + Choi external/ADMET/KG validation",
        "- 최종 추천은 screened drug-response 데이터에 한정한다.",
        "",
        "## 1. 실행 요약",
        "",
        f"- Numeric policy: {cfg['features']['numeric_base_policy']}",
        f"- SMILES dimension: {cfg['features']['smiles_svd_dim']}",
        f"- Strong context dimension: {cfg['features']['strong_context_dim']}",
        f"- CV: random sample {cfg['training']['n_splits']}-fold OOF",
        "",
    ]
    if not metrics.empty:
        md_lines += ["## 2. 모델 성능", "", metrics_to_markdown(metrics), ""]
    if ensemble:
        md_lines += [
            "## 3. 앙상블 및 Diversity",
            "",
            f"- 앙상블 OOF Spearman: {ensemble.get('ensemble_metrics', {}).get('spearman')}",
            f"- 앙상블 OOF RMSE: {ensemble.get('ensemble_metrics', {}).get('rmse')}",
            f"- 모델 가중치: {ensemble.get('weights', {})}",
            "",
        ]
    if final is not None and not final.empty:
        md_lines += ["## 4. 최종 후보", "", dataframe_to_markdown(final.head(20)), ""]
    md_lines += [
        "## 5. QC 요약",
        "",
        "- 각 단계별 상세 QC는 `reports/qc_step*.json`에 저장했다.",
        "- 외부검증 결과는 `external_validation/`, ADMET 결과는 `admet/`, knowledge validation 결과는 `phase5_final_results/`에 저장했다.",
        "",
        "## 6. 산출물 경로",
        "",
        "- `data/X_numeric.npy`",
        "- `data/X_numeric_smiles.npy`",
        "- `data/X_numeric_strong_context_smiles.npy`",
        "- `results/random3/`",
        "- `results/ensemble/thyroid_ensemble_top30_drugs.csv`",
        "- `external_validation/top15_validated.csv`",
        "- `admet/final_drug_candidates.csv`",
        "- `phase5_final_results/final_comprehensive_candidates.csv`",
        "",
    ]
    md = "\n".join(md_lines)
    report_md = paths.phase5_dir / "FINAL_REPORT.md"
    report_md.write_text(md, encoding="utf-8")

    report_html = paths.phase5_dir / "FINAL_REPORT.html"
    report_html.write_text(markdownish_to_html(md, qcs), encoding="utf-8")
    return report_md


def metrics_to_markdown(metrics: pd.DataFrame) -> str:
    cols = ["input_set", "model", "spearman", "pearson", "rmse", "mae", "r2", "ndcg_at_20", "train_oof_spearman_gap"]
    available = [c for c in cols if c in metrics.columns]
    table = metrics[available].copy()
    for col in table.select_dtypes(include=[float]).columns:
        table[col] = table[col].map(lambda x: round(float(x), 4) if pd.notna(x) else "")
    return dataframe_to_markdown(table)


def dataframe_to_markdown(table: pd.DataFrame) -> str:
    if table.empty:
        return ""
    rendered = []
    headers = [str(c) for c in table.columns]
    rendered.append("| " + " | ".join(headers) + " |")
    rendered.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in table.itertuples(index=False):
        rendered.append("| " + " | ".join(html.escape(str(v)) for v in row) + " |")
    return "\n".join(rendered)


def markdownish_to_html(md: str, qcs: dict[str, Any]) -> str:
    body = []
    for line in md.splitlines():
        if line.startswith("# "):
            body.append(f"<h1>{html.escape(line[2:])}</h1>")
        elif line.startswith("## "):
            body.append(f"<h2>{html.escape(line[3:])}</h2>")
        elif line.startswith("- "):
            body.append(f"<p class='bullet'>{html.escape(line)}</p>")
        elif line.startswith("|"):
            body.append(f"<pre>{html.escape(line)}</pre>")
        elif line.strip():
            body.append(f"<p>{html.escape(line)}</p>")
        else:
            body.append("")
    qc_block = html.escape(json.dumps(_json_safe(qcs), ensure_ascii=False, indent=2))
    return f"""<!doctype html>
<html lang="ko">
<head>
  <meta charset="utf-8">
  <title>Thyroid Hybrid Pipeline Final Report</title>
  <style>
    body {{ margin: 0; padding: 40px; background: #07111f; color: #e8f1ff; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
    main {{ max-width: 1180px; margin: 0 auto; }}
    h1, h2 {{ color: #6ee7ff; }}
    p, pre {{ color: #c7d4e6; line-height: 1.6; }}
    pre {{ background: #101d31; border: 1px solid #223550; border-radius: 12px; padding: 12px; overflow-x: auto; }}
    .bullet {{ margin-left: 12px; }}
    .card {{ background: #0d1b2e; border: 1px solid #223550; border-radius: 18px; padding: 20px; margin: 20px 0; }}
  </style>
</head>
<body>
  <main>
    <div class="card">{''.join(body)}</div>
    <h2>QC JSON Snapshot</h2>
    <pre>{qc_block}</pre>
  </main>
</body>
</html>
"""


def write_simple_html_qc(path: Path, title: str, qc: dict[str, Any]) -> None:
    escaped = html.escape(json.dumps(_json_safe(qc), ensure_ascii=False, indent=2))
    path.write_text(
        f"<!doctype html><html lang='ko'><head><meta charset='utf-8'><title>{html.escape(title)}</title></head><body><h1>{html.escape(title)}</h1><pre>{escaped}</pre></body></html>",
        encoding="utf-8",
    )


def run_pipeline(config_path: str | Path, demo: bool = False, steps: list[str] | None = None) -> dict[str, Any]:
    cfg = load_config(config_path)
    paths = pipeline_paths(cfg)
    if demo:
        seed_demo_data(config_path)

    steps = steps or [
        "prepare",
        "numeric",
        "smiles",
        "context",
        "assemble",
        "train",
        "groupcv",
        "ensemble",
        "external",
        "admet",
        "knowledge",
        "report",
    ]
    results: dict[str, Any] = {}
    rows = None
    if "prepare" in steps:
        rows = prepare_response_subset(cfg)
        results["prepare_rows"] = int(len(rows))
    if "numeric" in steps:
        build_numeric_base(cfg, rows)
    if "smiles" in steps:
        build_smiles_svd(cfg, rows)
    if "context" in steps:
        build_strong_context(cfg, rows)
    input_paths = None
    if "assemble" in steps:
        input_paths = assemble_model_inputs(cfg)
        results["input_paths"] = {k: str(v) for k, v in input_paths.items()}
    if "train" in steps:
        summaries = train_models(cfg, input_paths)
        results["metrics"] = {k: v.to_dict(orient="records") for k, v in summaries.items()}
    if "groupcv" in steps:
        results["groupcv"] = groupcv_stress_test(cfg)
    top = None
    if "ensemble" in steps:
        top, diversity = ensemble_and_diversity(cfg)
        results["top_candidates"] = top.head(10).to_dict(orient="records")
        results["diversity_rows"] = int(len(diversity))
    validated = None
    if "external" in steps:
        validated = external_validation(cfg, top)
        results["validated_rows"] = int(len(validated))
    admet = None
    if "admet" in steps:
        admet = admet_assessment(cfg, validated)
        results["admet_rows"] = int(len(admet))
    if "knowledge" in steps:
        final = knowledge_validation(cfg, admet)
        results["final_rows"] = int(len(final))
    if "report" in steps:
        report = build_final_report(cfg)
        results["final_report"] = str(report)
    write_json(paths.reports_dir / "run_summary.json", results)
    return results


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run thyroid cancer hybrid drug repurposing pipeline")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--demo", action="store_true", help="Seed demo data before running")
    parser.add_argument("--steps", nargs="*", default=None)
    args = parser.parse_args(argv)
    result = run_pipeline(args.config, demo=args.demo, steps=args.steps)
    print(json.dumps(_json_safe(result), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
