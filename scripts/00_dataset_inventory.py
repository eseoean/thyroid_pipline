#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import load_config, pipeline_paths, write_json, write_table


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def inspect_local(path: Path) -> dict:
    if path.is_dir():
        files = [p for p in path.rglob("*") if p.is_file()]
        return {
            "exists": True,
            "type": "directory",
            "file_count": len(files),
            "total_bytes": sum(p.stat().st_size for p in files),
            "examples": [str(p) for p in files[:10]],
        }
    if path.exists():
        return {
            "exists": True,
            "type": "file",
            "file_count": 1,
            "total_bytes": path.stat().st_size,
            "examples": [str(path)],
        }
    return {"exists": False, "type": "missing", "file_count": 0, "total_bytes": 0, "examples": []}


def inspect_s3(prefix: str) -> dict:
    if not prefix:
        return {"checked": False}
    try:
        proc = subprocess.run(["aws", "s3", "ls", prefix, "--recursive", "--summarize"], check=False, capture_output=True, text=True)
    except FileNotFoundError:
        return {"checked": True, "status": "aws_cli_missing"}
    output = proc.stdout + proc.stderr
    total_objects = None
    total_size = None
    for line in output.splitlines():
        if "Total Objects:" in line:
            total_objects = int(line.split(":", 1)[1].strip())
        if "Total Size:" in line:
            total_size = int(line.split(":", 1)[1].strip())
    return {
        "checked": True,
        "returncode": proc.returncode,
        "total_objects": total_objects,
        "total_size": total_size,
        "preview": output.splitlines()[:20],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Inventory thyroid pipeline datasets")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--s3-prefix", default="", help="Optional S3 prefix to inspect, e.g. s3://say2-4team/Lung_raw/")
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = pipeline_paths(cfg)
    rows = []
    payload = {"sources": [], "s3": inspect_s3(args.s3_prefix)}
    for source in cfg.get("data_sources", []):
        dest = resolve(paths.root, source["destination"])
        info = inspect_local(dest)
        row = {
            "name": source["name"],
            "required": bool(source.get("required", False)),
            "destination": str(dest),
            "exists": info["exists"],
            "type": info["type"],
            "file_count": info["file_count"],
            "total_bytes": info["total_bytes"],
            "description": source.get("description", ""),
        }
        rows.append(row)
        payload["sources"].append({**row, "examples": info["examples"]})
    missing_required = [r["name"] for r in rows if r["required"] and not r["exists"]]
    payload["missing_required"] = missing_required
    payload["ready_for_pipeline"] = len(missing_required) == 0
    write_table(pd.DataFrame(rows), paths.reports_dir / "dataset_inventory.csv")
    write_json(paths.reports_dir / "dataset_inventory.json", payload)
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0 if payload["ready_for_pipeline"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

