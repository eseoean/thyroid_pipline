#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import load_config, pipeline_paths


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def copy_local(src: Path, dest: Path, dry_run: bool) -> str:
    if dry_run:
        return f"DRY_RUN local copy {src} -> {dest}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        dest.mkdir(parents=True, exist_ok=True)
        for item in src.iterdir():
            target = dest / item.name
            if item.is_dir():
                shutil.copytree(item, target, dirs_exist_ok=True)
            else:
                shutil.copy2(item, target)
    else:
        shutil.copy2(src, dest)
    return f"copied {src} -> {dest}"


def copy_s3(s3_uri: str, dest: Path, dry_run: bool) -> str:
    if dry_run:
        return f"DRY_RUN aws s3 cp {s3_uri} {dest}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["aws", "s3", "cp", s3_uri, str(dest)]
    if s3_uri.endswith("/") or str(dest).endswith("/"):
        cmd.append("--recursive")
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr or proc.stdout)
    return proc.stdout.strip() or f"copied {s3_uri} -> {dest}"


def main() -> int:
    parser = argparse.ArgumentParser(description="Acquire datasets into local thyroid pipeline data directories")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    parser.add_argument("--manifest", default="config/data_manifest.json")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = pipeline_paths(cfg)
    manifest_path = resolve(paths.root, args.manifest)
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    logs = []
    for source in manifest.get("sources", []):
        dest = resolve(paths.root, source["destination"])
        local_path = source.get("local_path", "").strip()
        s3_uri = source.get("s3_uri", "").strip()
        if local_path:
            logs.append({"name": source["name"], "status": copy_local(Path(local_path).expanduser(), dest, args.dry_run)})
        elif s3_uri:
            logs.append({"name": source["name"], "status": copy_s3(s3_uri, dest, args.dry_run)})
        else:
            logs.append({"name": source["name"], "status": "skipped_no_source"})
    print(json.dumps(logs, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

