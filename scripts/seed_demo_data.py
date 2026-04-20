#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thyroid_pipeline.core import seed_demo_data


def main() -> int:
    parser = argparse.ArgumentParser(description="Create deterministic demo thyroid pipeline inputs")
    parser.add_argument("--config", default="config/thyroid_pipeline_config.json")
    args = parser.parse_args()
    outputs = seed_demo_data(args.config)
    print(json.dumps(outputs, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

