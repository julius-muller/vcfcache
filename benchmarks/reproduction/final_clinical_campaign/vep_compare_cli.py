#!/usr/bin/env python3
"""Run the VEP semantic comparator outside a timed benchmark cell."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from benchmarks.run_pilot import semantic_compare


def main() -> None:
    """Run VEP semantic comparison and write its structured report."""
    parser = argparse.ArgumentParser()
    parser.add_argument("direct", type=Path)
    parser.add_argument("cached", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    report = semantic_compare(args.cached, args.direct)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not report["semantic_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
