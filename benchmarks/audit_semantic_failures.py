#!/usr/bin/env python3
"""Re-evaluate preserved failed Slurm pairs with the current comparator."""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from benchmarks.run_pilot import semantic_compare, write_json_atomic


def audit_one(report_path: str) -> dict[str, Any]:
    """Re-run comparison for one preserved failed attempt."""
    report = Path(report_path)
    run_root = report.parent
    original = json.loads(report.read_text())
    current = semantic_compare(
        run_root / "cached_r01/output.bcf",
        run_root / "uncached_r01/output.bcf",
    )
    return {
        "task_id": int(
            next(part for part in report.parts if part.startswith("task-")).split("-")[
                1
            ]
        ),
        "source_report": str(report),
        "original_annotation_mismatches": original.get("annotation_mismatches"),
        "rechecked_at": datetime.now(timezone.utc).isoformat(),
        "current_comparison": current,
    }


def find_reports(campaign: Path, phase: str) -> list[Path]:
    """Find one semantic report for each preserved failed attempt."""
    reports: list[Path] = []
    for failure in (campaign / phase / "attempts").glob("task-*/job-*/failure.json"):
        matches = list(
            failure.parent.glob("run/pilot/*/*/semantic_comparison_r01.json")
        )
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one semantic report beside {failure}: {matches}"
            )
        reports.append(matches[0])
    return sorted(
        reports,
        key=lambda path: int(
            next(part for part in path.parts if part.startswith("task-")).split("-")[1]
        ),
    )


def main() -> None:
    """Audit all failed pairs and retain a machine-readable summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True, type=Path)
    parser.add_argument("--phase", default="warmup")
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    campaign = args.campaign.expanduser().resolve()
    reports = find_reports(campaign, args.phase)
    if not reports:
        raise RuntimeError("No failed semantic reports found")
    output_dir = campaign / "emergency_review/hgnc_aware"
    output_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(audit_one, str(path)): path for path in reports}
        for future in as_completed(futures):
            value = future.result()
            write_json_atomic(output_dir / f"task-{value['task_id']}.json", value)
            results.append(value)
            print(
                json.dumps(
                    {
                        "task_id": value["task_id"],
                        "semantic_pass": value["current_comparison"]["semantic_pass"],
                        "remaining_mismatches": value["current_comparison"][
                            "annotation_mismatches"
                        ],
                        "ignored_hgnc_mismatches": value["current_comparison"][
                            "ignored_annotation_mismatches"
                        ],
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
    results.sort(key=lambda value: value["task_id"])
    summary = {
        "campaign": str(campaign),
        "phase": args.phase,
        "audited_at": datetime.now(timezone.utc).isoformat(),
        "failed_pairs_audited": len(results),
        "semantic_passes": sum(
            value["current_comparison"]["semantic_pass"] for value in results
        ),
        "remaining_failures": [
            value["task_id"]
            for value in results
            if not value["current_comparison"]["semantic_pass"]
        ],
        "tasks": results,
    }
    write_json_atomic(output_dir / "summary.json", summary)
    print(
        json.dumps(
            {
                key: summary[key]
                for key in (
                    "failed_pairs_audited",
                    "semantic_passes",
                    "remaining_failures",
                )
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
