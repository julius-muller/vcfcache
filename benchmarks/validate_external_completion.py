#!/usr/bin/env python3
"""Fail closed unless every external-campaign task is publication complete."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

EXPECTED_STRATEGIES = {
    "gnomad_af_0.1",
    "gnomad_af_0.01",
    "cohort_3_genomes",
}


class MissingCompletedTaskError(RuntimeError):
    """A completed task is not yet visible through the shared filesystem."""


def validate(campaign_root: Path, phase: str) -> dict[str, object]:
    """Validate task coverage, uniqueness, and semantic equality."""
    manifest = campaign_root / "manifests" / f"{phase}.tsv"
    with manifest.open(newline="") as handle:
        tasks = list(csv.DictReader(handle, delimiter="\t"))
    if not tasks:
        raise RuntimeError(f"Empty task manifest: {manifest}")
    seen_samples: set[str] = set()
    seen_conditions: set[tuple[str, str]] = set()
    for task in tasks:
        task_id = int(task["task_id"])
        sample = task["sample"]
        if sample in seen_samples:
            raise RuntimeError(f"Repeated sample in publication phase: {sample}")
        seen_samples.add(sample)
        summary_path = (
            campaign_root
            / phase
            / "tasks"
            / f"task-{task_id}"
            / "external_summary.json"
        )
        if not summary_path.exists():
            raise MissingCompletedTaskError(f"Missing completed task: {summary_path}")
        summary = json.loads(summary_path.read_text())
        if summary["task"]["sample"] != sample:
            raise RuntimeError(f"Task/sample mismatch at {task_id}")
        rows = summary["rows"]
        strategies = {row["strategy"] for row in rows}
        if len(rows) != 3 or strategies != EXPECTED_STRATEGIES:
            raise RuntimeError(f"Incomplete strategies at task {task_id}: {strategies}")
        for row in rows:
            if row.get("semantic_pass") is not True:
                raise RuntimeError(f"Semantic failure at task {task_id}")
            key = (sample, row["strategy"])
            if key in seen_conditions:
                raise RuntimeError(f"Repeated sample/strategy result: {key}")
            seen_conditions.add(key)
    return {
        "campaign_root": str(campaign_root),
        "phase": phase,
        "samples": len(seen_samples),
        "strategy_results": len(seen_conditions),
        "semantic_passes": len(seen_conditions),
        "complete": True,
    }


def validate_with_visibility_retry(
    campaign_root: Path,
    phase: str,
    *,
    attempts: int = 1,
    poll_seconds: float = 5.0,
) -> dict[str, object]:
    """Retry only transient missing-file failures from the shared NFS mount."""
    if attempts < 1:
        raise ValueError("attempts must be at least one")
    if poll_seconds < 0:
        raise ValueError("poll_seconds must not be negative")
    for attempt in range(1, attempts + 1):
        try:
            return validate(campaign_root, phase)
        except MissingCompletedTaskError:
            if attempt == attempts:
                raise
            time.sleep(poll_seconds)
    raise AssertionError("visibility retry loop did not return")


def main() -> None:
    """Run command-line validation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--phase", default="measured")
    parser.add_argument("--visibility-attempts", type=int, default=1)
    parser.add_argument("--poll-seconds", type=float, default=5.0)
    args = parser.parse_args()
    print(
        json.dumps(
            validate_with_visibility_retry(
                args.campaign_root,
                args.phase,
                attempts=args.visibility_attempts,
                poll_seconds=args.poll_seconds,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
