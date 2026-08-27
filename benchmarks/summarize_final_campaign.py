#!/usr/bin/env python3
"""Validate and summarize the final paired annotation campaign.

The input contains one direct and one cached row per sample/assay cell.  A cell
is eligible for timing summaries only when both rows are present, their record
counts agree, and the cached row passed the external content comparison.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


PAIR_FIELDS = ("sample", "cohort", "tool", "level", "assembly", "assay")
GROUP_FIELDS = ("tool", "level", "assembly", "assay")
SUMMARY_FIELDS = (
    *GROUP_FIELDS,
    "n_valid",
    "n_failed",
    "median_variants",
    "median_direct_seconds",
    "median_cached_seconds",
    "median_speedup",
    "median_hit_rate",
)
FAILURE_FIELDS = (*PAIR_FIELDS, "reason", "direct_digest", "cached_digest")


def _median(values: list[float]) -> float:
    return float(statistics.median(values))


def read_and_pair(path: Path) -> dict[tuple[str, ...], dict[str, dict[str, str]]]:
    pairs: dict[tuple[str, ...], dict[str, dict[str, str]]] = defaultdict(dict)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = set(PAIR_FIELDS) | {
            "condition",
            "wall_seconds",
            "input_records",
            "output_records",
            "missing_variants",
            "digest",
            "concordant",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"missing columns: {', '.join(sorted(missing))}")
        for line_number, row in enumerate(reader, start=2):
            condition = row["condition"]
            if condition not in {"direct", "cached"}:
                raise ValueError(f"line {line_number}: invalid condition {condition!r}")
            key = tuple(row[field] for field in PAIR_FIELDS)
            if condition in pairs[key]:
                raise ValueError(
                    f"line {line_number}: duplicate {condition} row for {key}"
                )
            pairs[key][condition] = row
    incomplete = [
        key for key, rows in pairs.items() if set(rows) != {"direct", "cached"}
    ]
    if incomplete:
        raise ValueError(
            f"{len(incomplete)} incomplete direct/cached pair(s): {incomplete[:3]}"
        )
    return dict(pairs)


def summarize(
    pairs: dict[tuple[str, ...], dict[str, dict[str, str]]]
) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    valid: dict[tuple[str, ...], list[dict[str, float]]] = defaultdict(list)
    failed_by_group: dict[tuple[str, ...], int] = defaultdict(int)
    failures: list[dict[str, str]] = []

    for pair_key, rows in pairs.items():
        direct, cached = rows["direct"], rows["cached"]
        reasons: list[str] = []
        direct_input = int(direct["input_records"])
        cached_input = int(cached["input_records"])
        direct_output = int(direct["output_records"])
        cached_output = int(cached["output_records"])
        if direct_input != cached_input:
            reasons.append("input_record_count")
        if direct_output != cached_output:
            reasons.append("output_record_count")
        if direct_input != direct_output or cached_input != cached_output:
            reasons.append("input_output_record_count")
        if cached["concordant"] != "PASS":
            reasons.append(cached["concordant"] or "missing_concordance_status")

        group_key = tuple(cached[field] for field in GROUP_FIELDS)
        if reasons:
            failed_by_group[group_key] += 1
            failure = {field: cached[field] for field in PAIR_FIELDS}
            failure.update(
                {
                    "reason": ";".join(reasons),
                    "direct_digest": direct["digest"],
                    "cached_digest": cached["digest"],
                }
            )
            failures.append(failure)
            continue

        direct_seconds = float(direct["wall_seconds"])
        cached_seconds = float(cached["wall_seconds"])
        missing_variants = int(cached["missing_variants"])
        if direct_seconds <= 0 or cached_seconds <= 0 or direct_input <= 0:
            raise ValueError(f"non-positive timing or input count for {pair_key}")
        if not 0 <= missing_variants <= direct_input:
            raise ValueError(f"invalid missing-variant count for {pair_key}")
        valid[group_key].append(
            {
                "variants": float(direct_input),
                "direct_seconds": direct_seconds,
                "cached_seconds": cached_seconds,
                "speedup": direct_seconds / cached_seconds,
                "hit_rate": 1.0 - missing_variants / direct_input,
            }
        )

    summaries: list[dict[str, object]] = []
    for group_key in sorted(set(valid) | set(failed_by_group)):
        cells = valid[group_key]
        if not cells:
            continue
        summary: dict[str, object] = dict(zip(GROUP_FIELDS, group_key, strict=True))
        summary.update(
            {
                "n_valid": len(cells),
                "n_failed": failed_by_group[group_key],
                "median_variants": _median([cell["variants"] for cell in cells]),
                "median_direct_seconds": _median(
                    [cell["direct_seconds"] for cell in cells]
                ),
                "median_cached_seconds": _median(
                    [cell["cached_seconds"] for cell in cells]
                ),
                "median_speedup": _median([cell["speedup"] for cell in cells]),
                "median_hit_rate": _median([cell["hit_rate"] for cell in cells]),
            }
        )
        summaries.append(summary)
    return summaries, failures


def write_tsv(
    path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("metrics", type=Path)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--failures", type=Path, required=True)
    args = parser.parse_args()

    pairs = read_and_pair(args.metrics)
    summaries, failures = summarize(pairs)
    write_tsv(args.summary, SUMMARY_FIELDS, summaries)
    write_tsv(args.failures, FAILURE_FIELDS, failures)
    print(
        f"validated {len(pairs)} paired cells: "
        f"{len(pairs) - len(failures)} accepted, {len(failures)} excluded"
    )


if __name__ == "__main__":
    main()
