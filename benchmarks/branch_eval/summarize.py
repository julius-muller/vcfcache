#!/usr/bin/env python3
"""Summarize the branch evaluation: measured floors and two-point slopes.

The one-variant cell is reported as a measured floor, not as a fitted startup
constant: it is a complete invocation whose per-variant work is negligible, so
it bounds fixed cost from below without assuming the cost model is linear.
Slopes are differences between adjacent measured scales and are labelled as
such; a least-squares line over every available scale is reported beside them
so that a reader can see when the two disagree.
"""

import argparse
import json
from pathlib import Path

SCALES = ("single", "minimal", "panel", "wes", "wgs")
ARMS = ("direct", "main-legacy", "legacy", "original-fused", "auto", "native")
ARM_LABELS = {
    "direct": "uncached",
    "main-legacy": "main",
    "legacy": "legacy@head",
    "original-fused": "fused",
    "auto": "adaptive",
    "native": "native",
}


def load(roots):
    """Collect every timed cell, keyed by annotator, scale and arm."""
    cells, counts, audits = {}, {}, []
    for root in roots:
        for measurements in sorted(Path(root).rglob("measurements.json")):
            conditions = measurements.parent.parent / "CONDITIONS.txt"
            if not conditions.is_file():
                continue
            fields = dict(
                line.split("=", 1)
                for line in conditions.read_text().split()
                if "=" in line
            )
            annotator = fields["annotator"]
            for row in json.loads(measurements.read_text()):
                scale = Path(row["input"]).stem
                cells[(annotator, scale, row["arm"])] = row["wall_seconds"]
                audits.append((annotator, scale, row["arm"], row["semantic_pass"]))
    inventory = Path(roots[0]).parent / "inputs/INVENTORY.tsv"
    if inventory.is_file():
        for line in inventory.read_text().splitlines()[1:]:
            name, records, _ = line.split("\t")
            counts[name] = int(records)
    return cells, counts, audits


def fit(points):
    """Least-squares slope and intercept over (records, seconds) pairs."""
    if len(points) < 2:
        return None, None
    n = len(points)
    mean_x = sum(x for x, _ in points) / n
    mean_y = sum(y for _, y in points) / n
    spread = sum((x - mean_x) ** 2 for x, _ in points)
    if spread == 0:
        return None, None
    slope = sum((x - mean_x) * (y - mean_y) for x, y in points) / spread
    return slope, mean_y - slope * mean_x


def main():
    """Print exploratory summaries; see audit_engines.py for backend-aware data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", nargs="+")
    args = parser.parse_args()
    cells, counts, audits = load(args.roots)
    annotators = sorted({key[0] for key in cells})

    failed = [row for row in audits if not row[3]]
    print(f"# timed cells: {len(cells)}   semantic audits failed: {len(failed)}")
    for row in failed:
        print(f"#   FAILED {row[0]} {row[1]} {row[2]}")
    print()

    print("## Wall seconds per cell")
    header = ["annotator", "scale", "records"] + [ARM_LABELS[a] for a in ARMS]
    print("\t".join(header))
    for annotator in annotators:
        for scale in SCALES:
            if not any((annotator, scale, arm) in cells for arm in ARMS):
                continue
            row = [annotator, scale, str(counts.get(scale, ""))]
            for arm in ARMS:
                value = cells.get((annotator, scale, arm))
                row.append(f"{value:.2f}" if value is not None else "-")
            print("\t".join(row))
    print()

    print("## Fixed floor (one variant) and per-variant cost")
    print(
        "\t".join(
            [
                "annotator",
                "arm",
                "floor_1var_s",
                "slope_panel_wes_us",
                "slope_wes_wgs_us",
                "lsq_slope_us",
                "lsq_intercept_s",
            ]
        )
    )
    for annotator in annotators:
        for arm in ARMS:
            # Fit only the assay scales. The one- and 22-variant cells exist to
            # measure the fixed floor; including them lets a fixed-cost
            # pathology dominate the regression and produce slopes that are not
            # per-variant costs at all, including negative ones.
            available = [
                (counts[scale], cells[(annotator, scale, arm)])
                for scale in ("panel", "wes", "wgs")
                if scale in counts and (annotator, scale, arm) in cells
            ]
            if not available:
                continue
            floor = cells.get((annotator, "single", arm))
            pairs = {}
            for low, high in (("panel", "wes"), ("wes", "wgs")):
                if (annotator, low, arm) in cells and (annotator, high, arm) in cells:
                    delta_t = cells[(annotator, high, arm)] - cells[(annotator, low, arm)]
                    delta_n = counts[high] - counts[low]
                    pairs[(low, high)] = delta_t / delta_n * 1e6
            slope, intercept = fit(available)
            print(
                "\t".join(
                    [
                        annotator,
                        ARM_LABELS[arm],
                        f"{floor:.2f}" if floor is not None else "-",
                        f"{pairs.get(('panel', 'wes'), float('nan')):.1f}",
                        f"{pairs.get(('wes', 'wgs'), float('nan')):.1f}",
                        f"{slope * 1e6:.1f}" if slope is not None else "-",
                        f"{intercept:.2f}" if intercept is not None else "-",
                    ]
                )
            )
    print()

    print("## Speedup over uncached annotation")
    print("\t".join(["annotator", "scale"] + [ARM_LABELS[a] for a in ARMS[1:]]))
    for annotator in annotators:
        for scale in SCALES:
            direct = cells.get((annotator, scale, "direct"))
            if direct is None:
                continue
            row = [annotator, scale]
            for arm in ARMS[1:]:
                value = cells.get((annotator, scale, arm))
                row.append(f"{direct / value:.2f}x" if value else "-")
            print("\t".join(row))


if __name__ == "__main__":
    main()
