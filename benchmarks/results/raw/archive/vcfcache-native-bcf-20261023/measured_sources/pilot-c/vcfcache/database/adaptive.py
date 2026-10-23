# SPDX-License-Identifier: GPL-3.0-only
"""Index-estimated lookup planning; no variant payloads are parsed in Python."""

import math
from pathlib import Path

WINDOW = 16384


def select_workers(prepared, cpu_budget, lookup):
    """Trade pipeline concurrency for BGZF decoding on cache-heavy work.

    Near one-to-one dense overlays spend relatively more time on target records;
    the original two-foreground-process allocation is appropriate there. Small
    inputs and larger candidate caches retain codec-priority allocation.
    """
    records = sum(
        row["records"] for _, metadata in prepared for row in metadata["windows"]
    )
    if records <= 2000:
        return 1
    candidates = 0
    for _, metadata in prepared:
        for contig in metadata["contigs"]:
            if "cache_records" not in contig:
                return min(8, max(1, cpu_budget // 4))
            full = contig["estimated_cache_bytes"]
            requested = contig["requested_cache_bytes"]
            sparse = lookup == "regions" or (lookup == "auto" and requested < full / 2)
            fraction = min(1, requested / full) if sparse and full else 1
            candidates += contig["cache_records"] * fraction
    divisor = 2 if records > 2000 and candidates <= records * 1.25 else 4
    return min(8, max(1, cpu_budget // divisor))


def plan_tasks(prepared, lengths, workers, lookup):
    """Return coordinate-ordered tasks, weighted by candidate cache bytes.

    `prepared` is a sequence of (normalized_path, native_metadata). All regions
    are zero-based, half-open internally and written as BED. Estimates are not
    measured I/O; CSI chunks and the final-block allowance are conservative.
    """
    selected = []
    for source, metadata in prepared:
        windows = metadata["windows"]
        for contig in metadata["contigs"]:
            sparse = lookup == "regions" or (
                lookup == "auto"
                and contig["requested_cache_bytes"]
                < contig["estimated_cache_bytes"] / 2
            )
            rows = [w for w in windows if w["contig"] == contig["name"]]
            selected.append((source, contig, rows, "regions" if sparse else "contigs"))
    work = sum(
        (
            contig["requested_cache_bytes"]
            if mode == "regions"
            else contig["estimated_cache_bytes"]
        )
        + sum(row["records"] for row in rows) * 64
        for _, contig, rows, mode in selected
    )
    # Bound process startup on tiny caches while retaining multiple tasks per
    # worker on expensive inputs. Task order is independent of completion order.
    target = max(1024 * 1024, work / max(1, workers * 2))
    if (
        workers == 1
        and sum(row["records"] for _, meta in prepared for row in meta["windows"])
        <= 2000
    ):
        target = math.inf
    units = []
    strategies = []
    for source, contig, rows, mode in selected:
        name = contig["name"]
        strategies.append({**contig, "strategy": mode})
        if mode == "regions":
            units.extend(
                (
                    source,
                    name,
                    row["start"],
                    row["end"],
                    row["estimated_cache_bytes"] + row["records"] * 64,
                )
                for row in rows
            )
        else:
            length = max(lengths.get(name) or 0, rows[-1]["end"])
            splits = min(
                max(1, math.ceil(contig["estimated_cache_bytes"] / target)), workers * 2
            )
            step = max(WINDOW, math.ceil(length / splits / WINDOW) * WINDOW)
            buckets = {}
            for row in rows:
                bucket = row["start"] // step
                buckets[bucket] = buckets.get(bucket, 0) + row["records"]
            for bucket, records in buckets.items():
                start, end = bucket * step, min((bucket + 1) * step, length)
                weight = (
                    contig["estimated_cache_bytes"] * (end - start) / length
                    + records * 64
                )
                units.append((source, name, start, end, weight))
    tasks = []
    current = None
    for source, name, start, end, weight in units:
        if (
            current is None
            or current["source"] != source
            or current["weight"] >= target
        ):
            current = {"source": source, "regions": [], "weight": 0}
            tasks.append(current)
        regions = current["regions"]
        if regions and regions[-1][0] == name and regions[-1][2] == start:
            regions[-1] = (name, regions[-1][1], end)
        else:
            regions.append((name, start, end))
        current["weight"] += weight
    if not tasks:
        # An empty target still needs the annotation header. Restrict overlay to
        # a single base instead of walking the cache to EOF. It emits no records.
        source = prepared[0][0]
        name = next(iter(lengths), "__vcfcache_empty__")
        tasks = [{"source": source, "regions": [(name, 0, 1)], "weight": 0}]
    return tasks, strategies


def write_regions(path: Path, regions):
    """BED extension makes coordinate conventions unambiguous to bcftools."""
    path.write_text(
        "".join(f"{name}\t{start}\t{end}\n" for name, start, end in regions)
    )


def preprocess_groups(contigs, indexed_counts, workers):
    """Preprocessing cost follows input counts, unlike annotation lookup cost."""
    if not contigs:
        return [None]
    if sum(indexed_counts.values()) <= 2000:
        return [",".join(contigs)]
    target = sum(indexed_counts.values()) / min(workers, len(contigs))
    groups, group, size = [], [], 0
    for name in contigs:
        if group and size >= target and len(groups) < workers - 1:
            groups.append(",".join(group))
            group, size = [], 0
        group.append(name)
        size += indexed_counts[name]
    if group:
        groups.append(",".join(group))
    return groups
