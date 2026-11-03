"""Approximate core-WGS command-wall attribution from existing logs, not new runs."""

import argparse
import csv
import hashlib
import re
import tarfile
from pathlib import Path


def components(text):
    """Take maxima for parallel shard commands; never sum concurrent walls."""
    stages = {i: [] for i in range(5)}
    stage = None
    for line in text.splitlines():
        match = re.search(r"Step ([0-4])/4:", line)
        if match:
            stage = int(match[1])
        match = re.search(r"Command completed in ([0-9.]+)s: (.*)", line)
        if match and stage is not None:
            stages[stage].append((match[2], float(match[1])))
    assert all(stages.values()), "incomplete log"

    def durations(stage, name):
        return [seconds for command, seconds in stages[stage] if command == name]

    # In step 2/4, shard filters are parallel; concat runs after all finish.
    # In step 4/4, shard merges precede one concat/index and final view/index.
    return dict(
        normalization_seconds=max(v for _, v in stages[0]),
        cache_overlay_seconds=max(durations(1, "bcftools annotate")),
        select_misses_seconds=max(durations(2, "bcftools filter"))
        + sum(durations(2, "bcftools concat")),
        annotate_misses_seconds=sum(v for _, v in stages[3]),
        merge_overlay_seconds=max(durations(4, "bcftools annotate")),
        gather_index_seconds=sum(durations(4, "bcftools concat")),
        final_write_index_seconds=sum(durations(4, "bcftools view")),
    )


def run(snapshot, output):
    """Export twelve auditable log-derived decompositions with residuals."""
    with (snapshot / "tables/cells.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    result = []
    for row in rows:
        if (row["tool"], row["level"], row["assay"]) != ("veprs", "core", "wgs"):
            continue
        relative = Path(row["source_evidence"]).relative_to("attempts")
        path = output / "source/logs" / relative / "run.log"
        parts = components(path.read_text())
        residual = float(row["cached_seconds"]) - sum(parts.values())
        assert residual >= -0.5, "command critical-path estimate exceeds total"
        result.append(
            dict(
                condition_id=row["condition_id"],
                sample=row["sample"],
                cohort=row["cohort"],
                direct_seconds=row["direct_seconds"],
                cached_seconds=row["cached_seconds"],
                **parts,
                unassigned_seconds=residual,
                log_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope="command walls; concurrent shard maxima; not instrumented phase timings",
            )
        )
    assert len(result) == 12
    path = output / "tables/core_WGS_log_headroom.tsv"
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, result[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result)
    with tarfile.open(output / "source/core_WGS_logs.tar.gz", "w:gz") as archive:
        for row in result:
            paths = list(
                (output / "source/logs").glob("*/" + row["condition_id"] + "/run.log")
            )
            assert len(paths) == 1
            archive.add(
                paths[0], arcname=str(paths[0].relative_to(output / "source/logs"))
            )
    import statistics

    for key in parts:
        print(key, round(statistics.median(r[key] for r in result), 2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.snapshot, args.output)
