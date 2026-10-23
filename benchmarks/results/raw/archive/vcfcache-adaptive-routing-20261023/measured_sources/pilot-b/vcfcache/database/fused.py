# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2026 Julius Müller

"""Experimental native routing with bounded contig-parallel subprocess pipelines.

bcftools owns normalization and INFO overlay semantics. The annotator runs once,
on the concatenated misses, irrespective of the number of routing workers.
"""

import concurrent.futures
import json
import os
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import pysam

from vcfcache.database.adaptive import plan_tasks, preprocess_groups, write_regions


def native_binary() -> str:
    """Resolve the opt-in executable; never compile during annotation."""
    override = os.environ.get("VCFCACHE_ENGINE")
    if override:
        candidate = Path(override).expanduser().resolve()
        if not candidate.is_file() or not os.access(candidate, os.X_OK):
            raise FileNotFoundError(f"VCFCACHE_ENGINE is not executable: {candidate}")
        return str(candidate)
    installed = shutil.which("vcfcache-engine")
    if installed:
        return installed
    candidate = (
        Path(__file__).resolve().parents[2] / "native/target/release/vcfcache-engine"
    )
    if candidate.is_file():
        return str(candidate)
    raise FileNotFoundError(
        "--engine fused requires vcfcache-engine. Build with "
        "cargo build --release --locked --manifest-path native/Cargo.toml, "
        "or set VCFCACHE_ENGINE to the installed executable."
    )


def validate_info_headers(target, annotations):
    """Reject conflicting typed INFO definitions before an overlay can coerce data."""
    # Header-only reads do not need CSI indexes. Opening file objects avoids
    # HTSlib's misleading missing-index diagnostic on deliberately plain spools.
    with (
        Path(target).open("rb") as target_stream,
        Path(annotations).open("rb") as annotation_stream,
        pysam.VariantFile(target_stream) as original,
        pysam.VariantFile(annotation_stream) as incoming,
    ):
        for name in original.header.info:
            if name in incoming.header.info:
                left, right = original.header.info[name], incoming.header.info[name]
                if (left.type, left.number) != (right.type, right.number):
                    raise ValueError(
                        f"Incompatible INFO header for {name}: "
                        f"{left.type}/{left.number} versus {right.type}/{right.number}"
                    )


class Pipelines:
    """Track process groups and stop siblings on failure or interruption."""

    def __init__(self, logger):
        self.logger = logger
        self.lock = threading.Lock()
        self.active = set()
        self.cancelled = False
        self.timings = []

    def stop(self):
        with self.lock:
            self.cancelled = True
            for process in self.active:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass

    def run(self, command, directory, label):
        start = time.perf_counter()
        directory.mkdir(parents=True, exist_ok=True)
        (directory / f"{label}.sh").write_text(
            "#!/bin/bash\nset -euo pipefail\n" + command + "\n"
        )
        # Binary record streams remain inside the shell pipeline. Log files are
        # direct file descriptors, not pipes that can fill or capture a WGS.
        with (
            (directory / f"{label}.stdout").open("wb") as out,
            (directory / f"{label}.stderr").open("wb") as err,
        ):
            with self.lock:
                if self.cancelled:
                    raise RuntimeError("fused pipeline cancelled")
                process = subprocess.Popen(
                    ["/bin/bash", "-euo", "pipefail", "-c", command],
                    stdout=out,
                    stderr=err,
                    cwd=directory,
                    start_new_session=True,
                )
                self.active.add(process)
            try:
                code = process.wait()
                if code:
                    self.stop()
                    detail = (directory / f"{label}.stderr").read_text(
                        errors="replace"
                    )[-4000:]
                    raise RuntimeError(f"fused {label} failed ({code}): {detail}")
            finally:
                if process.poll() is None:
                    self.stop()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                with self.lock:
                    self.active.discard(process)
        duration = time.perf_counter() - start
        with self.lock:
            self.timings.append(
                {"stage": label, "task": directory.name, "wall_seconds": duration}
            )
        self.logger.info(f"Fused {label} completed in {duration:.3f}s")


def run_fused(
    manager,
    cache,
    preserve_unannotated=False,
    skip_split_multiallelic=False,
    workers=0,
    lookup="auto",
):
    """Route contigs concurrently, annotate once, and reconstruct original order."""
    started = time.perf_counter()
    phase_timings = []
    binary = native_binary()
    quote = shlex.quote
    bcftools = manager.params_file_content["bcftools_cmd"]
    native = quote(binary)
    tag = manager.nfa_config_content["must_contain_info_tag"]
    total_threads = int(manager.params_file_content.get("threads", 1))
    if total_threads < 1 or workers < 0:
        raise ValueError("threads must be positive and engine_workers non-negative")
    cpu_budget = min(total_threads, len(os.sched_getaffinity(0)))
    if lookup not in {"auto", "regions", "contigs"}:
        raise ValueError("engine_lookup must be auto, regions or contigs")
    validate_info_headers(manager.input_file, cache)
    with pysam.VariantFile(str(manager.input_file)) as variant_file:
        contigs = list(variant_file.header.contigs)
        lengths = {
            name: value.length for name, value in variant_file.header.contigs.items()
        }
    with pysam.VariantFile(str(cache)) as cache_file:
        for name, value in cache_file.header.contigs.items():
            lengths[name] = max(lengths.get(name) or 0, value.length or 0)
    # Keep the original header order; no alphabetical/numeric contig reordering.
    index = subprocess.run(
        [bcftools, "index", "--stats", str(manager.input_file)],
        capture_output=True,
        text=True,
        check=True,
    )
    indexed_counts = {
        line.split("\t")[0]: int(line.split("\t")[-1])
        for line in index.stdout.splitlines()
        if int(line.split("\t")[-1]) > 0
    }
    contigs = [name for name in contigs if name in indexed_counts]
    max_workers = max(1, cpu_budget // 2)
    if workers > max_workers:
        raise ValueError(
            f"engine_workers={workers} exceeds the {cpu_budget}-CPU pipeline budget; "
            f"use at most {max_workers} (two cooperating processes per worker)"
        )
    workers = workers or min(8, max(1, cpu_budget // 4))
    # Native routing and bcftools each need a foreground CPU. Rust writes final
    # compressed shards directly; remaining slots go to its BGZF pool.
    compression_threads = max(0, cpu_budget // workers - 2)
    root = manager.work_dir / "fused"
    root.mkdir(parents=True, exist_ok=True)
    pipelines = Pipelines(manager.logger)
    version = subprocess.run(
        [binary, "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    manager.logger.info(
        f"Fused engine: {version}; workers={workers}; CPU budget={cpu_budget}"
    )

    def parallel(function, items):
        phase_start = time.perf_counter()
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=workers)
        futures = {
            executor.submit(function, *item): index for index, item in enumerate(items)
        }
        try:
            results = [None] * len(items)
            for future in concurrent.futures.as_completed(futures):
                results[futures[future]] = future.result()
            return results
        except BaseException:
            pipelines.stop()
            for future in futures:
                future.cancel()
            raise
        finally:
            executor.shutdown(wait=True, cancel_futures=True)
            phase_timings.append(
                {
                    "phase": function.__name__,
                    "wall_seconds": time.perf_counter() - phase_start,
                }
            )

    def preprocess(shard, contig):
        shard.mkdir(parents=True, exist_ok=True)
        source = quote(str(manager.input_file))
        selection = ""
        if contig is not None:
            selection = f"-r {quote(contig)} --regions-overlap 0"
        if skip_split_multiallelic:
            preparation = f"{bcftools} view {selection} {source} "
        else:
            preparation = (
                f"{bcftools} norm {selection} {source} -m- -Ou | " f"{bcftools} view "
            )
        # Indexed target is required by bcftools annotate for BCF annotations.
        pipelines.run(
            f"{preparation}-e 'ALT=\"*\"' -Ou | "
            f"{native} prepare --output normalized.bcf --metadata plan.json --cache {quote(str(cache))}",
            shard,
            "preprocess",
        )
        return shard / "normalized.bcf", json.loads((shard / "plan.json").read_text())

    groups = preprocess_groups(contigs, indexed_counts, workers)
    prepared = parallel(
        preprocess,
        [(root / f"prepare-{i:04d}", group) for i, group in enumerate(groups)],
    )
    tasks, strategies = plan_tasks(prepared, lengths, workers, lookup)
    workers = min(workers, len(tasks))
    compression_threads = max(0, cpu_budget // workers - 2)
    shards = [root / f"shard-{i:04d}" for i in range(len(tasks))]
    for shard, task in zip(shards, tasks, strict=True):
        shard.mkdir(parents=True, exist_ok=True)
        write_regions(shard / "regions.bed", task["regions"])

    def route(shard, task):
        pipelines.run(
            f"{bcftools} annotate -R regions.bed --regions-overlap 0 -a {quote(str(cache))} "
            f"{quote(str(task['source']))} -c INFO -Ou --threads {compression_threads} | "
            f"{native} route --hits hits.bcf --misses misses.raw.bcf "
            f"--order order.bin --stats route.json --tag {quote(tag)}",
            shard,
            "route",
        )
        counts = json.loads((shard / "route.json").read_text())
        return counts

    counts = parallel(route, list(zip(shards, tasks, strict=True)))
    reasons = sorted(
        {value["fallback_reason"] for value in counts if value["fallback_reason"]}
    )
    if reasons:
        reason = "; ".join(reasons)
        manager.logger.warning(f"Fused engine falling back to legacy: {reason}")
        original_output = manager.output_file
        if original_output == "-":
            manager.output_file = root / "legacy-fallback.bcf"
        try:
            result = manager._run_annotate(
                cache, preserve_unannotated, skip_split_multiallelic
            )
            if original_output == "-":
                with manager.output_file.open("rb") as stream:
                    shutil.copyfileobj(stream, sys.stdout.buffer)
        finally:
            manager.output_file = original_output
        manager.last_run_stats.update(
            engine="legacy", requested_engine="fused", fallback_reason=reason
        )
        return result
    missing = sum(value["missing_variants"] for value in counts)
    total = sum(value["input_variants"] for value in counts)
    annotated = root / "annotated.bcf"
    indexed_misses = root / "misses.bcf"
    recipe_output_count = 0
    if missing:
        files = " ".join(quote(str(shard / "misses.raw.bcf")) for shard in shards)
        pipelines.run(
            f"{bcftools} concat --no-version {files} -Ob0 -W --threads {max(0, cpu_budget - 1)} "
            f"-o {quote(str(indexed_misses))}",
            root,
            "collect-misses",
        )
        auxiliary = manager.output_dir / "auxiliary"
        auxiliary.mkdir(exist_ok=True)
        command = manager._substitute_variables(
            manager.nfa_config_content["annotation_cmd"],
            extra_vars={
                "INPUT_BCF": str(indexed_misses),
                "OUTPUT_BCF": str(annotated),
                "AUXILIARY_DIR": str(auxiliary),
            },
        )
        pipelines.run(command, root, "annotate-misses")
        # Preserve existing recipe log locations for diagnostics and provenance.
        shutil.copyfile(
            root / "annotate-misses.stdout", manager.output_dir / "annotation_tool.log"
        )
        shutil.copyfile(
            root / "annotate-misses.stderr",
            manager.output_dir / "annotation_tool_err.log",
        )
        if not annotated.exists():
            raise RuntimeError("annotation recipe did not produce OUTPUT_BCF")
        # Index validates coordinate ordering. Only sort when it is not indexable.
        check = subprocess.run(
            [bcftools, "index", "-f", str(annotated)], capture_output=True
        )
        if check.returncode:
            pipelines.run(
                f"{bcftools} sort {quote(str(annotated))} -Ob -o sorted.bcf && "
                f"{bcftools} index sorted.bcf",
                root,
                "sort-annotator-output",
            )
            annotated = root / "sorted.bcf"
        recipe_output_count = int(
            subprocess.run(
                [bcftools, "index", "-n", str(annotated)],
                capture_output=True,
                text=True,
                check=True,
            ).stdout.strip()
        )
        validate_info_headers(shards[0] / "misses.raw.bcf", annotated)

    def assemble(shard, task):
        count = counts[shards.index(shard)]
        preserve = " --preserve-unannotated" if preserve_unannotated else ""
        if missing:
            # The overlay must operate on every original miss record, including
            # records dropped by the recipe. Never substitute recipe FORMAT/GT.
            prepare = f"{bcftools} view misses.raw.bcf -Ob0 -W -o misses.bcf && "
            overlay = f"{bcftools} annotate -R regions.bed --regions-overlap 0 -a {quote(str(annotated))} misses.bcf -c INFO -Ou --threads {compression_threads} | "
            miss_source = "-"
        else:
            prepare, overlay, miss_source = "", "", "misses.raw.bcf"
        pipelines.run(
            f"{prepare}{overlay}{native} assemble --hits hits.bcf --misses {miss_source} "
            f"--order order.bin --stats assemble.json --tag {quote(tag)}{preserve} "
            f"--output output.bcf --threads {compression_threads}",
            shard,
            "assemble",
        )
        result = json.loads((shard / "assemble.json").read_text())
        if result["input_variants"] != count["input_variants"]:
            raise RuntimeError("native routing/reconstruction count mismatch")
        return result

    assembled = parallel(assemble, list(zip(shards, tasks, strict=True)))
    output = manager._resolve_output_bcf(manager.input_file.stem)
    target = root / "final.bcf" if output == "-" else Path(output)
    # Stage beside the destination so publication uses same-filesystem rename.
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".vcfcache-fused-", dir=target.parent
    ) as staging:
        staged = Path(staging) / "output.bcf"
        files = " ".join(quote(str(shard / "output.bcf")) for shard in shards)
        # Naive concat checks header dictionary compatibility and copies BGZF
        # blocks without decompressing/recompressing the complete WGS again.
        pipelines.run(
            f"{bcftools} concat --no-version --naive {files} -Ob "
            f"-o {quote(str(staged))}",
            root,
            "gather",
        )
        pipelines.run(
            f"{bcftools} index --threads {max(0, cpu_budget - 1)} {quote(str(staged))}",
            root,
            "index",
        )
        os.replace(str(staged) + ".csi", str(target) + ".csi")
        os.replace(staged, target)
    if output == "-":
        with target.open("rb") as stream:
            shutil.copyfileobj(stream, sys.stdout.buffer)
    manager.last_run_stats = {
        "engine": "fused",
        "requested_engine": "fused",
        "native_version": version,
        "workers": workers,
        "cpu_budget": cpu_budget,
        "codec_threads_per_worker": compression_threads,
        "engine_lookup": lookup,
        "lookup_strategies": strategies,
        "lookup_tasks": len(tasks),
        "lookup_regions": sum(len(task["regions"]) for task in tasks),
        "task_estimated_work": [task["weight"] for task in tasks],
        "estimate_note": "CSI candidate byte estimates, not observed cache I/O",
        "phase_timings": phase_timings,
        "engine_wall_seconds": time.perf_counter() - started,
        "stage_timings": pipelines.timings,
        "input_variants": total,
        "missing_variants": missing,
        "output_variants": sum(value["output_variants"] for value in assembled),
        "annotated_output": sum(value["annotated_output"] for value in assembled),
        "missing_annotated": recipe_output_count,
        "dropped_variants": max(0, missing - recipe_output_count),
    }
    (root / "engine.json").write_text(
        json.dumps(manager.last_run_stats, indent=2) + "\n"
    )
    return subprocess.CompletedProcess([binary], 0)
