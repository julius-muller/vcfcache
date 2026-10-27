# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2026 Julius Müller

"""Experimental adaptive routing with bounded coordinate-parallel pipelines.

bcftools owns normalization and INFO overlay semantics. The annotator runs once,
on the concatenated misses, irrespective of the number of routing workers.
"""

import concurrent.futures
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pysam

from vcfcache.database.adaptive import (
    plan_tasks,
    preprocess_groups,
    select_workers,
    write_regions,
)
from vcfcache.database.engine_runtime import (
    Pipelines,
    cpu_affinity,
    helper_version,
    native_binary,
    staged_output,
    validate_info_headers,
)
from vcfcache.database.profiling import profiled


@profiled("adaptive.attempt")
def run_fused(
    manager,
    cache,
    preserve_unannotated=False,
    skip_split_multiallelic=False,
    workers=0,
    lookup="auto",
):
    """Bound the calling thread and its newly spawned descendants, then restore.

    Linux affinity is thread-local. Existing unrelated application threads are
    untouched; our new worker threads and subprocesses inherit the chosen mask.
    This avoids unsafe preexec_fn hooks in a multithreaded Python orchestrator.
    """
    threads = int(manager.params_file_content.get("threads", 1))
    with cpu_affinity(threads):
        return _run_fused(
            manager,
            cache,
            preserve_unannotated,
            skip_split_multiallelic,
            workers,
            lookup,
        )


def _run_fused(
    manager, cache, preserve_unannotated, skip_split_multiallelic, workers, lookup
):
    """Route coordinate tasks concurrently, annotate once, and gather in order."""
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
    explicit_workers = workers
    workers = workers or min(8, max(1, cpu_budget // 4))
    # Native routing and bcftools each need a foreground CPU. Rust writes final
    # compressed shards directly; remaining slots go to its BGZF pool.
    compression_threads = max(0, cpu_budget // workers - 2)
    root = manager.work_dir / "fused"
    root.mkdir(parents=True, exist_ok=True)
    pipelines = Pipelines(manager.logger)
    version = helper_version(binary, (0, 2, 0))
    indexed_miss_spools = os.environ.get("VCFCACHE_EXPERIMENT_INDEXED_MISSES", "0")
    if indexed_miss_spools not in {"0", "1"}:
        raise ValueError("VCFCACHE_EXPERIMENT_INDEXED_MISSES must be 0 or 1")
    indexed_miss_spools = indexed_miss_spools == "1"
    if indexed_miss_spools:
        helper_version(binary, (0, 5, 0))
    manager.logger.info(
        f"Fused engine: {version}; preparation workers={workers}; CPU budget={cpu_budget}"
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
    if not explicit_workers:
        workers = select_workers(prepared, cpu_budget, lookup)
    tasks, strategies = plan_tasks(prepared, lengths, workers, lookup)
    workers = min(workers, len(tasks))
    compression_threads = max(0, cpu_budget // workers - 2)
    manager.logger.info(
        f"Adaptive lookup: workers={workers}; codec threads={compression_threads}; tasks={len(tasks)}; mode={lookup}"
    )
    shards = [root / f"shard-{i:04d}" for i in range(len(tasks))]
    for shard, task in zip(shards, tasks, strict=True):
        shard.mkdir(parents=True, exist_ok=True)
        write_regions(shard / "regions.bed", task["regions"])

    def route(shard, task):
        pipelines.run(
            f"{bcftools} annotate -R regions.bed --regions-overlap 0 -a {quote(str(cache))} "
            f"{quote(str(task['source']))} -c INFO -Ou --threads {compression_threads} | "
            f"{native} route --hits hits.bcf --misses misses.raw.bcf "
            f"--order order.bin --stats route.json --tag {quote(tag)}"
            + (" --indexed-misses" if indexed_miss_spools else ""),
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
    # A miss decoder and the final encoder coexist during assembly. Sharing
    # the codec budget avoids assigning a full pool to BOTH pipeline halves.
    assembly_decode_threads = min(
        compression_threads // 2, int(compression_threads * missing / max(1, total))
    )
    assembly_encode_threads = (
        (compression_threads - assembly_decode_threads)
        if missing
        else max(0, cpu_budget // workers - 1)
    )
    if not total:
        assembly_encode_threads = 0
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
            prepare = (
                ""
                if indexed_miss_spools
                else f"{bcftools} view misses.raw.bcf -Ob0 -W -o misses.bcf && "
            )
            original_misses = "misses.raw.bcf" if indexed_miss_spools else "misses.bcf"
            overlay = f"{bcftools} annotate -R regions.bed --regions-overlap 0 -a {quote(str(annotated))} {original_misses} -c INFO -Ou --threads {assembly_decode_threads} | "
            miss_source = "-"
        else:
            prepare, overlay, miss_source = "", "", "misses.raw.bcf"
        pipelines.run(
            f"{prepare}{overlay}{native} assemble --hits hits.bcf --misses {miss_source} "
            f"--order order.bin --stats assemble.json --tag {quote(tag)}{preserve} "
            f"--output output.bcf --threads {assembly_encode_threads}",
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
    with staged_output(target, ".vcfcache-fused-") as stage:
        staged = stage.path
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
    if output == "-":
        with target.open("rb") as stream:
            shutil.copyfileobj(stream, sys.stdout.buffer)
    manager.last_run_stats = {
        "engine": "fused",
        "experimental_indexed_misses": indexed_miss_spools,
        "requested_engine": "fused",
        "native_version": version,
        "workers": workers,
        "cpu_budget": cpu_budget,
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "codec_threads_per_worker": compression_threads,
        "assembly_decode_threads": assembly_decode_threads,
        "assembly_encode_threads": assembly_encode_threads,
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
