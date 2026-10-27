# SPDX-License-Identifier: GPL-3.0-only
"""Opt-in resident native overlay; arbitrary annotation recipes remain Python-owned."""

import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

from vcfcache.database.engine_runtime import (
    Pipelines,
    cpu_affinity,
    helper_version,
    native_binary,
    staged_output,
    terminate_process,
    validate_info_headers,
)
from vcfcache.database.fused import run_fused
from vcfcache.database.memory_preflight import predict_buffer
from vcfcache.database.profiling import phase, profiled


@profiled("native.attempt")
def run_native(
    manager,
    cache,
    preserve=False,
    skip_split=False,
    workers=0,
    lookup="auto",
    memory_mb=512,
):
    """Keep buffered records in one child; no background service or format change."""
    if memory_mb < 1 or workers < 0:
        raise ValueError("native memory budget must be positive; workers non-negative")
    cpus = int(manager.params_file_content.get("threads", 1))
    if cpus < 1:
        raise ValueError("threads must be positive")
    with cpu_affinity(cpus) as selected:
        return _run(
            manager,
            cache,
            preserve,
            skip_split,
            workers,
            lookup,
            memory_mb,
            len(selected),
        )


def _run(manager, cache, preserve, skip_split, workers, lookup, memory_mb, cpus):
    start = time.perf_counter()
    binary = native_binary()
    version = helper_version(binary, (0, 3, 0))
    window_bp = int(os.environ.get("VCFCACHE_EXPERIMENT_WINDOW_BP", "16384"))
    if window_bp not in {1024, 4096, 16384}:
        raise ValueError("VCFCACHE_EXPERIMENT_WINDOW_BP must be 1024, 4096 or 16384")
    if window_bp != 16384:
        helper_version(binary, (0, 5, 0))
    if workers > cpus:
        raise ValueError("native workers exceed configured CPU budget")
    workers = workers or min(8, max(1, cpus // 2))
    validate_info_headers(manager.input_file, cache)
    bcftools = manager.params_file_content["bcftools_cmd"]
    with phase("native.index_preflight"):
        count = subprocess.run(
            [bcftools, "index", "-n", str(manager.input_file)],
            capture_output=True,
            text=True,
        )
    if count.returncode == 0 and int(count.stdout) * 1024 > memory_mb * 1024 * 1024:
        # Even the minimum retained-record allowance cannot fit. Avoid reading
        # part of a WGS only to discover a predictable buffer-budget fallback.
        result = run_fused(manager, cache, preserve, skip_split, 0, lookup)
        manager.last_run_stats.update(
            requested_engine="native",
            fallback_reason="native memory buffer budget exceeded (indexed record-count lower bound)",
        )
        return result
    prediction = None
    # The paired 36-cell screen confirmed that bounded preflight avoids wasted
    # resident loads on over-budget WGS inputs. Keep an explicit engineering
    # override for replaying the old control, but enable the supported guard.
    predict = os.environ.get("VCFCACHE_EXPERIMENT_MEMORY_PREFLIGHT", "1")
    if predict not in {"0", "1"}:
        raise ValueError("VCFCACHE_EXPERIMENT_MEMORY_PREFLIGHT must be 0 or 1")
    if (
        predict == "1"
        and memory_mb >= 2048
        and count.returncode == 0
        and int(count.stdout) * 1024 > memory_mb * 1024 * 1024 / 2
    ):
        prediction = predict_buffer(manager.input_file, cache, int(count.stdout))
        if (
            prediction
            and prediction["estimated_buffer_bytes"] > memory_mb * 1024 * 1024 * 1.1
        ):
            result = run_fused(manager, cache, preserve, skip_split, 0, lookup)
            manager.last_run_stats.update(
                requested_engine="native",
                memory_prediction=prediction,
                fallback_reason="native memory buffer predicted to exceed budget (bounded payload sample)",
            )
            return result
    root = manager.work_dir / "native"
    root.mkdir(parents=True, exist_ok=True)
    pipelines = Pipelines(manager.logger)
    quote = shlex.quote
    # Normalization writes a FIFO, not an indexed on-disk sample. The resident
    # process reads it once; run normalization asynchronously to avoid deadlock.
    fifo = root / "input.bcf.pipe"
    os.mkfifo(fifo)
    source = quote(str(manager.input_file))
    normalize = (
        f"{bcftools} view {source} -Ou"
        if skip_split
        else f"{bcftools} norm {source} -m- -Ou"
    )
    stderr = (root / "native.stderr").open("wb")
    normalizer = None
    process = None
    output = manager._resolve_output_bcf(manager.input_file.stem)
    target = root / "stdout.bcf" if output == "-" else Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        with staged_output(target, ".vcfcache-native-") as stage:
            staged = stage.path
            spec = dict(
                input=str(fifo),
                cache=str(cache),
                misses=str(root / "misses.bcf"),
                annotated=str(root / "annotated.bcf"),
                output=str(staged),
                stats=str(root / "stats.json"),
                tag=manager.nfa_config_content["must_contain_info_tag"],
                workers=workers,
                cpus=cpus,
                memory_bytes=memory_mb * 1024 * 1024,
                preserve=preserve,
                lookup=lookup,
                window_bp=window_bp,
            )
            job = root / "job.json"
            job.write_text(json.dumps(spec))
            normalizer = subprocess.Popen(
                [
                    "/bin/bash",
                    "-euo",
                    "pipefail",
                    "-c",
                    f"{normalize} > {quote(str(fifo))}",
                ],
                stderr=stderr,
                start_new_session=True,
            )
            process = subprocess.Popen(
                [binary, "resident", "--job", str(job)],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=stderr,
                text=True,
                start_new_session=True,
            )

            def receive():
                line = process.stdout.readline()
                if not line:
                    code = process.wait()
                    detail = (root / "native.stderr").read_text(errors="replace")[
                        -4000:
                    ]
                    raise RuntimeError(f"native engine stopped ({code}): {detail}")
                return json.loads(line)

            with phase("native.await_misses", child_pid=process.pid):
                state = receive()
            if state["event"] == "fallback":
                stage.aborted = True
                process.wait()
                # A budget fallback can stop reading the FIFO early (SIGPIPE).
                # Stop only our normalization group; never replay an annotator.
                terminate_process(normalizer)
                result = run_fused(manager, cache, preserve, skip_split, 0, lookup)
                manager.last_run_stats.update(
                    requested_engine="native", fallback_reason=state["reason"]
                )
                return result
            if state["event"] != "misses-ready":
                raise RuntimeError(f"unexpected native event: {state}")
            if normalizer.wait():
                raise RuntimeError("native normalization failed")
            if state["missing_variants"]:
                auxiliary = manager.output_dir / "auxiliary"
                auxiliary.mkdir(exist_ok=True)
                command = manager._substitute_variables(
                    manager.nfa_config_content["annotation_cmd"],
                    extra_vars={
                        "INPUT_BCF": spec["misses"],
                        "OUTPUT_BCF": spec["annotated"],
                        "AUXILIARY_DIR": str(auxiliary),
                    },
                )
                pipelines.run(command, root, "annotate-misses")
                for source_name, dest in (
                    ("stdout", "annotation_tool.log"),
                    ("stderr", "annotation_tool_err.log"),
                ):
                    shutil.copyfile(
                        root / f"annotate-misses.{source_name}",
                        manager.output_dir / dest,
                    )
                if not Path(spec["annotated"]).exists():
                    raise RuntimeError("annotation recipe did not produce OUTPUT_BCF")
                with phase("annotation.indexability_check"):
                    check = subprocess.run(
                        [bcftools, "index", "-f", spec["annotated"]],
                        capture_output=True,
                    )
                if check.returncode:
                    pipelines.run(
                        f"{bcftools} sort {quote(spec['annotated'])} -Ob -o sorted.bcf && {bcftools} index sorted.bcf",
                        root,
                        "sort-annotator-output",
                    )
                    spec["annotated"] = str(root / "sorted.bcf")
                validate_info_headers(Path(spec["misses"]), Path(spec["annotated"]))
            process.stdin.write(
                json.dumps({"action": "resume", "annotated": spec["annotated"]}) + "\n"
            )
            process.stdin.flush()
            state = receive()
            compatibility_reason = None
            if state["event"] == "merge-fallback":
                compatibility_reason = state["reason"]
                resolved = root / "merge-resolved.bcf"
                pipelines.run(
                    f"{bcftools} annotate -a {quote(spec['annotated'])} "
                    f"{quote(state['input'])} -c INFO -Ob0 -W -o {quote(str(resolved))}",
                    root,
                    "compatibility-merge",
                )
                process.stdin.write(
                    json.dumps({"action": "resume", "annotated": str(resolved)}) + "\n"
                )
                process.stdin.flush()
                state = receive()
            if state["event"] != "complete" or process.wait():
                raise RuntimeError(f"native completion failed: {state}")
        stats = json.loads((root / "stats.json").read_text())
        stats.update(
            engine="native",
            requested_engine="native",
            native_version=version,
            workers=workers,
            cpu_budget=cpus,
            cpu_affinity=sorted(os.sched_getaffinity(0)),
            engine_lookup=lookup,
            experimental_window_bp=window_bp,
            memory_prediction=prediction,
            engine_memory_mb=memory_mb,
            stage_timings=pipelines.timings,
            engine_wall_seconds=time.perf_counter() - start,
            dropped_variants=max(
                0, stats["missing_variants"] - stats["missing_annotated"]
            ),
        )
        manager.last_run_stats = stats
        if compatibility_reason:
            stats["compatibility_merge_reason"] = compatibility_reason
        if output == "-":
            with target.open("rb") as stream:
                shutil.copyfileobj(stream, sys.stdout.buffer)
        return subprocess.CompletedProcess([binary, "resident"], 0)
    finally:
        pipelines.stop()
        for child in (process, normalizer):
            terminate_process(child)
        if process is not None:
            process.stdin.close()
            process.stdout.close()
        stderr.close()
        fifo.unlink(missing_ok=True)
