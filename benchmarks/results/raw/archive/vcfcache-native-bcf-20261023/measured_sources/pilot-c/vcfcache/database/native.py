# SPDX-License-Identifier: GPL-3.0-only
"""Opt-in resident native overlay; arbitrary annotation recipes remain Python-owned."""

import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from vcfcache.database.fused import (
    Pipelines,
    native_binary,
    run_fused,
    validate_info_headers,
)


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
    previous = os.sched_getaffinity(0)
    selected = set(sorted(previous)[:cpus])
    try:
        os.sched_setaffinity(0, selected)
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
    finally:
        os.sched_setaffinity(0, previous)


def _run(manager, cache, preserve, skip_split, workers, lookup, memory_mb, cpus):
    start = time.perf_counter()
    binary = native_binary()
    version = subprocess.check_output([binary, "--version"], text=True).strip()
    match = re.search(r"vcfcache-engine (\d+)\.(\d+)\.(\d+)", version)
    if not match or tuple(map(int, match.groups())) < (0, 3, 0):
        raise RuntimeError(
            "--engine native requires vcfcache-engine >=0.3.0; rebuild the helper"
        )
    if workers > cpus:
        raise ValueError("native workers exceed configured CPU budget")
    workers = workers or min(8, max(1, cpus // 2))
    validate_info_headers(manager.input_file, cache)
    bcftools = manager.params_file_content["bcftools_cmd"]
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
        with tempfile.TemporaryDirectory(
            prefix=".vcfcache-native-", dir=target.parent
        ) as staging:
            staged = Path(staging) / "output.bcf"
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

            state = receive()
            if state["event"] == "fallback":
                process.wait()
                # A budget fallback can stop reading the FIFO early (SIGPIPE).
                # Stop only our normalization group; never replay an annotator.
                if normalizer.poll() is None:
                    os.killpg(normalizer.pid, signal.SIGTERM)
                normalizer.wait()
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
                check = subprocess.run(
                    [bcftools, "index", "-f", spec["annotated"]], capture_output=True
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
            os.replace(str(staged) + ".csi", str(target) + ".csi")
            os.replace(staged, target)
        stats = json.loads((root / "stats.json").read_text())
        stats.update(
            engine="native",
            requested_engine="native",
            native_version=version,
            workers=workers,
            cpu_budget=cpus,
            cpu_affinity=sorted(os.sched_getaffinity(0)),
            engine_lookup=lookup,
            engine_memory_mb=memory_mb,
            stage_timings=pipelines.timings,
            engine_wall_seconds=time.perf_counter() - start,
            dropped_variants=stats["input_variants"] - stats["output_variants"],
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
            if child is not None:
                if child.poll() is None:
                    os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
        if process is not None:
            process.stdin.close()
            process.stdout.close()
        stderr.close()
        fifo.unlink(missing_ok=True)
