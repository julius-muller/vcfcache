# SPDX-License-Identifier: GPL-3.0-only
"""Shared process, helper and publication contracts for optional routing engines."""

import os
import re
import shutil
import signal
import subprocess
import tempfile
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import pysam


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
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate)
    raise FileNotFoundError(
        "Adaptive/native routing requires vcfcache-engine. Build with "
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


@contextmanager
def cpu_affinity(threads):
    """Bound this thread and new descendants, then restore its original mask."""
    if threads < 1:
        raise ValueError("threads must be positive")
    if not hasattr(os, "sched_getaffinity") or not hasattr(os, "sched_setaffinity"):
        raise RuntimeError(
            "Adaptive/native engines require Linux affinity support; use --engine legacy"
        )
    previous = os.sched_getaffinity(0)
    selected = set(sorted(previous)[:threads])
    if not selected:
        raise RuntimeError("no CPUs available to annotation engine")
    try:
        os.sched_setaffinity(0, selected)
        yield selected
    finally:
        os.sched_setaffinity(0, previous)


def helper_version(binary, minimum):
    """Validate an explicitly resolved helper before starting any annotation."""
    version = subprocess.check_output([binary, "--version"], text=True).strip()
    match = re.fullmatch(r"vcfcache-engine (\d+)\.(\d+)\.(\d+)(?:[-+].*)?", version)
    if not match or tuple(map(int, match.groups())) < minimum:
        required = ".".join(map(str, minimum))
        raise RuntimeError(
            f"Routing requires vcfcache-engine >={required}; rebuild native/Cargo.toml or set VCFCACHE_ENGINE"
        )
    return version


def terminate_process(process):
    """Terminate only a process group created by our workflow, then reap it."""
    if process is None:
        return
    if process.poll() is None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


@contextmanager
def staged_output(target, prefix):
    """Publish completed BCF/CSI files beside their destination, on success only."""
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=prefix, dir=target.parent) as directory:
        stage = OutputStage(Path(directory) / "output.bcf")
        yield stage
        if stage.aborted:
            return
        staged = stage.path
        if not staged.is_file() or not Path(str(staged) + ".csi").is_file():
            raise RuntimeError("routing did not produce both BCF and CSI")
        os.replace(str(staged) + ".csi", str(target) + ".csi")
        os.replace(staged, target)


@dataclass
class OutputStage:
    """Explicitly distinguish an intentional fallback from a missing output."""

    path: Path
    aborted: bool = False


def select_engine(requested):
    """Auto chooses adaptive only when supported; explicit configuration errors fail."""
    if requested != "auto":
        return requested, None
    # Never conceal an invalid explicit override, even on an unsupported host.
    override = native_binary() if os.environ.get("VCFCACHE_ENGINE") else None
    if not hasattr(os, "sched_getaffinity") or not hasattr(os, "sched_setaffinity"):
        return "legacy", "optional routing engines unsupported on this platform"
    try:
        binary = override or native_binary()
    except FileNotFoundError:
        if os.environ.get("VCFCACHE_ENGINE"):
            raise
        return "legacy", "optional vcfcache-engine is not installed"
    helper_version(binary, (0, 2, 0))
    return (
        "adaptive",
        "optional helper available; adaptive selected (native remains opt-in)",
    )
