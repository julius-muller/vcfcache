# SPDX-License-Identifier: GPL-3.0-only
"""Opt-in engineering traces; no public CLI or annotation-semantic changes."""

import functools
import json
import os
import threading
import time
from contextlib import contextmanager
from pathlib import Path

_lock = threading.Lock()


@contextmanager
def phase(name, **details):
    """Append a completed span, including failures, only when explicitly enabled."""
    directory = os.environ.get("VCFCACHE_PROFILE_DIR")
    if not directory:
        yield
        return
    start, epoch = time.monotonic_ns(), time.time_ns()
    outcome = "complete"
    try:
        yield
    except BaseException as error:
        outcome = type(error).__name__
        raise
    finally:
        duration = time.monotonic_ns() - start
        record = dict(
            schema_version=1,
            name=name,
            pid=os.getpid(),
            tid=threading.get_native_id(),
            start_unix_ns=epoch,
            duration_ns=duration,
            outcome=outcome,
            **details,
        )
        path = Path(directory) / f"python-{os.getpid()}.jsonl"
        # The harness creates the destination before timing. Missing trace
        # destinations fail visibly instead of silently losing profiling data.
        with _lock, path.open("a") as stream:
            stream.write(json.dumps(record) + "\n")


def profiled(name):
    """Trace a coarse function boundary, never the per-record hot path."""

    def decorate(function):
        @functools.wraps(function)
        def call(*args, **kwargs):
            with phase(name):
                return function(*args, **kwargs)

        return call

    return decorate
