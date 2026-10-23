#!/usr/bin/env python3
"""Capture resource counters inside one independently allocated benchmark step."""

import argparse
import json
import resource
import subprocess
import time
from pathlib import Path


def snapshot():
    """Read current cgroup counters; absence is explicit, never a fabricated zero."""
    memberships = Path("/proc/self/cgroup").read_text().splitlines()
    relative = next((line[3:] for line in memberships if line.startswith("0::")), None)
    if relative is None:
        return None
    root = Path("/sys/fs/cgroup") / relative.lstrip("/")
    return {name: (root / name).read_text().strip() for name in
            ("cpu.stat", "io.stat", "memory.peak") if (root / name).exists()}


def main():
    """Run a complete CLI and save wall/CPU/step-aggregate memory counters."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    before = snapshot()
    start = time.perf_counter()
    result = subprocess.run(command, check=False)
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    args.metrics.write_text(json.dumps(dict(wall_seconds=time.perf_counter()-start,
        cpu_seconds=usage.ru_utime+usage.ru_stime,largest_child_rss_kib=usage.ru_maxrss,
        cgroup_before=before,cgroup_after=snapshot(),exit_code=result.returncode),indent=2)+"\n")
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
