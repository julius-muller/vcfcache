"""One timed Slurm step, with durable aggregate resource samples."""

import argparse
import json
import os
import resource
import signal
import subprocess
import time
from pathlib import Path

from benchmarks.run_pilot import cgroup_v2_snapshot, write_json_atomic


def main():
    """Keep sampling out of the CLI, including for profiling-disabled controls."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    start, epoch = time.monotonic(), time.time_ns()
    before = cgroup_v2_snapshot()
    process = subprocess.Popen(command, start_new_session=True)

    def stop(signum, frame):
        try:
            os.killpg(process.pid, signum)
        except ProcessLookupError:
            pass

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    with (args.directory / "resources.jsonl").open("w") as stream:
        while process.poll() is None:
            stream.write(
                json.dumps(dict(unix_ns=time.time_ns(), cgroup=cgroup_v2_snapshot()))
                + "\n"
            )
            stream.flush()
            try:
                process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                pass
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    write_json_atomic(
        args.directory / "resources.json",
        dict(
            start_unix_ns=epoch,
            wall_seconds=time.monotonic() - start,
            cpu_seconds=usage.ru_utime + usage.ru_stime,
            largest_child_rss_kib=usage.ru_maxrss,
            affinity=sorted(os.sched_getaffinity(0)),
            cgroup_before=before,
            cgroup_after=cgroup_v2_snapshot(),
            exit_code=process.returncode,
        ),
    )
    raise SystemExit(process.returncode)


if __name__ == "__main__":
    main()
