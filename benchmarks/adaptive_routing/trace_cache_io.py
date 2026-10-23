#!/usr/bin/env python3
"""Untimed diagnostic: count cache-file reads, not physical storage traffic."""

import argparse
import json
import re
import subprocess
from pathlib import Path


def main():
    """Trace one invocation; retain raw evidence and a small reproducible summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-bcf", type=Path, required=True)
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.trace.exists() or args.summary.exists():
        parser.error("trace and summary paths must be new")
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    result = subprocess.run(
        [
            "strace",
            "-f",
            "-qq",
            "-yy",
            "-s",
            "0",
            "-e",
            "trace=read,pread64,lseek",
            "-o",
            str(args.trace),
            *command,
        ],
        check=False,
    )
    marker = f"<{args.cache_bcf.resolve()}>"
    calls = reads = seeks = 0
    pending = set()
    with args.trace.open() as stream:
        for line in stream:
            tid = line.split(maxsplit=1)[0]
            if marker in line:
                if "lseek(" in line:
                    seeks += 1
                if re.search(r"\b(?:read|pread64)\(", line):
                    calls += 1
                    if "<unfinished ...>" in line:
                        pending.add(tid)
                    else:
                        count = re.search(r"= (\d+)\s*$", line)
                        if count:
                            reads += int(count[1])
            elif tid in pending and re.search(
                r"<\.\.\. (?:read|pread64) resumed>", line
            ):
                count = re.search(r"= (\d+)\s*$", line)
                if count:
                    reads += int(count[1])
                pending.discard(tid)
    args.summary.write_text(
        json.dumps(
            dict(
                command=command,
                cache=str(args.cache_bcf),
                exit_code=result.returncode,
                cache_read_calls=calls,
                cache_read_bytes=reads,
                cache_seek_calls=seeks,
                unresolved_reads=len(pending),
                note="strace diagnostic, not a timed benchmark; file-read bytes include page-cache hits and repeated reads, not decoded bytes or physical disk traffic",
            ),
            indent=2,
        )
        + "\n"
    )
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
