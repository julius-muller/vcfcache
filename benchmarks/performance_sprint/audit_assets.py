"""Read-only verification of frozen annotation assets, with durable audit evidence."""

import argparse
import json
import time
from pathlib import Path

from benchmarks.run_pilot import sha256sum, write_json_atomic


def main():
    """Hash existing resources; missing/mismatched assets are never repaired in place."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    rows = []
    with (args.output / "assets.jsonl").open("x") as log:
        for manifest in args.manifest:
            for line in manifest.read_text().splitlines():
                if not line.strip():
                    continue
                expected, name = line.split(maxsplit=1)
                path = Path(name.lstrip("*"))
                row = dict(
                    path=str(path),
                    expected_sha256=expected,
                    manifest=str(manifest),
                    checked_unix=time.time(),
                )
                try:
                    row["observed_sha256"] = sha256sum(path)
                    row["state"] = (
                        "match" if row["observed_sha256"] == expected else "mismatch"
                    )
                except OSError as error:
                    row.update(state="unavailable", error=str(error))
                rows.append(row)
                log.write(json.dumps(row) + "\n")
                log.flush()
    summary = dict(
        manifests={str(p): sha256sum(p) for p in args.manifest},
        finished_unix=time.time(),
        matched=sum(r["state"] == "match" for r in rows),
        mismatched=sum(r["state"] == "mismatch" for r in rows),
        unavailable=sum(r["state"] == "unavailable" for r in rows),
        note="Asset identity evidence only: inputs, recipe/CPU settings and retained output evidence require separate reuse approval.",
    )
    write_json_atomic(args.output / "summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
