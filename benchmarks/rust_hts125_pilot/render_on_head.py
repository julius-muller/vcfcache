"""Bounded one-shot R postprocessor on the controller, where R is installed."""

import argparse
import json
import os
import subprocess
import time
from pathlib import Path


def main():
    """Wait for terminal collected evidence; never submit or rerun benchmarks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    state_path = args.root / "render_state.json"
    state = dict(
        state="waiting",
        pid=os.getpid(),
        started=time.time(),
        deadline=time.time() + 5 * 3600,
    )
    state_path.write_text(json.dumps(state, indent=2) + "\n")
    while time.time() < state["deadline"]:
        review = args.root / "REVIEW.json"
        if review.exists():
            data = json.loads(review.read_text())
            if all(
                c["state"] in {"complete", "failed", "interrupted"}
                for c in data["cases"]
            ):
                script = Path(__file__).with_name("render.R")
                tools = json.loads((args.root / "manifest.json").read_text())["tools"]
                library = (
                    tools.replace("/results/", "/mnt/data/slurm-results/", 1)
                    + "/R-library"
                )
                with (args.root / "render_on_head.log").open("w") as log:
                    result = subprocess.run(
                        [
                            "Rscript",
                            "--vanilla",
                            str(script),
                            str(args.root),
                            "figures-head-v1",
                        ],
                        env=dict(os.environ, R_LIBS_USER=library),
                        stdout=log,
                        stderr=log,
                    )
                state.update(
                    state="complete" if result.returncode == 0 else "failed",
                    exit_code=result.returncode,
                    finished=time.time(),
                )
                state_path.write_text(json.dumps(state, indent=2) + "\n")
                return
        time.sleep(30)
    state.update(
        state="expired",
        finished=time.time(),
        reason="terminal review not ready within five hours",
    )
    state_path.write_text(json.dumps(state, indent=2) + "\n")


if __name__ == "__main__":
    main()
