"""Finish missing semantic audits from retained files; never rerun annotation."""

import argparse
import json
import os
import socket
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from benchmarks.pass_fusion_pilot.run import MODES, manifest
from benchmarks.run_pilot import sha256sum, strict_semantic_compare, write_json_atomic


def recover(root, identifier, source_commit):
    """Bind original timings and unchanged same-node outputs before recovery."""
    data = manifest(root)
    case = next(c for c in data["cases"] if c["id"] == identifier)
    folder = root / "cases" / identifier
    state_path = folder / "state.json"
    original_hash = sha256sum(state_path)
    state = json.loads(state_path.read_text())
    assert state["node"] == case["node"] == socket.gethostname()
    assert state["state"] == "running" and len(state["rows"]) == 2 * len(MODES)
    log = root / "logs" / f"{identifier}-{state['job_id']}.out"
    assert "DUE TO TIME LIMIT" in log.read_text()
    assert not (folder / "AUDIT_RECOVERY.json").exists()
    assert len(os.sched_getaffinity(0)) == 8
    task = next(t for t in data["tasks"] if t["sample"] == case["sample"])
    bindings = []
    outputs = {}
    missing = []
    for row in state["rows"]:
        repeat, mode = row["repeat"], row["mode"]
        directory = folder / f"r{repeat:02d}" / mode
        command = json.loads((directory / "command.json").read_text())
        metrics = json.loads((directory / "resources.json").read_text())
        stats = json.loads((directory / "statistics.json").read_text())
        assert metrics["exit_code"] == 0 and len(metrics["affinity"]) == 8
        assert metrics["wall_seconds"] == row["seconds"]
        assert stats["variant_counts"]["total_output"] == task["input"]["records"]
        assert command["source_commit"] == data["commit"]
        bindings.append(
            dict(
                repeat=repeat,
                mode=mode,
                command_sha256=sha256sum(directory / "command.json"),
                resources_sha256=sha256sum(directory / "resources.json"),
            )
        )
        argv = command["command"]
        assert argv.count("-o") == 1
        output = Path(argv[argv.index("-o") + 1])
        if repeat == 2:
            scratch = Path(command["parameters"]["temp_dir"]).parent.resolve()
            assert str(scratch).startswith("/mnt/data/tmp/pass-fusion-")
            assert output.resolve(strict=True).is_relative_to(scratch)
            index = Path(str(output) + ".csi")
            assert index.exists()
            outputs[mode] = dict(
                path=str(output),
                sha256=sha256sum(output),
                index_sha256=sha256sum(index),
            )
        if mode == "direct":
            continue
        comparison = directory / "against_direct.json"
        if comparison.exists():
            value = json.loads(comparison.read_text())
            assert value["semantic_pass"] and value["ignored_fields"] == []
            assert value["records_compared"] == task["input"]["records"]
        else:
            assert repeat == 2, "earlier repeat output has already been cleaned"
            missing.append((mode, comparison))
    assert missing, "no missing comparisons to recover"
    with ProcessPoolExecutor(max_workers=2) as pool:
        pending = {
            pool.submit(
                strict_semantic_compare,
                Path(outputs[mode]["path"]),
                Path(outputs["direct"]["path"]),
            ): (mode, path)
            for mode, path in missing
        }
        for future in as_completed(pending):
            mode, path = pending[future]
            value = future.result()
            assert value["semantic_pass"] and value["ignored_fields"] == []
            assert (
                value["records_compared"] == task["input"]["records"]
                and not path.exists()
            )
            write_json_atomic(path, value)
    for value in outputs.values():
        assert sha256sum(Path(value["path"])) == value["sha256"]
        assert sha256sum(Path(value["path"] + ".csi")) == value["index_sha256"]
    assert sha256sum(state_path) == original_hash
    comparisons = {
        str(p.relative_to(folder)): sha256sum(p)
        for p in folder.glob("r*/*/against_direct.json")
    }
    assert len(comparisons) == 8
    write_json_atomic(
        folder / "AUDIT_RECOVERY.json",
        dict(
            accepted=True,
            audit_only=True,
            no_timing_reruns=True,
            original_state_sha256=original_hash,
            source_manifest_sha256=sha256sum(root / "manifest.json"),
            audit_source_commit=source_commit,
            audit_job_id=os.environ.get("SLURM_JOB_ID"),
            node=case["node"],
            timed_bindings=bindings,
            retained_repeat2_outputs=outputs,
            comparison_sha256=comparisons,
            recovered_modes=[mode for mode, _ in missing],
            reason="wall limit during untimed audit after all ten timed invocations",
        ),
    )
    print(
        json.dumps(
            dict(
                accepted=True,
                recovered_modes=[mode for mode, _ in missing],
                timings_changed=False,
            )
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--case", required=True)
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args()
    recover(args.root, args.case, args.source_commit)
