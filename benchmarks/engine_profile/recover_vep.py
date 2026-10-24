"""Create a separate eight-cell VEP WGS recovery; never mutate the frozen parent."""

import argparse
import copy
import time
from pathlib import Path

from benchmarks.engine_profile.run import claim, controller, latest, verify
from benchmarks.run_pilot import sha256sum, write_json_atomic


def create(parent, root, harness, commit):
    """Retain original engines/fixtures and bind explicit lower-memory VEP settings."""
    original = verify(parent)
    manifest = copy.deepcopy(original)
    groups, previous = [], {}
    for group in manifest["groups"]:
        cells = [
            c for c in group["cells"] if c["recipe"] == "vep" and c["assay"] == "wgs"
        ]
        if not cells:
            continue
        for cell in cells:
            state = latest(parent, cell["id"])
            if not state or state["state"] not in {"failed", "skipped"}:
                raise ValueError(f"not an unsuccessful parent cell: {cell['id']}")
            previous[cell["id"]] = state
            cell["annotation_overrides"] = {"vep_forks": 4, "vep_buffer": 50000}
        group["cells"] = cells
        groups.append(group)
    if sum(len(g["cells"]) for g in groups) != 8:
        raise ValueError("expected exactly eight VEP WGS recovery cells")
    harness_path = controller(harness)
    files = {
        str(p.relative_to(harness_path)): sha256sum(p)
        for p in sorted(harness_path.rglob("*"))
        if p.is_file() and "__pycache__" not in p.parts
    }
    if not files or not (harness_path / "benchmarks/engine_profile/job.sh").is_file():
        raise ValueError("incomplete recovery harness")
    claim(root)
    for name in ("logs", "attempts", "references"):
        (root / name).mkdir()
    manifest.update(
        root=str(root).replace("/mnt/data/slurm-results/", "/results/", 1),
        created_unix=time.time(),
        groups=groups,
        max_cells=8,
        harness_root=harness,
        harness_files=files,
        harness_commit=commit,
        recovery=dict(
            parent=original["root"],
            parent_manifest_sha256=sha256sum(parent / "manifest.json"),
            parent_states=previous,
            reason="VEP direct WGS OOM; dependent cells skipped",
            comparability="Separate lower-memory configuration: four VEP forks and 50000-record buffer; 16 CPU engine allocation unchanged. Do not pool with original VEP timings.",
        ),
    )
    write_json_atomic(root / "manifest.json", manifest)
    (root / "READY").write_text(sha256sum(root / "manifest.json") + "\n")
    return manifest


def main():
    """Prepare the recovery manifest on the controller before explicit submission."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent", type=Path)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--harness", required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    create(args.parent.resolve(), args.campaign.resolve(), args.harness, args.commit)
    print(args.campaign)


if __name__ == "__main__":
    main()
