"""Freeze 14 cached-only screening conditions against verified retained references."""

import argparse
import copy
import json
import time
from pathlib import Path

from benchmarks.engine_profile.run import claim, latest, read, verify
from benchmarks.run_pilot import sha256sum, write_json_atomic


def design(original, screen="initial"):
    """Screen each change separately, with opposite per-sample condition order."""
    groups = []
    for identifier in ("g00", "g01", "g03", "g04"):
        group = copy.deepcopy(
            next(g for g in original["groups"] if g["id"] == identifier)
        )
        if screen == "memory-preflight" and group["recipe"] != "veprs-core":
            continue
        assay = "wes" if group["recipe"] == "fastvep" else "wgs"
        template = next(
            c for c in group["cells"] if c["assay"] == assay and c["role"] == "measure"
        )
        candidates = [
            ("control-adaptive", "adaptive", "baseline", {}),
            (
                "indexed-misses",
                "adaptive",
                "profile",
                {"VCFCACHE_EXPERIMENT_INDEXED_MISSES": "1"},
            ),
        ]
        if screen == "trusted-hits":
            candidates = [
                ("control-adaptive", "adaptive", "baseline", {}),
                (
                    "trust-routed-hits",
                    "adaptive",
                    "profile",
                    {"VCFCACHE_EXPERIMENT_TRUST_ROUTED_HITS": "1"},
                ),
            ]
        elif screen == "native-route":
            candidates = [
                ("control-adaptive", "adaptive", "baseline", {}),
                (
                    "native-route",
                    "adaptive",
                    "profile",
                    {"VCFCACHE_EXPERIMENT_NATIVE_ROUTE": "1"},
                ),
            ]
        elif screen == "memory-preflight":
            candidates = [
                ("control-high-memory", "native", "baseline", {}),
                (
                    "memory-preflight",
                    "native",
                    "profile",
                    {"VCFCACHE_EXPERIMENT_MEMORY_PREFLIGHT": "1"},
                ),
            ]
        elif screen != "initial":
            raise ValueError("unknown screening phase")
        elif assay == "wes":
            candidates += [
                ("control-native", "native", "baseline", {}),
                (
                    "window-4096",
                    "native",
                    "profile",
                    {"VCFCACHE_EXPERIMENT_WINDOW_BP": "4096"},
                ),
                (
                    "window-1024",
                    "native",
                    "profile",
                    {"VCFCACHE_EXPERIMENT_WINDOW_BP": "1024"},
                ),
            ]
        if identifier in {"g03", "g04"}:
            candidates.reverse()
        group["cells"] = []
        for index, (label, engine, build, environment) in enumerate(candidates):
            cell = dict(
                template,
                id=f"{identifier}-s{index:03d}",
                engine=engine,
                build=build,
                purpose=label,
                profile=True,
                environment={
                    "VCFCACHE_EXPERIMENT_INDEXED_MISSES": "0",
                    "VCFCACHE_EXPERIMENT_WINDOW_BP": "16384",
                    "VCFCACHE_EXPERIMENT_TRUST_ROUTED_HITS": "0",
                    "VCFCACHE_EXPERIMENT_MEMORY_PREFLIGHT": "0",
                    "VCFCACHE_EXPERIMENT_NATIVE_ROUTE": "0",
                    **environment,
                },
            )
            if screen == "memory-preflight":
                cell["memory_mb"] = 8192
            group["cells"].append(cell)
        groups.append(group)
    assert (
        sum(len(g["cells"]) for g in groups)
        == {"initial": 14, "trusted-hits": 8, "memory-preflight": 4, "native-route": 8}[
            screen
        ]
    )
    return groups


def prepare(parent, root, request):
    """Reuse immutable assets and baselines; never reannotate direct samples."""
    original = verify(parent)
    manifest = copy.deepcopy(original)
    manifest.update(
        root=str(root).replace("/mnt/data/slurm-results/", "/results/", 1),
        created_unix=time.time(),
        groups=design(original, request.get("screen", "initial")),
        builds={
            "baseline": original["builds"]["profile"],
            "profile": request["profile"],
        },
        reference_campaign=original["root"],
        reference_manifest_sha256=sha256sum(parent / "manifest.json"),
        max_cells=36,
        planned_timed_cells={
            "initial": 14,
            "trusted-hits": 8,
            "memory-preflight": 4,
            "native-route": 8,
        }[request.get("screen", "initial")],
        deadline="2026-11-02T18:00:00+01:00",
        references_policy="Retained direct and cached outputs; checksum verification outside timing; no direct reruns",
        sprint_policy="Experimental switches only; no default change before measured gates",
    )
    manifest.pop("harness_root", None)
    manifest.pop("harness_files", None)
    for build in manifest["builds"].values():
        for relative, digest in build["files"].items():
            path = Path(build["root"]).joinpath(relative)
            if not Path("/results").exists():
                path = Path(
                    str(path).replace("/results/", "/mnt/data/slurm-results/", 1)
                )
            if sha256sum(path) != digest:
                raise RuntimeError(f"deployment drift: {path}")
    claim(root)
    for name in ("logs", "attempts", "references"):
        (root / name).mkdir()
    write_json_atomic(root / "manifest.json", manifest)
    (root / "READY").write_text(sha256sum(root / "manifest.json") + "\n")
    return manifest


def report(root):
    """Report comparisons without automatically promoting a single-run winner."""
    manifest = verify(root)
    rows = []
    for group in manifest["groups"]:
        controls = {}
        for cell in group["cells"]:
            state = latest(root, cell["id"]) or {}
            if (
                cell["purpose"].startswith("control-")
                and state.get("state") == "complete"
            ):
                controls[cell["engine"]] = state["metrics"]["wall_seconds"]
        for cell in group["cells"]:
            state = latest(root, cell["id"]) or {}
            wall = state.get("metrics", {}).get("wall_seconds")
            control = controls.get(cell["engine"])
            rows.append(
                dict(
                    cell=cell["id"],
                    sample=cell["sample"],
                    recipe=cell["recipe"],
                    condition=cell["purpose"],
                    state=state.get("state", "not_started"),
                    wall_seconds=wall,
                    control_seconds=control,
                    reduction=1 - wall / control if wall and control else None,
                    direct_pass=state.get("direct_pass"),
                    cached_reference_pass=state.get("instrumentation_equivalent"),
                )
            )
    print(json.dumps(rows, indent=2))


def main():
    """Prepare immutable pilot input or print its durable results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare", "report"])
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--parent", type=Path)
    parser.add_argument("--request", type=Path)
    args = parser.parse_args()
    if args.command == "prepare":
        if not args.parent or not args.request:
            parser.error("prepare requires --parent and --request")
        prepare(args.parent, args.campaign, read(args.request))
    else:
        report(args.campaign)


if __name__ == "__main__":
    main()
