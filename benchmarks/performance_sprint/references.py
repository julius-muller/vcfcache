"""Bind historical manuscript rows to retained baselines without rerunning annotation."""

import argparse
import json
import time
from pathlib import Path

from benchmarks.engine_profile.run import controller, read
from benchmarks.performance_sprint.ledger import inventory, read_tsv
from benchmarks.run_pilot import REPO_ROOT, sha256sum, write_json_atomic

RAW = Path("benchmarks/results/raw/campaigns")
FULL52 = {
    "vep": "external-vep-light-full52-bbef41d7db8e",
    "fastvep": "external-fastvep-804963323e11-light-strict",
}
ASSAY = "external-assay-v3-independent-af001-warm"
SPECTRA = {
    "AF ≥ 10% cache": "wgs-pipeline-spectrum-4f854e5-kpgp00319",
    "AF ≥ 1% cache": "wgs-pipeline-spectrum-af001-cd87186-kpgp00319",
}


def request(repo):
    """Resolve every historical row, preserving separate recipes and exclusions."""
    full = {
        tool: {
            r["sample"]: r
            for r in read_tsv(repo / RAW / name / "manifests/measured.tsv")
        }
        for tool, name in FULL52.items()
    }
    assays = {
        (r["tool"], r["sample"], r["assay"]): r
        for r in read_tsv(repo / RAW / ASSAY / "manifests/benchmark_tasks.tsv")
    }
    spectra = {
        strategy: {
            r["pipeline"]: r
            for r in read_tsv(repo / RAW / name / "manifests/tasks.tsv")
        }
        for strategy, name in SPECTRA.items()
    }
    rows = inventory(repo)
    for row in rows:
        if row["family"] == "native_clinical_assays":
            # The frozen C2 runner explicitly removed scratch on EXIT. Keep the
            # digest limitation visible, rather than inventing full BCF proof.
            row.update(
                reference_kind="historical_annotation_key_digest_only",
                expected_direct_digest=row["historical_direct_digest"],
                cpu_limit=8,
                normalization="BED projection when applicable; norm -m-; remove ALT=*; norm --rm-dup exact",
                full_direct_output_retained=False,
            )
            continue
        if row["family"] == "full52_wgs_strategies" or (
            row["family"] == "vep_matched_assays" and row["assay"] == "wgs"
        ):
            task = full[row["tool"]][row["sample"]]
            campaign = FULL52[row["tool"]]
            relative = f"measured/tasks/task-{task['task_id']}/runs"
            direct = "uncached/pilot/*/*/uncached_r01/metrics.json"
            cached = f"{row['strategy']}/pilot/*/*/cached_r01/metrics.json"
            row["expected_input_sha256"] = task["input_sha256"]
        elif row["family"] == "vep_matched_assays":
            task = assays[(row["tool"], row["sample"], row["assay"])]
            campaign, relative = ASSAY, f"tasks/task-{task['task_id']}"
            direct, cached = (
                "pilot/*/*/uncached_r01/metrics.json",
                "pilot/*/*/cached_r01/metrics.json",
            )
        elif row["family"] == "synthetic_pipeline_spectrum":
            task = spectra[row["strategy"]][row["level"]]
            campaign = SPECTRA[row["strategy"]]
            relative = f"tasks/task-{task['task_id']}/run"
            direct, cached = (
                "pilot/*/*/uncached_r01/metrics.json",
                "pilot/*/*/cached_r01/metrics.json",
            )
            row["expected_input_sha256"] = task["input_sha256"]
            # The 1% load campaign intentionally reused the same six direct
            # measurements. There is no AF-specific extra direct invocation.
            if row["strategy"] == "AF ≥ 1% cache":
                row["direct_root"] = (
                    f"/results/campaigns/{SPECTRA['AF ≥ 10% cache']}/{relative}"
                )
        else:
            raise ValueError(f"unresolved family {row['family']}")
        root = f"/results/campaigns/{campaign}/{relative}"
        row.update(
            reference_kind="retained_complete_record_output",
            direct_root=row.get("direct_root", root),
            cached_root=root,
            direct_metrics_glob=direct,
            cached_metrics_glob=cached,
            original_task=task,
            cpu_limit=8,
        )
    return dict(
        schema=1,
        rows=rows,
        sources={
            str(p): sha256sum(repo / p)
            for p in sorted({Path(r["source"]) for r in rows})
        },
        note="Reference identity audit only. No annotation calls and no automatic baseline reuse approval.",
    )


def audit(request_file, output):
    """Hash each unique retained output once and checkpoint every result."""
    document = read(request_file)
    output.mkdir(parents=True, exist_ok=False)
    cache = {}
    rows = []

    def check(root, pattern, expected_seconds):
        key = (root, pattern)
        if key not in cache:
            matches = sorted(controller(root).glob(pattern))
            if len(matches) != 1:
                raise RuntimeError(
                    f"expected one reference, found {len(matches)}: {root}/{pattern}"
                )
            path = matches[0]
            metrics = read(path)
            for name, digest_key in [
                ("output.bcf", "output_sha256"),
                ("output.bcf.csi", "output_index_sha256"),
            ]:
                artifact = path.parent / name
                if sha256sum(artifact) != metrics[digest_key]:
                    raise RuntimeError(f"retained output drift: {artifact}")
            cache[key] = dict(
                metrics_path=str(path),
                metrics_sha256=sha256sum(path),
                output_path=str(path.parent / "output.bcf"),
                output_sha256=metrics["output_sha256"],
                output_index_sha256=metrics["output_index_sha256"],
                wall_seconds=metrics["wall_seconds"],
                output_records=metrics["output_records"],
                command=metrics["command"],
                input_vcf=metrics["input_vcf"],
            )
        result = cache[key]
        if abs(float(result["wall_seconds"]) - float(expected_seconds)) > 0.05:
            raise RuntimeError(
                f"reference timing does not match source row: {root}/{pattern}"
            )
        return result

    with (output / "references.jsonl").open("x") as checkpoint:
        for row in document["rows"]:
            result = dict(
                id=row["id"], family=row["family"], reference_kind=row["reference_kind"]
            )
            try:
                if row["reference_kind"] == "historical_annotation_key_digest_only":
                    result.update(
                        state="digest_only",
                        expected_direct_digest=row["expected_direct_digest"],
                    )
                else:
                    result["direct"] = check(
                        row["direct_root"],
                        row["direct_metrics_glob"],
                        row["historical_direct_seconds"],
                    )
                    result["cached"] = check(
                        row["cached_root"],
                        row["cached_metrics_glob"],
                        row["historical_cached_seconds"],
                    )
                    for role in ["direct", "cached"]:
                        command = result[role]["command"]
                        if ("--uncached" in command) != (role == "direct"):
                            raise RuntimeError("reference command role mismatch")
                    if result["direct"]["input_vcf"] != result["cached"]["input_vcf"]:
                        raise RuntimeError("direct/cached inputs differ")
                    result["state"] = "retained_outputs_verified"
            except (OSError, ValueError, KeyError, RuntimeError) as error:
                result.update(state="unresolved", error=str(error))
            rows.append(result)
            checkpoint.write(json.dumps(result) + "\n")
            checkpoint.flush()
    summary = dict(
        request_sha256=sha256sum(request_file),
        finished_unix=time.time(),
        rows=len(rows),
        retained_outputs_verified=sum(
            r["state"] == "retained_outputs_verified" for r in rows
        ),
        digest_only=sum(r["state"] == "digest_only" for r in rows),
        unresolved=sum(r["state"] == "unresolved" for r in rows),
        unique_output_pairs_hashed=len(cache),
        note="Output retention evidence only; CPU, preprocessing, recipes and resources still gate baseline reuse.",
    )
    write_json_atomic(output / "summary.json", summary)
    return summary


def main():
    """Prepare a request locally, or run its read-only audit on a worker."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare", "audit"])
    parser.add_argument("output", type=Path)
    parser.add_argument("--request", type=Path)
    args = parser.parse_args()
    if args.command == "prepare":
        if args.output.exists():
            raise FileExistsError(args.output)
        write_json_atomic(args.output, request(REPO_ROOT))
    else:
        if not args.request:
            parser.error("audit requires --request")
        print(json.dumps(audit(args.request, args.output), indent=2))


if __name__ == "__main__":
    main()
