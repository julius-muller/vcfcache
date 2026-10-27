"""Bind historical manuscript rows to retained baselines without rerunning annotation."""

import argparse
import json
import time
from pathlib import Path

from benchmarks.engine_profile.run import controller, read
from benchmarks.performance_sprint.ledger import EXTERNAL, inventory, read_tsv
from benchmarks.run_pilot import REPO_ROOT, sha256sum, write_json_atomic

RAW = Path("benchmarks/results/raw/campaigns")
ASSAY = "external-assay-v3-independent-af001-warm"
SPECTRA = {
    "AF ≥ 10% cache": "wgs-pipeline-spectrum-4f854e5-kpgp00319",
    "AF ≥ 1% cache": "wgs-pipeline-spectrum-af001-cd87186-kpgp00319",
}


def request(repo):
    """Resolve every historical row, preserving separate recipes and exclusions."""
    snapshot_path = EXTERNAL.parent / "SNAPSHOT.json"
    snapshot = read(repo / snapshot_path)
    if sha256sum(repo / EXTERNAL) != snapshot["files"][EXTERNAL.name]["sha256"]:
        raise RuntimeError("historical metrics do not match their source snapshot")
    campaigns = {tool: snapshot[tool]["campaign_id"] for tool in ("vep", "fastvep")}
    full = {
        tool: {
            r["sample"]: r
            for r in read_tsv(repo / RAW / name / "manifests/measured.tsv")
        }
        for tool, name in campaigns.items()
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
            campaign = campaigns[row["tool"]]
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
        if row["family"] == "synthetic_pipeline_spectrum":
            # These archives retain timings/comparisons, not the expensive raw
            # per-delay BCFs. Do not report them as complete-record references.
            baseline = spectra["AF ≥ 1% cache"][row["level"]]
            base = f"/results/campaigns/{SPECTRA['AF ≥ 10% cache']}/tasks/task-{task['task_id']}"
            current = f"/results/campaigns/{campaign}/tasks/task-{task['task_id']}"
            row.update(
                reference_kind="historical_semantic_metadata_only",
                baseline_summary=dict(
                    path=base + "/wgs_pipeline_summary.json",
                    sha256=baseline["baseline_summary_sha256"],
                ),
                baseline_comparison=base + "/semantic_comparison.json",
                cached_summary=current
                + (
                    "/wgs_pipeline_cached_summary.json"
                    if row["strategy"] == "AF ≥ 1% cache"
                    else "/wgs_pipeline_summary.json"
                ),
                composed_validation=f"/results/campaigns/{campaign}/publication/semantic_validation.json"
                if row["strategy"] == "AF ≥ 1% cache"
                else None,
            )
    return dict(
        schema=1,
        rows=rows,
        sources={
            str(p): sha256sum(repo / p)
            for p in sorted(
                {Path(r["source"]) for r in rows}
                | {snapshot_path}
                | {RAW / name / "manifests/measured.tsv" for name in campaigns.values()}
                | {RAW / ASSAY / "manifests/benchmark_tasks.tsv"}
                | {RAW / name / "manifests/tasks.tsv" for name in SPECTRA.values()}
            )
        },
        note="Reference identity audit only. No annotation calls and no automatic baseline reuse approval.",
    )


def historical_metadata(row):
    """Verify timing and semantic metadata, explicitly not missing raw outputs."""
    summary_path = controller(row["baseline_summary"]["path"])
    if sha256sum(summary_path) != row["baseline_summary"]["sha256"]:
        raise RuntimeError("historical synthetic baseline summary changed")
    baseline = read(summary_path)
    cached_path = controller(row["cached_summary"])
    cached = read(cached_path)
    comparison_path = controller(row["baseline_comparison"])
    comparison = read(comparison_path)
    if (
        baseline.get("semantic_pass") is not True
        or comparison.get("semantic_pass") is not True
    ):
        raise RuntimeError("historical synthetic semantic comparison did not pass")
    for value, expected in [
        (baseline["uncached_wall_seconds"], row["historical_direct_seconds"]),
        (cached["cached_wall_seconds"], row["historical_cached_seconds"]),
    ]:
        if abs(float(value) - float(expected)) > 0.05:
            raise RuntimeError("synthetic metadata timing does not match source table")
    if int(cached["input_records"]) != int(row["input_records"]):
        raise RuntimeError(
            "synthetic metadata record count does not match source table"
        )
    paths = [summary_path, cached_path, comparison_path]
    fingerprint = None
    if row.get("composed_validation"):
        validation_path = controller(row["composed_validation"])
        validation = read(validation_path)
        fingerprint = cached["semantic_fingerprint"]
        if (
            validation.get("semantic_pass") is not True
            or cached.get("source_semantic_pass") is not True
            or fingerprint["semantic_sha256"] != validation["canonical_semantic_sha256"]
            or fingerprint["csq_header_sha256"] != validation["csq_header_sha256"]
            or int(fingerprint["records"]) != int(row["input_records"])
        ):
            raise RuntimeError("composed synthetic fingerprint validation differs")
        paths.append(validation_path)
    return dict(
        state="historical_semantic_metadata_only",
        evidence={str(p): sha256sum(p) for p in paths},
        fingerprint=fingerprint,
        full_output_validation=False,
        note="Historical timings and semantic evidence verified; per-delay direct BCFs are not retained. Full-record reuse is NOT approved by this check.",
    )


def audit(request_file, output, unresolved_from=None):
    """Hash each unique retained output once and checkpoint every result."""
    document = read(request_file)
    selected = None
    if unresolved_from:
        selected = {
            r["id"]
            for line in unresolved_from.read_text().splitlines()
            if (r := json.loads(line))["state"] == "unresolved"
        }
        if selected - {r["id"] for r in document["rows"]}:
            raise ValueError("previous unresolved rows do not match the new request")
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
            if selected is not None and row["id"] not in selected:
                continue
            result = dict(
                id=row["id"], family=row["family"], reference_kind=row["reference_kind"]
            )
            try:
                if row["reference_kind"] == "historical_annotation_key_digest_only":
                    result.update(
                        state="digest_only",
                        expected_direct_digest=row["expected_direct_digest"],
                    )
                elif row["reference_kind"] == "historical_semantic_metadata_only":
                    result.update(historical_metadata(row))
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
        historical_metadata_only=sum(
            r["state"] == "historical_semantic_metadata_only" for r in rows
        ),
        unresolved=sum(r["state"] == "unresolved" for r in rows),
        unique_output_pairs_hashed=len(cache),
        full_request_rows=len(document["rows"]),
        unresolved_checkpoint_sha256=sha256sum(unresolved_from)
        if unresolved_from
        else None,
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
    parser.add_argument("--only-unresolved-from", type=Path)
    args = parser.parse_args()
    if args.command == "prepare":
        if args.output.exists():
            raise FileExistsError(args.output)
        write_json_atomic(args.output, request(REPO_ROOT))
    else:
        if not args.request:
            parser.error("audit requires --request")
        print(
            json.dumps(
                audit(args.request, args.output, args.only_unresolved_from), indent=2
            )
        )


if __name__ == "__main__":
    main()
