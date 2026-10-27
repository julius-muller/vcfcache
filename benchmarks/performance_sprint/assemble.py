"""Materialize the complete manuscript matrix from audited historical evidence.

Run on a worker with /results and the established annotator assets mounted.
No annotation is performed here. Baseline approvals bind explicit evidence and
retain digest-only or composed-reference limitations rather than inventing BCFs.
"""

import argparse
import copy
import csv
import hashlib
import json
import shlex
import shutil
import time
from pathlib import Path

import yaml

from benchmarks.engine_profile.run import read
from benchmarks.performance_sprint.manuscript import prepare, validate_task
from benchmarks.run_pilot import sha256sum, write_json_atomic


def rows(path):
    """Load flushed JSONL audit rows."""
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line]


def merge_audits(original, supplement):
    """Replace only explicitly unresolved rows; never silently drop coverage."""
    combined = {r["id"]: r for r in rows(original)}
    old_summary = Path(original).parent / "summary.json"
    old_time = read(old_summary)["finished_unix"] if old_summary.exists() else None
    for row in combined.values():
        row["audited_unix"] = old_time
    new_summary = Path(supplement).parent / "summary.json"
    new_time = read(new_summary)["finished_unix"] if new_summary.exists() else None
    for row in rows(supplement):
        if row["id"] not in combined or combined[row["id"]]["state"] != "unresolved":
            raise ValueError("supplement does not resolve an original unresolved row")
        combined[row["id"]] = row
        combined[row["id"]]["audited_unix"] = new_time
    if any(r["state"] == "unresolved" for r in combined.values()):
        raise ValueError("unresolved historical references remain")
    return combined


def argument(command, flag):
    """Read an exact frozen CLI argument; missing values are errors."""
    return command[command.index(flag) + 1]


def c2_input(row, grch38_inputs):
    """C2's frozen PGP subset is not necessarily the full52 PGP selection."""
    if row["assembly"] == "GRCh38":
        return grch38_inputs[row["sample"]]["input_vcf"]
    return f"/mnt/data/vcfcache_benchmarks/external_wgs/samples/{row['assembly']}/{row['cohort']}/{row['sample']}.{row['assembly']}.small_variants.vcf.gz"


class Files:
    """Hash new cache/resources once, and pin previously audited file metadata."""

    def __init__(self):
        self.known = {}

    def get(self, path, expected=None, checked_unix=None):
        """Return a SHA-bound file and register its unchanged-metadata attestation."""
        path = str(path)
        if path not in self.known:
            actual = Path(path)
            before = actual.stat()
            # Reuse only a successful SHA audit predating no file modification.
            if expected and checked_unix and before.st_mtime_ns <= checked_unix * 1e9:
                digest = expected
                mode = "successful_sha_audit_with_unchanged_metadata"
            else:
                digest = sha256sum(actual)
                mode = "fresh_sha256"
            after = actual.stat()
            if (before.st_size, before.st_mtime_ns) != (
                after.st_size,
                after.st_mtime_ns,
            ):
                raise RuntimeError(f"file changed during freezing: {path}")
            self.known[path] = dict(
                path=path,
                sha256=digest,
                bytes=after.st_size,
                mtime_ns=after.st_mtime_ns,
                mode=mode,
            )
        result = self.known[path]
        if expected and result["sha256"] != expected:
            raise RuntimeError(f"frozen file differs: {path}")
        return dict(path=path, sha256=result["sha256"])


def historical_tasks(request, audits, files, proof, source_root):
    """Build every historical condition with its own matched baseline settings."""
    tasks = []
    cache_hashes = {}
    all_inputs = {
        r["sample"]: r["original_task"]
        for r in request["rows"]
        if r["family"] == "full52_wgs_strategies" and r["tool"] == "vep"
    }
    c2_rows = csv.DictReader(
        (
            source_root
            / "benchmarks/figures/output/final/2026-08-25T224623Z-native-annotation-final-v3/source/campaign_metrics.tsv"
        ).open(),
        delimiter="\t",
    )
    c2_cached = {
        (r["sample"], r["tool"], r["level"], r["assay"]): r
        for r in c2_rows
        if r["condition"] == "cached"
    }
    blocksort = files.get(
        source_root
        / "benchmarks/reproduction/final_clinical_campaign/blocksort.frozen.awk"
    )
    for number, row in enumerate(request["rows"], 1):
        print(
            f"Freezing historical task {number}/540: {row['id']} {row['tool']} {row['assay']}",
            flush=True,
        )
        evidence = audits[row["id"]]
        task = dict(
            id=row["id"],
            family=row["family"],
            sample=row["sample"],
            cohort=row["cohort"],
            assembly=row["assembly"],
            tool=row["tool"],
            level=row["level"],
            assay=row["assay"],
            strategy=row["strategy"],
            cpus=8,
            engine="adaptive",
            historical_cached_seconds=row["historical_cached_seconds"],
            environment={},
            annotation_assets=[files.get("/usr/local/bin/bcftools")],
            baseline=dict(
                reuse_decision="approved",
                cpus=8,
                direct_seconds=row["historical_direct_seconds"],
                approval_evidence=proof,
                historical_acceptance=row["historical_acceptance"],
            ),
        )
        if row["family"] == "native_clinical_assays":
            name = f"{row['tool']}_{row['level']}_{row['assembly'].lower()}_af001"
            prefix = (
                "/results/richness/caches"
                if row["tool"] == "fastvep"
                else "/results/snpeff/caches"
            )
            cache = Path(prefix) / name / "cache" / name
            frozen_recipe = (
                source_root
                / "benchmarks/reproduction/final_clinical_campaign/cache_recipes"
                / f"{row['tool']}_{row['level']}_{row['assembly']}.annotation.yaml"
            )
            files.get(cache / "annotation.yaml", sha256sum(frozen_recipe))
            files.get(
                cache / "params.snapshot.yaml",
                sha256sum(
                    frozen_recipe.with_name(
                        f"{row['tool']}_{row['level']}_{row['assembly']}.params.yaml"
                    )
                ),
            )
            params = yaml.safe_load((cache / "params.snapshot.yaml").read_text())
            if row["tool"] == "snpeff":
                params["snpeff_clinical_run"] = (
                    "/results/snpeff/dist/snpeff_clinical_run.sh"
                )
                task["environment"]["SNPEFF_HEAP"] = "8g"
                query = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN"
                if row["level"] == "clinical":
                    query += "\t%INFO/CLNSIG\t%INFO/AF"
                    task["environment"].update(
                        SNPSIFT_JAR="/results/snpeff/dist/SnpSift.jar",
                        CLINVAR_VCF=f"/results/clinical_sources/clinvar_{row['assembly']}.vcf.gz",
                        GNOMAD_VCF=f"/results/clinical_sources/gnomad_exomes_{row['assembly']}.vcf.gz",
                    )
            else:
                query = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\t%INFO/FV_CLINVAR\t%INFO/FV_GNOMAD"
            source = c2_input(row, all_inputs)
            # Use exactly C2's original source, including its explicit duplicate
            # normalization, rather than substituting a differently prepared VCF.
            regions = {
                (
                    "GRCh38",
                    "panel",
                ): "acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed",
                (
                    "GRCh38",
                    "wes",
                ): "twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed",
                (
                    "GRCh37",
                    "panel",
                ): "acmg_sf_v3.3.ensembl115_lifted.mane_select_cds_pad20_GRCh37.bed",
                (
                    "GRCh37",
                    "wes",
                ): "twist_human_core_exome_hg19.chr1-22-X.pad125.merged_GRCh37.bed",
            }
            task["preprocess"] = dict(
                bed=files.get(
                    "/results/richness/regions/"
                    + regions[(row["assembly"], row["assay"])]
                )
                if row["assay"] != "wgs"
                else None
            )
            cached_row = c2_cached[
                (row["sample"], row["tool"], row["level"], row["assay"])
            ]
            task["baseline"].update(
                kind="historical_annotation_digest",
                direct_digest=row["expected_direct_digest"],
                cached_digest=cached_row["digest"],
                query_format=query + "\n",
                blocksort=blocksort,
            )
        elif row["family"] == "synthetic_pipeline_spectrum":
            source = row["original_task"]["input_vcf"]
            cache = Path(row["original_task"]["cache_dir"])
            params = yaml.safe_load(
                Path(row["original_task"]["params_file"]).read_text()
            )
            # SyntheticDelay adds no annotation. Its measured direct timing is
            # retained; compare against the already-retained vanilla recipe
            # reference using the original documented VEP semantic projection.
            reference_row = next(
                r
                for r in request["rows"]
                if r["family"] == "full52_wgs_strategies"
                and r["tool"] == "vep"
                and r["sample"] == row["sample"]
                and r["strategy"]
                == ("gnomad_af_0.01" if "1%" in row["strategy"] else "gnomad_af_0.1")
            )
            reference = audits[reference_row["id"]]
            task["baseline"].update(
                kind="composed_vep_reference",
                direct=dict(
                    path=reference["direct"]["output_path"],
                    sha256=reference["direct"]["output_sha256"],
                ),
                cached=dict(
                    path=reference["cached"]["output_path"],
                    sha256=reference["cached"]["output_sha256"],
                ),
                inherited_direct_difference=True,
                composition="SyntheticDelay emits no annotations; historical paired/composed semantic metadata and unchanged recipe/payload gate timing reuse",
            )
        else:
            command = evidence["cached"]["command"]
            source = evidence["cached"]["input_vcf"]
            cache = Path(argument(command, "-a"))
            param_path = Path(argument(command, "-y"))
            if row["family"] == "vep_matched_assays" and row["assay"] != "wgs":
                original_params = read(
                    Path(
                        "/results/campaigns/external-vep-light-full52-bbef41d7db8e/manifests/strategies.json"
                    )
                )["runtime_params_by_assembly"][row["assembly"]]
                files.get(param_path, original_params["sha256"])
            if row["family"] == "full52_wgs_strategies" or row["assay"] == "wgs":
                campaign = Path(row["cached_root"].split("/measured/")[0])
                strategies = read(campaign / "manifests/strategies.json")
                runtime = strategies["runtime_params_by_assembly"][row["assembly"]]
                files.get(param_path, runtime["sha256"])
                candidates = list(
                    strategies["bundled_strategies_by_assembly"][row["assembly"]]
                ) + [strategies["cohort_strategies"][row["cohort"]]]
                strategy = next(s for s in candidates if s["name"] == row["strategy"])
                if Path(strategy["cache_dir"]) != cache:
                    raise RuntimeError("historical strategy cache directory differs")
                files.get(cache / "annotation.yaml", strategy["annotation_yaml_sha256"])
                provenance = read(Path(strategy["provenance_path"]))
                if any(
                    provenance.get(k) != v
                    for k, v in strategy["provenance_expected"].items()
                ):
                    raise RuntimeError("historical cache provenance differs")
                task["annotation_assets"].append(files.get(strategy["provenance_path"]))
            params = yaml.safe_load(param_path.read_text())
            task["baseline"].update(
                kind="complete_records",
                direct=dict(
                    path=evidence["direct"]["output_path"],
                    sha256=evidence["direct"]["output_sha256"],
                ),
                cached=dict(
                    path=evidence["cached"]["output_path"],
                    sha256=evidence["cached"]["output_sha256"],
                ),
                allow_documented_vep_exception=row["tool"] == "vep",
            )
            task["baseline"]["direct_seconds"] = evidence["direct"]["wall_seconds"]
        if int(params["threads"]) != 8:
            raise RuntimeError(f"historical runtime CPU mismatch: {row['id']}")
        params.pop("temp_dir", None)
        task["parameters"] = params
        if str(cache) not in cache_hashes:
            cache_hashes[str(cache)] = {
                name: files.get(cache / name)["sha256"]
                for name in [
                    "annotation.yaml",
                    "params.snapshot.yaml",
                    "vcfcache_annotated.bcf",
                    "vcfcache_annotated.bcf.csi",
                ]
            }
        task["cache"] = dict(directory=str(cache), files=cache_hashes[str(cache)])
        expected_input = row.get("expected_input_sha256")
        task["input"] = dict(
            **files.get(source, expected_input), records=int(row["input_records"])
        )
        for suffix in [".csi", ".tbi"]:
            if Path(source + suffix).is_file():
                task["annotation_assets"].append(files.get(source + suffix))
        if evidence.get("direct"):
            for role in ["direct", "cached"]:
                reference = evidence[role]
                files.get(
                    reference["output_path"],
                    reference["output_sha256"],
                    evidence.get("audited_unix"),
                )
                task["annotation_assets"].append(
                    files.get(
                        reference["output_path"] + ".csi",
                        reference["output_index_sha256"],
                    )
                )
        task["annotation_asset_group"] = row["tool"].lower()
        # An old WGS copy also appears in the assay-extension plot. A later
        # deduplication step maps it, never times the same condition twice.
        tasks.append(task)
    return tasks


def extension_tasks(tasks, files, proof, output, groups):
    """Add pinned core/clinical vep-rs on the exact twelve GRCh38 assay genomes."""
    templates = [
        t
        for t in tasks
        if t["family"] == "native_clinical_assays"
        and t["tool"] == "fastvep"
        and t["assembly"] == "GRCh38"
    ]
    if len(templates) != 36 or len({t["sample"] for t in templates}) != 12:
        raise ValueError(
            "vep-rs extension must cover twelve genomes and all three assays"
        )
    original = read(
        Path("/results/engine-profile/20261023T190736Z-97545d739/manifest.json")
    )
    for resource in original["resources"]:
        if not resource["path"].startswith("/results/veprs/"):
            continue
        if "sha256" in resource:
            groups["veprs"].append(files.get(resource["path"], resource["sha256"]))
        else:
            for relative, metadata in resource["files"].items():
                groups["veprs"].append(
                    files.get(Path(resource["path"]) / relative, metadata["sha256"])
                )
    for path in [
        "/results/veprs/data/gnomad_joint_gw_GRCh38.vcf.gz",
        "/results/veprs/data/gnomad_joint_gw_GRCh38.vcf.gz.tbi",
        "/results/native-bcf/veprs-caches-20261023b/clinvar.tsv.gz",
        "/results/native-bcf/veprs-caches-20261023b/clinvar.tsv.gz.tbi",
    ]:
        groups["veprs"].append(files.get(path))
    result = []
    for level in ["core", "clinical"]:
        name = f"veprs_{level}"
        origin = Path("/results/native-bcf/veprs-caches-20261023b") / name
        root = output / "cache-overlays" / name
        cache = root / "cache" / name
        cache.mkdir(parents=True)
        (root / "workflow").mkdir()
        (root / "blueprint").mkdir()
        shutil.copy2(origin / "workflow/init.yaml", root / "workflow/init.yaml")
        for filename in [
            "vcfcache_annotated.bcf",
            "vcfcache_annotated.bcf.csi",
            "params.snapshot.yaml",
        ]:
            (cache / filename).symlink_to(origin / "cache" / name / filename)
        recipe = yaml.safe_load(
            (origin / "cache" / name / "annotation.yaml").read_text()
        )
        recipe.pop("tool_version_command", None)
        (cache / "annotation.yaml").write_text(yaml.safe_dump(recipe, sort_keys=False))
        write_json_atomic(
            cache / "STAGING_PROVENANCE.json",
            dict(
                original_annotation=files.get(
                    origin / "cache" / name / "annotation.yaml"
                ),
                repair="remove params-only tool_version_command from staged annotation schema; annotation command and original source cache unchanged",
            ),
        )
        params = yaml.safe_load((cache / "params.snapshot.yaml").read_text())
        params.update(threads=8)
        params.pop("temp_dir", None)
        payload = {
            filename: files.get(cache / filename)["sha256"]
            for filename in [
                "annotation.yaml",
                "params.snapshot.yaml",
                "vcfcache_annotated.bcf",
                "vcfcache_annotated.bcf.csi",
            ]
        }
        for template in templates:
            task = copy.deepcopy(template)
            task.update(
                id=f"veprs-{level}-{template['sample'].replace('_', '-')}-{template['assay']}",
                family="veprs_extension",
                tool="veprs",
                level=level,
                parameters=params,
                environment={},
                annotation_asset_group="veprs",
                cache=dict(directory=str(cache), files=payload),
                historical_cached_seconds=1200 if template["assay"] == "wgs" else 60,
                wall_minutes=360 if template["assay"] == "wgs" else 120,
            )
            task["baseline"] = dict(
                kind="new_paired_reference",
                reuse_decision="new_direct_required",
                cpus=8,
                direct_seconds=None,
                approval_evidence=proof,
                reason="New pinned vep-rs configuration on the matched assay cohort; no compatible historical 8-CPU direct measurement exists",
                known_limit="Clinical CSQ discrepancies are timing-only evidence, never a validated speedup claim",
            )
            result.append(task)
    return result


def first_cells(tasks):
    """Use six actual primary cells to validate launch paths before auto-dispatch."""
    selectors = [
        lambda t: t["tool"] == "vep" and t["assay"] == "panel",
        lambda t: (
            t["family"] == "full52_wgs_strategies"
            and t["tool"] == "fastvep"
            and t["assembly"] == "GRCh37"
        ),
        lambda t: (
            t["family"] == "native_clinical_assays"
            and t["tool"] == "fastvep"
            and t["assay"] == "panel"
        ),
        lambda t: (
            t["tool"] == "snpeff"
            and t["level"] == "clinical"
            and t["assembly"] == "GRCh38"
            and t["assay"] == "panel"
        ),
        lambda t: (
            t["tool"] == "snpeff"
            and t["level"] == "core"
            and t["assembly"] == "GRCh37"
            and t["assay"] == "panel"
        ),
        lambda t: (
            t["tool"] == "veprs" and t["level"] == "core" and t["assay"] == "panel"
        ),
    ]
    result = [next(t["id"] for t in tasks if select(t)) for select in selectors]
    if len(set(result)) != 6:
        raise ValueError(
            "first-stage selectors did not choose six unique primary cells"
        )
    return result


def recipe_assets(tasks, files, groups):
    """Freeze absolute resources resolved by each exact recipe, including SA/plugins."""
    seen = set()
    for task in tasks:
        cache = Path(task["cache"]["directory"])
        command = yaml.safe_load((cache / "annotation.yaml").read_text())[
            "annotation_cmd"
        ]
        for key, value in task["parameters"].items():
            command = command.replace("${params." + key + "}", str(value))
        key = (task["annotation_asset_group"], command)
        if key in seen:
            continue
        seen.add(key)
        paths = set()
        for token in shlex.split(command.replace("\\\n", "")):
            for value in token.split(","):
                path = Path(value)
                if not value.startswith("/") or not path.exists():
                    continue
                if str(path) in {"/", "/results", "/mnt/data", "/mnt/data/tmp", "/tmp"}:
                    continue
                if path.is_dir():
                    paths.update(p for p in path.rglob("*") if p.is_file())
                else:
                    paths.add(path)
                    for suffix in [".fai", ".csi", ".tbi"]:
                        if Path(str(path) + suffix).is_file():
                            paths.add(Path(str(path) + suffix))
        for path in sorted(paths):
            expected = None
            if path.name == "SyntheticDelay.pm":
                expected = read(
                    Path(
                        "/results/campaigns/wgs-pipeline-spectrum-af001-cd87186-kpgp00319/campaign.json"
                    )
                )["source_files"]["SyntheticDelay.pm"]
            groups[task["annotation_asset_group"]].append(files.get(path, expected))
    for group, records in groups.items():
        groups[group] = list({r["path"]: r for r in records}.values())


def deduplicate(tasks):
    """Deduplicate only exact execution/baseline bindings, retaining source IDs."""
    canonical, aliases = {}, {}
    for task in tasks:
        identity = {
            k: task.get(k)
            for k in [
                "input",
                "cache",
                "parameters",
                "environment",
                "preprocess",
                "engine",
                "cpus",
            ]
        }
        identity["baseline"] = {
            k: v
            for k, v in task["baseline"].items()
            if k not in {"historical_acceptance"}
        }
        key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        if key in canonical:
            aliases[task["id"]] = canonical[key]["id"]
            canonical[key].setdefault("source_aliases", []).append(task["id"])
        else:
            canonical[key] = task
    return list(canonical.values()), aliases


def main():
    """Freeze approved historical tasks and asset attestations on shared storage."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    files = Files()
    root = Path("/results/performance-sprint")
    old = root / "reference-audit-20261027/references.jsonl"
    new = root / "reference-audit-v2-20261027/references.jsonl"
    audited = merge_audits(old, new)
    proofs = [files.get(old), files.get(new), files.get(args.request)]
    assets = rows(root / "asset-audit-20261027/assets.jsonl")
    groups = {"vep": [], "fastvep": [], "snpeff": [], "veprs": []}
    for asset in assets:
        if asset["state"] != "match":
            raise RuntimeError("annotation asset audit did not pass")
        record = files.get(
            asset["path"], asset["expected_sha256"], asset["checked_unix"]
        )
        path = asset["path"]
        group = (
            "vep"
            if "ensembl-vep" in path
            else "fastvep"
            if "fastvep" in path
            else "snpeff"
        )
        groups[group].append(record)
    profile = read(
        root.parent / "engine-profile/20261023T190736Z-97545d739/manifest.json"
    )
    for resource in profile["resources"]:
        if "/richness/sa/" not in resource["path"]:
            continue
        if "sha256" in resource:
            groups["fastvep"].append(files.get(resource["path"], resource["sha256"]))
        else:
            for relative, metadata in resource["files"].items():
                groups["fastvep"].append(
                    files.get(Path(resource["path"]) / relative, metadata["sha256"])
                )
    for folder in Path("/results/richness/sa").glob("clinical_GRCh*"):
        for path in sorted(folder.rglob("*")):
            if path.is_file():
                groups["fastvep"].append(files.get(path))
    for filename in [
        "/results/snpeff/dist/jre/bin/java",
        "/results/snpeff/dist/jre/lib/modules",
        "/results/snpeff/dist/snpeff_parallel.sh",
        "/results/snpeff/dist/snpeff_clinical_run.sh",
        "/results/richness/scripts/fastvep_run.sh",
    ]:
        groups["snpeff" if "/snpeff/" in filename else "fastvep"].append(
            files.get(filename)
        )
    tasks = historical_tasks(read(args.request), audited, files, proofs, args.source)
    tasks, aliases = deduplicate(tasks)
    tasks += extension_tasks(tasks, files, proofs, args.output, groups)
    recipe_assets(tasks, files, groups)
    build = read(args.build)
    for task in tasks:
        validate_task(task)
    write_json_atomic(args.output / "tasks.json", tasks)
    freeze = dict(
        approved=True,
        approved_decision="User selected A: stop optimizing and run the full manuscript benchmark",
        tests_passed=True,
        build=build,
        tasks_sha256=sha256sum(args.output / "tasks.json"),
        required_task_ids=[t["id"] for t in tasks],
        node_exclusions={},
        annotation_asset_groups=groups,
        asset_attestations=files.known,
        aliases=aliases,
        first_cells=first_cells(tasks),
        created_unix=time.time(),
        engine_policy="adaptive primary; native remains opt-in; failed experiments disabled",
        metadata_scope="Historical direct timings reused only within exact recipe/input/CPU bindings; digest/composed limitations preserved",
    )
    write_json_atomic(args.output / "freeze.json", freeze)
    prepare(args.output / "campaign", args.output / "tasks.json", freeze)
    write_json_atomic(
        args.output / "ASSEMBLY.json",
        dict(
            source_rows=540,
            tasks=len(tasks),
            aliases=aliases,
            state="ready_not_submitted",
        ),
    )
    if args.submit:
        from benchmarks.performance_sprint.manuscript import start

        receipt = start(args.output / "campaign")
        write_json_atomic(args.output / "SUBMITTED.json", receipt)


if __name__ == "__main__":
    main()
