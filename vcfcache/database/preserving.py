"""Input-preserving annotation with exact biallelic cache matches.

The sample is the immutable backbone. Cache construction may still split
alleles; sample annotation never does. No native routing engine is required.
"""

import hashlib
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HIT = "VCFCACHE_INTERNAL_HIT"
IDENTITY = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%FORMAT\n"


def count(bcftools, path):
    """Fail closed instead of treating an unreadable index as zero misses."""
    value = subprocess.check_output([bcftools, "index", "-n", str(path)], text=True)
    try:
        result = int(value.strip())
    except ValueError as error:
        raise RuntimeError(f"Cannot count indexed records in {path}") from error
    if result < 0:
        raise RuntimeError(f"Invalid record count in {path}")
    return result


def identity(bcftools, path):
    """Hash unchanged columns and all FORMAT data without transferring INFO text.

    Visit contigs in name order so different numeric header IDs cannot create
    false differences. No normalization, allele reordering or field exclusions.
    """
    samples = subprocess.check_output([bcftools, "query", "-l", str(path)])
    index = subprocess.check_output(
        [bcftools, "index", "--stats", str(path)], text=True
    )
    contigs = sorted(line.split("\t")[0] for line in index.splitlines() if line)
    command = [bcftools, "query", "-f", IDENTITY]
    if contigs:
        command += ["-r", ",".join(contigs)]
    command.append(str(path))
    # A file-backed stderr avoids a pipe deadlock on malformed inputs.
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors)
        assert process.stdout is not None
        with process.stdout:
            digest = hashlib.file_digest(process.stdout, "sha256").hexdigest()
        if process.wait():
            errors.seek(0)
            raise RuntimeError(errors.read().decode(errors="replace"))
    return samples, digest


def validate_recipe(bcftools, source, output):
    """Reject truncation and altered keys/genotypes before output is published."""
    expected, actual = count(bcftools, source), count(bcftools, output)
    if actual != expected:
        raise RuntimeError(
            f"Annotation record count changed: {expected} input, {actual} output; "
            "recipes must retain all records, including unannotated records"
        )
    with ThreadPoolExecutor(max_workers=2) as pool:
        before = pool.submit(identity, bcftools, source)
        after = pool.submit(identity, bcftools, output)
        if before.result() != after.result():
            raise RuntimeError(
                "Annotation recipe changed variant identity, record order, "
                "samples or FORMAT data; refusing to publish"
            )
    return actual


def inspect_input(bcftools, source, destination):
    """Check duplicate keys and optionally emit occupied windows in one read.

    bcftools' VCF-to-VCF pairing is one-to-one: duplicate input keys could
    otherwise consume one cached annotation and incorrectly overlay a miss.
    Those uncommon inputs use direct annotation, never silently lose records.
    """
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(
            [bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(source)],
            stdout=subprocess.PIPE,
            stderr=errors,
            text=True,
        )
        current = None
        locus = None
        alleles = set()
        duplicate = False
        assert process.stdout is not None
        with process.stdout, (destination or Path(os.devnull)).open("w") as stream:
            for line in process.stdout:
                contig, position, ref, alt = line.rstrip().split("\t")
                if locus != (contig, position):
                    locus, alleles = (contig, position), set()
                if (ref, alt) in alleles:
                    duplicate = True
                alleles.add((ref, alt))
                start = ((int(position) - 1) // 16384) * 16384
                if current and current[0] == contig and start <= current[2]:
                    current[2] = max(current[2], start + 16384)
                else:
                    if current:
                        stream.write("\t".join(map(str, current)) + "\n")
                    current = [contig, start, start + 16384]
            if current:
                stream.write("\t".join(map(str, current)) + "\n")
        if process.wait():
            errors.seek(0)
            raise RuntimeError(errors.read().decode(errors="replace"))
    return duplicate


def run(manager, cache=None):
    """Run the recipe once, on the full input or the union of cache misses."""
    from vcfcache.database.workflow_manager import BcftoolsCommand

    work = manager.work_dir / ("annotate" if cache else "annotate-nocache")
    work.mkdir(parents=True, exist_ok=True)
    bcftools = str(manager.params_file_content["bcftools_cmd"])
    threads = manager.params_file_content.get("threads", 1)
    if threads < 1:
        raise ValueError("threads must be positive")
    tool = shlex.quote(bcftools)
    source = Path(manager.input_file)
    requested_output = manager._resolve_output_bcf(source.stem)
    if requested_output != "-":
        destination = Path(requested_output)
        if destination.resolve() == source.resolve() or (
            destination.exists() and destination.samefile(source)
        ):
            raise ValueError("Output must not overwrite the input file")
    total = count(bcftools, source)
    tag = manager.nfa_config_content["must_contain_info_tag"]
    mode = os.environ.get("VCFCACHE_PILOT_REGION_MODE", "off")
    if mode not in {"off", "windows"}:
        raise ValueError("VCFCACHE_PILOT_REGION_MODE must be off or windows")

    def command(text):
        return BcftoolsCommand(text, manager.logger, work).run()

    def quoted(path):
        return shlex.quote(str(path))

    restriction = ""
    if cache and total:
        bed = work / "occupied.bed" if mode == "windows" else None
        inspection_start = time.monotonic()
        if inspect_input(bcftools, source, bed):
            manager.logger.warning(
                "Duplicate variant keys: use direct annotation to preserve every occurrence"
            )
            result = run(manager)
            manager.last_run_stats.update(
                missing_variants=total,
                missing_annotated=total,
                cache_matching="bypassed-duplicate-keys",
            )
            return result
        manager.logger.info(
            f"Input key inspection completed in {time.monotonic() - inspection_start:.3f}s"
        )
        if bed:
            restriction = f"-R {quoted(bed)} --regions-overlap 0 "

    recipe_input = source
    missing = total
    if cache:
        for path in (source, cache):
            header = subprocess.check_output(
                [bcftools, "view", "-h", str(path)], text=True
            )
            if f"##INFO=<ID={HIT}," in header:
                raise ValueError(f"Reserved INFO/{HIT} already exists in {path}")
        overlay = work / "cache_overlay.bcf"
        manager.logger.info("Step 1/4: Exact biallelic cache overlay; input unchanged")
        # Exact compares allele SETS, not ALT order. Route every multiallelic
        # record to the recipe rather than reusing allele-indexed CSQ values.
        command(
            f"{tool} annotate --no-version {restriction}--pair-logic exact -i 'N_ALT=1' -k "
            f"-m +{HIT} -a {quoted(cache)} {quoted(source)} -c INFO "
            f"-o {quoted(overlay)} -Ob -W --threads {threads}"
        )
        if count(bcftools, overlay) != total:
            raise RuntimeError("Cache overlay changed input record count")
        recipe_input = work / "misses.bcf"
        manager.logger.info("Step 2/4: Select actual misses, not pre-existing CSQ")
        command(
            f"{tool} view --no-version -i 'INFO/{HIT}!=1 || INFO/{tag}=\"\"' "
            f"{quoted(overlay)} -Ou | {tool} annotate --no-version -x INFO/{HIT} "
            f"-o {quoted(recipe_input)} -Ob -W --threads {threads}"
        )
        missing = count(bcftools, recipe_input)

    # Stage on the destination filesystem so failed validation leaves an old
    # deliverable untouched. The original input is never opened for writing.
    target = manager._resolve_output_bcf(source.stem)
    stdout = target == "-"
    target = work / "stdout.bcf" if stdout else Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".vcfcache-publish-", dir=target.parent
    ) as temporary:
        staged = Path(temporary) / "output.bcf"
        annotated = work / "recipe.bcf" if cache else staged
        result = None
        if not cache or missing:
            manager.logger.info(f"Step 3/4: Annotating {missing} unchanged records")
            auxiliary = manager.output_dir / "auxiliary"
            auxiliary.mkdir(exist_ok=True)
            recipe = manager._substitute_variables(
                manager.nfa_config_content["annotation_cmd"],
                extra_vars={
                    "INPUT_BCF": str(recipe_input),
                    "OUTPUT_BCF": str(annotated),
                    "AUXILIARY_DIR": str(auxiliary),
                },
            )
            result = command(recipe)
            for filename, output in (
                ("stdout.txt", "annotation_tool.log"),
                ("stderr.txt", "annotation_tool_err.log"),
            ):
                path = work / filename
                if path.exists() and path.stat().st_size:
                    shutil.copyfile(path, manager.output_dir / output)
            if not any(auxiliary.iterdir()):
                auxiliary.rmdir()
            try:
                validation_start = time.monotonic()
                manager._validate_info_tag(annotated, tag)
                validate_recipe(bcftools, recipe_input, annotated)
                manager.logger.info(
                    f"Recipe identity/FORMAT validation completed in {time.monotonic() - validation_start:.3f}s"
                )
            except Exception:
                if annotated.exists():
                    shutil.copyfile(annotated, work / "failed_recipe.bcf")
                    if Path(str(annotated) + ".csi").exists():
                        shutil.copyfile(
                            Path(str(annotated) + ".csi"),
                            work / "failed_recipe.bcf.csi",
                        )
                raise
        if cache:
            manager.logger.info("Step 4/4: Overlay annotations; retain every record")
            if missing:
                merged = work / "merged.bcf"
                command(
                    f"{tool} annotate --no-version {restriction}--pair-logic exact "
                    f"-a {quoted(annotated)} {quoted(overlay)} -c INFO "
                    f"-o {quoted(merged)} -Ob -W --threads {threads}"
                )
            else:
                merged = overlay
            result = command(
                f"{tool} annotate --no-version -x INFO/{HIT} {quoted(merged)} "
                f"-o {quoted(staged)} -Ob -W --threads {threads}"
            )
        if count(bcftools, staged) != total:
            raise RuntimeError("Final output record count differs from input")
        if stdout:
            with staged.open("rb") as stream:
                shutil.copyfileobj(stream, sys.stdout.buffer)
        else:
            os.replace(staged, target)
            os.replace(Path(str(staged) + ".csi"), Path(str(target) + ".csi"))
    manager.last_run_stats = {
        "input_variants": total,
        "output_variants": total,
        "missing_variants": missing if cache else None,
        "missing_annotated": missing if cache else None,
        "dropped_variants": 0,
        "input_preprocessing": "none",
        "cache_matching": "biallelic-exact" if cache else None,
        "region_mode": mode if cache else "off",
    }
    return result
