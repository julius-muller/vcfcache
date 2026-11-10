"""Input-preserving annotation with exact biallelic cache matches.

The sample is the immutable backbone. Cache construction may still split
alleles; sample annotation never does. The cache answers only the records it
can answer exactly; every other record reaches the recipe unchanged, as direct
annotation would see it. An unusual input therefore costs hits, not correctness.

The cached path runs coordinate-disjoint shards concurrently and publishes
their outputs by BGZF block copy, so no record is decoded or recompressed after
it has been annotated. Cache access is restricted to windows around the input
when the cache is much larger than the input, and to whole contigs otherwise.
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
KEY = "%CHROM\t%POS\t%REF\t%ALT\n"
WINDOW = 16384
# Window lookups pay per-window index seeks; they win once the cache holds far
# more records than the input (measured: wins at 250x, loses at 4x).
WINDOW_RATIO = 20
# Upper bound on shards; contigs are grouped into shards of comparable size.
SHARDS_PER_THREAD = 4
UNBOUNDED = 2**31 - 1
CHUNK = 1 << 24
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


def index_stats(bcftools, path):
    """Return (contig, length, records) of occupied contigs from the index alone."""
    listing = subprocess.check_output(
        [bcftools, "index", "--stats", str(path)], text=True
    )
    occupied = []
    for line in listing.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3 and fields[2].isdigit() and int(fields[2]):
            length = int(fields[1]) if fields[1].isdigit() else UNBOUNDED
            occupied.append((fields[0], length, int(fields[2])))
    return occupied


def plan_shards(occupied, limit):
    """Group contigs, in index order, into shards of comparable record count.

    Every shard is a contiguous run of contigs, so the shard outputs
    concatenate into a correctly sorted file without decoding any record.
    """
    if not occupied:
        return []
    number = min(max(1, limit), len(occupied))
    target = sum(records for _, _, records in occupied) / number
    shards, current, load = [], [], 0
    for index, (contig, length, records) in enumerate(occupied):
        current.append((contig, length))
        load += records
        if load >= target and index < len(occupied) - 1 and len(shards) + 1 < number:
            shards.append(current)
            current, load = [], 0
    if current:
        shards.append(current)
    return shards


def inspect_keys(bcftools, source, bed, windowed):
    """Find duplicate variant keys and, optionally, the occupied windows.

    bcftools pairs records one-to-one: of two identical input keys only the
    first would take a cached annotation, and the recipe result of the second
    would later be paired with the first. A shard holding such keys bypasses
    the cache. Keys compare without case, as bcftools compares alleles.

    Returns:
        (duplicate, rows) where rows are merged BED windows, or None unless
        windowed and duplicate-free.
    """
    blocks: dict = {}
    order = []
    duplicate = False
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(
            [bcftools, "query", "-f", KEY, "-R", str(bed), "--regions-overlap", "0"]
            + [str(source)],
            stdout=subprocess.PIPE,
            stderr=errors,
        )
        assert process.stdout is not None
        tail = b""
        with process.stdout:
            while chunk := process.stdout.read(CHUNK):
                data = tail + chunk
                cut = data.rfind(b"\n")
                if cut < 0:
                    tail = data
                    continue
                lines = data[:cut].split(b"\n")
                if len(set(data[:cut].upper().split(b"\n"))) != len(lines):
                    duplicate = True
                    break
                if windowed:
                    for line in lines:
                        contig, position = line.split(b"\t", 2)[:2]
                        if contig not in blocks:
                            blocks[contig] = set()
                            order.append(contig)
                        blocks[contig].add((int(position) - 1) // WINDOW)
                # Keys of the last locus may continue in the next chunk.
                locus = b"\t".join(lines[-1].split(b"\t", 2)[:2]) + b"\t"
                keep = len(lines)
                while keep and lines[keep - 1].startswith(locus):
                    keep -= 1
                tail = b"".join(line + b"\n" for line in lines[keep:]) + data[cut + 1 :]
        if duplicate:
            process.kill()
            process.wait()
            return True, None
        if process.wait():
            errors.seek(0)
            raise RuntimeError(errors.read().decode(errors="replace"))
    if not windowed:
        return False, None
    rows = []
    for contig in order:
        previous = None
        for block in sorted(blocks[contig]):
            start, end = block * WINDOW, (block + 1) * WINDOW
            if previous and start == previous[1]:
                previous[1] = end
                continue
            if previous:
                rows.append(f"{contig.decode()}\t{previous[0]}\t{previous[1]}")
            previous = [start, end]
        if previous:
            rows.append(f"{contig.decode()}\t{previous[0]}\t{previous[1]}")
    return False, rows


def run_recipe(manager, work, command, recipe_input, annotated):
    """Run the annotation recipe once; its output is checked by check_recipe."""
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
    return result


def check_recipe(manager, work, recipe_input, annotated, tag, bcftools):
    """Refuse recipe output that lacks the tag or altered any record."""
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
                    Path(str(annotated) + ".csi"), work / "failed_recipe.bcf.csi"
                )
        raise


def run(manager, cache=None):
    """Run the recipe once, on the full input or the union of cache misses."""
    from vcfcache.database.workflow_manager import BcftoolsCommand

    work = manager.work_dir / ("annotate" if cache else "annotate-nocache")
    work.mkdir(parents=True, exist_ok=True)
    bcftools = str(manager.params_file_content["bcftools_cmd"])
    threads = manager.params_file_content.get("threads", 1)
    if threads < 1:
        raise ValueError("threads must be positive")
    level = getattr(manager, "output_level", 1)
    if level not in range(10):
        raise ValueError("Output compression level must be between 0 and 9")
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

    def command(text, folder=work):
        return BcftoolsCommand(text, manager.logger, folder).run()

    def quoted(path):
        return shlex.quote(str(path))

    # Stage on the destination filesystem so failed validation leaves an old
    # deliverable untouched. The original input is never opened for writing.
    stdout = requested_output == "-"
    target = work / "stdout.bcf" if stdout else Path(requested_output)
    target.parent.mkdir(parents=True, exist_ok=True)
    stats = {}
    missing = total
    with tempfile.TemporaryDirectory(
        prefix=".vcfcache-publish-", dir=target.parent
    ) as temporary:
        staged = Path(temporary) / "output.bcf"
        if not cache:
            manager.logger.info(f"Annotating {total} unchanged records")
            result = run_recipe(manager, work, command, source, staged)
            check_recipe(manager, work, source, staged, tag, bcftools)
        else:
            for path in (source, cache):
                header = subprocess.check_output(
                    [bcftools, "view", "-h", str(path)], text=True
                )
                if f"##INFO=<ID={HIT}," in header:
                    raise ValueError(f"Reserved INFO/{HIT} already exists in {path}")
            # A header-only cache gives bypassed shards the same header as
            # every other shard without reading any cache record.
            blank = work / "blank_cache.bcf"
            command(
                f"{tool} view --no-version -h {quoted(cache)} | "
                f"{tool} view --no-version -Ob -W -o {quoted(blank)}"
            )
            shards = plan_shards(
                index_stats(bcftools, source), SHARDS_PER_THREAD * threads
            )
            windowed = bool(total) and count(bcftools, cache) > WINDOW_RATIO * total
            folders = []
            for number in range(max(1, len(shards))):
                folder = work / f"shard-{number}"
                folder.mkdir(exist_ok=True)
                folders.append(folder)
            regions = [""] * len(folders)
            bypassed = [False] * len(folders)
            selection = f"-i 'INFO/{HIT}!=1 || INFO/{tag}=\"\"'"

            def lookup(number):
                """Overlay exact cache hits on one shard and select its misses."""
                folder = folders[number]
                lookup_cache = cache
                if shards:
                    bed = folder / "contigs.bed"
                    bed.write_text(
                        "".join(f"{c}\t0\t{length}\n" for c, length in shards[number])
                    )
                    bypassed[number], rows = inspect_keys(
                        bcftools, source, bed, windowed
                    )
                    if rows:
                        bed = folder / "windows.bed"
                        bed.write_text("".join(f"{row}\n" for row in rows))
                    regions[number] = f"-R {quoted(bed)} --regions-overlap 0 "
                    if bypassed[number]:
                        lookup_cache = blank
                # Exact compares allele SETS, not ALT order. Route every
                # multiallelic record to the recipe rather than reusing
                # allele-indexed annotation values.
                command(
                    f"{tool} annotate --no-version {regions[number]}--pair-logic exact "
                    f"-i 'N_ALT=1' -k -m +{HIT} -a {quoted(lookup_cache)} "
                    f"{quoted(source)} -c INFO -o overlay.bcf -Ob1 -W --threads 1",
                    folder,
                )
                command(
                    f"{tool} view --no-version {selection} overlay.bcf -Ou | "
                    f"{tool} annotate --no-version -x INFO/{HIT} -o misses.bcf -Ob",
                    folder,
                )
                return count(bcftools, folder / "overlay.bcf")

            manager.logger.info(
                f"Step 1/4: Exact biallelic cache overlay in {len(folders)} shard(s), "
                f"cache access by {'windows' if windowed else 'contigs'}; input unchanged"
            )
            with ThreadPoolExecutor(max_workers=threads) as pool:
                if sum(pool.map(lookup, range(len(folders)))) != total:
                    raise RuntimeError("Cache overlay changed input record count")
            if any(bypassed):
                manager.logger.warning(
                    f"Duplicate variant keys: {sum(bypassed)} shard(s) bypass the "
                    "cache so that every occurrence is annotated by the recipe"
                )
            manager.logger.info("Step 2/4: Select actual misses, not pre-existing tags")
            recipe_input = work / "misses.bcf"
            parts = " ".join(quoted(folder / "misses.bcf") for folder in folders)
            command(
                f"{tool} concat --no-version --naive {parts} -o {quoted(recipe_input)} && "
                f"{tool} index -f {quoted(recipe_input)}"
            )
            missing = count(bcftools, recipe_input)
            annotated = work / "recipe.bcf"
            # The recipe output is checked while the shards are merged; nothing
            # is published before that check has passed.
            checker = ThreadPoolExecutor(max_workers=1)
            checked = None
            if missing:
                manager.logger.info(f"Step 3/4: Annotating {missing} unchanged records")
                run_recipe(manager, work, command, recipe_input, annotated)
                checked = checker.submit(
                    check_recipe, manager, work, recipe_input, annotated, tag, bcftools
                )
            manager.logger.info("Step 4/4: Overlay annotations; retain every record")

            def merge(number):
                """Add recipe annotations to one shard and drop the marker."""
                overlay = (
                    f"{regions[number]}--pair-logic exact -a {quoted(annotated)} -c INFO "
                    if missing
                    else ""
                )
                command(
                    f"{tool} annotate --no-version {overlay}-x INFO/{HIT} overlay.bcf "
                    f"-o final.bcf -Ob{level} --threads 1",
                    folders[number],
                )

            try:
                with ThreadPoolExecutor(max_workers=threads) as pool:
                    list(pool.map(merge, range(len(folders))))
            finally:
                # A recipe that altered records explains a failed merge best.
                checker.shutdown(wait=True)
                if checked:
                    checked.result()
            # Gather by BGZF block copy: shards are coordinate-disjoint and in
            # index order, so nothing is decoded or recompressed here.
            parts = " ".join(quoted(folder / "final.bcf") for folder in folders)
            result = command(
                f"{tool} concat --no-version --naive {parts} -o {quoted(staged)} && "
                f"{tool} index -f --threads {threads} {quoted(staged)}"
            )
            stats = {
                "shards": len(folders),
                "bypassed_shards": sum(bypassed),
                "cache_access": "windows" if windowed else "contigs",
                "output_level": level,
            }
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
        **stats,
    }
    return result
