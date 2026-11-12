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
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HIT = "VCFCACHE_INTERNAL_HIT"
KEY = "%CHROM\t%POS\t%REF\t%ALT\n"
WINDOW = 16384
# Window lookups pay per-window index seeks; they win once the cache holds far
# more records than the input (measured: wins at 250x, loses at 4x).
WINDOW_RATIO = 20
# The input is cut into SHARDS_PER_THREAD * threads shards of equal record
# count, regardless of contig boundaries; more shards than cores lets the pool
# absorb uneven cache density. Small inputs use fewer, larger shards.
SHARDS_PER_THREAD = 2
MIN_SHARD = 5000
# Every SAMPLE-th record's position is kept to place the cuts.
SAMPLE = 1024
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


def scan_contig(bcftools, source, contig, windowed):
    """Read one contig's variant keys once, at C speed.

    bcftools pairs records one-to-one: of two identical input keys only the
    first would take a cached annotation, and the recipe result of the second
    would later be paired with the first. Every record at such a position is
    therefore annotated by the recipe; all other records still use the cache.
    Keys compare without case, as bcftools compares alleles.

    Returns:
        dict with the record count, (index, POS) of every SAMPLE-th record,
        positions holding duplicate keys, and occupied windows when windowed.
    """
    samples = []
    duplicates: dict = {}
    blocks: set = set()
    seen = 0
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(
            [bcftools, "query", "-f", KEY, "-r", contig, str(source)],
            stdout=subprocess.PIPE,
            stderr=errors,
        )
        assert process.stdout is not None
        tail = b""

        def scan(lines):
            """Count, sample, and record duplicates and windows of whole loci."""
            nonlocal seen
            for j in range((-seen) % SAMPLE, len(lines), SAMPLE):
                samples.append((seen + j, int(lines[j].split(b"\t", 2)[1])))
            seen += len(lines)
            keys = [line.upper() for line in lines]
            if len(set(keys)) != len(keys):
                counts = Counter(keys)
                for key, line in zip(keys, lines, strict=True):
                    if counts[key] > 1:
                        duplicates[int(line.split(b"\t", 2)[1])] = None
            if windowed:
                for line in lines:
                    blocks.add((int(line.split(b"\t", 2)[1]) - 1) // WINDOW)

        with process.stdout:
            while chunk := process.stdout.read(CHUNK):
                data = tail + chunk
                cut = data.rfind(b"\n")
                if cut < 0:
                    tail = data
                    continue
                lines = data[:cut].split(b"\n")
                # Keys of the last locus may continue in the next chunk.
                locus = b"\t".join(lines[-1].split(b"\t", 2)[:2]) + b"\t"
                keep = len(lines)
                while keep and lines[keep - 1].startswith(locus):
                    keep -= 1
                scan(lines[:keep])
                tail = b"".join(line + b"\n" for line in lines[keep:]) + data[cut + 1 :]
        if tail.strip():
            scan(tail.rstrip(b"\n").split(b"\n"))
        if process.wait():
            errors.seek(0)
            raise RuntimeError(errors.read().decode(errors="replace"))
    # Insertion order is file order, so the positions are already sorted.
    return dict(
        records=seen,
        samples=samples,
        duplicates=list(duplicates),
        blocks=blocks if windowed else None,
    )


def plan_ranges(order, scans, parts):
    """Cut the input into about equal record counts, between positions.

    A cut falls before a sampled record's POS, so records sharing a position
    always stay in one shard. Shards may span contigs; bcftools reads regions
    in the order they are listed, which here is file order.

    Returns:
        Shards, each a list of (contig, first POS, last POS) in file order.
    """
    total = sum(scans[c]["records"] for c in order)
    if not total:
        return []
    parts = max(1, min(parts, total))
    cuts = set()
    offset = 0
    targets = [round(k * total / parts) for k in range(1, parts)]
    for contig in order:
        records = scans[contig]["records"]
        for target in targets:
            if offset <= target < offset + records:
                local = target - offset
                chosen = None
                for index, position in scans[contig]["samples"]:
                    if index > local:
                        break
                    chosen = (index, position)
                if chosen is None or chosen[0] == 0:
                    cuts.add((contig, 1))
                else:
                    cuts.add((contig, chosen[1]))
        offset += records
    shards: list = [[]]
    for contig in order:
        starts = sorted(p for c, p in cuts if c == contig)
        if starts and starts[0] == 1 and shards[-1]:
            shards.append([])
        bounds = [p for p in starts if p > 1]
        first = 1
        for position in bounds:
            shards[-1].append((contig, first, position - 1))
            shards.append([])
            first = position
        shards[-1].append((contig, first, UNBOUNDED))
    return [shard for shard in shards if shard]


def shard_rows(shard, scans):
    """BED rows of one shard: whole ranges, or its windows clipped to them."""
    rows = []
    for contig, first, last in shard:
        blocks = scans[contig]["blocks"]
        if blocks is None:
            rows.append(f"{contig}\t{first - 1}\t{last}")
            continue
        previous = None
        for block in sorted(blocks):
            start = max(block * WINDOW, first - 1)
            end = min((block + 1) * WINDOW, last)
            if start >= end:
                continue
            if previous and start <= previous[1]:
                previous[1] = max(previous[1], end)
                continue
            if previous:
                rows.append(f"{contig}\t{previous[0]}\t{previous[1]}")
            previous = [start, end]
        if previous:
            rows.append(f"{contig}\t{previous[0]}\t{previous[1]}")
    return rows


def shard_duplicates(shard, scans):
    """Return CHROM-tab-POS lines of the duplicate positions in one shard."""
    return [
        f"{contig}\t{position}"
        for contig, first, last in shard
        for position in scans[contig]["duplicates"]
        if first <= position <= last
    ]


def contig_order(bcftools, source):
    """Occupied contigs in the order the output must follow."""
    return [contig for contig, _, _ in index_stats(bcftools, source)]


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
            headers = [
                subprocess.check_output([bcftools, "view", "-h", str(path)], text=True)
                for path in (source, cache)
            ]
            for path, header in zip((source, cache), headers, strict=True):
                if f"##INFO=<ID={HIT}," in header:
                    raise ValueError(f"Reserved INFO/{HIT} already exists in {path}")
            # An input that already declares the tag may carry annotation the
            # cache did not supply; only then must hits be marked explicitly.
            # Otherwise an empty tag identifies a miss, which saves the marker
            # and its removal on every record.
            marked = f"##INFO=<ID={tag}," in headers[0]
            windowed = bool(total) and count(bcftools, cache) > WINDOW_RATIO * total
            # One parallel pass over the input's keys places equal-record cuts
            # and finds duplicate positions and, for small inputs, windows.
            order = contig_order(bcftools, source)
            scan_start = time.monotonic()
            with ThreadPoolExecutor(max_workers=threads) as pool:
                scans = dict(
                    zip(
                        order,
                        pool.map(
                            lambda c: scan_contig(bcftools, source, c, windowed),
                            order,
                        ),
                        strict=True,
                    )
                )
            if sum(s["records"] for s in scans.values()) != total:
                raise RuntimeError("Key scan read a different number of records")
            parts = min(SHARDS_PER_THREAD * threads, -(-total // MIN_SHARD))
            shards = plan_ranges(order, scans, parts)
            manager.logger.info(
                f"Key scan of {len(order)} contig(s) and planning of {len(shards)} "
                f"shard(s) completed in {time.monotonic() - scan_start:.3f}s"
            )
            folders = []
            for number in range(max(1, len(shards))):
                folder = work / f"shard-{number}"
                folder.mkdir(exist_ok=True)
                folders.append(folder)
            regions = [""] * len(folders)
            # Positions with duplicate keys, per shard; see scan_contig.
            duplicates: list = [shard_duplicates(shard, scans) for shard in shards] or [
                []
            ]
            mark = f"-m +{HIT} " if marked else ""
            strip = f"-x INFO/{HIT} " if marked else ""

            def lookup(number):
                """Overlay exact cache hits on one shard and select its misses."""
                folder = folders[number]
                if shards:
                    bed = folder / "regions.bed"
                    rows = shard_rows(shards[number], scans)
                    bed.write_text("".join(f"{row}\n" for row in rows))
                    regions[number] = f"-R {quoted(bed)} --regions-overlap 0 "
                # Exact compares allele SETS, not ALT order. Route every
                # multiallelic record to the recipe rather than reusing
                # allele-indexed annotation values.
                command(
                    f"{tool} annotate --no-version {regions[number]}--pair-logic exact "
                    f"-i 'N_ALT=1' -k {mark}-a {quoted(cache)} "
                    f"{quoted(source)} -c INFO -o overlay.bcf -Ob1 -W --threads 1",
                    folder,
                )
                # Without duplicates the misses are written directly; with
                # them, the misses elsewhere and every record at a duplicate
                # position (taken from the input, never from the overlay) are
                # merged by position.
                dups = folder / "duplicates.tsv"
                if duplicates[number]:
                    dups.write_text("".join(f"{row}\n" for row in duplicates[number]))
                elsewhere = f"-T ^{quoted(dups)} " if duplicates[number] else ""
                target = "elsewhere.bcf -W" if duplicates[number] else "misses.bcf"
                if marked:
                    command(
                        f"{tool} view --no-version {elsewhere}"
                        f"-i 'INFO/{HIT}!=1 || INFO/{tag}=\"\"' overlay.bcf -Ou | "
                        f"{tool} annotate --no-version {strip}-o {target} -Ob",
                        folder,
                    )
                else:
                    command(
                        f"{tool} view --no-version {elsewhere}-i 'INFO/{tag}=\"\"' "
                        f"overlay.bcf -o {target} -Ob --threads 1",
                        folder,
                    )
                if duplicates[number]:
                    command(
                        f"{tool} view --no-version -R {quoted(dups)} --regions-overlap 0 "
                        f"{quoted(source)} -o duplicated.bcf -Ob -W --threads 1 && "
                        f"{tool} concat --no-version -a elsewhere.bcf duplicated.bcf "
                        "-o misses.bcf -Ob --threads 1",
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
            sites = sum(len(d) for d in duplicates)
            if sites:
                manager.logger.warning(
                    f"Duplicate variant keys at {sites} position(s): every record "
                    "there is annotated by the recipe"
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
            # At a duplicate position the overlay may hold a cached value that
            # the merge only replaces where the recipe wrote the tag. Should
            # the recipe leave such a record without it, merge onto a base
            # that holds the input's own records there instead.
            clean = False
            if sites and missing:
                every = work / "duplicates.tsv"
                every.write_text("".join(f"{row}\n" for d in duplicates for row in d))
                if checked:
                    checked.result()
                untagged = subprocess.run(
                    [bcftools, "view", "-H", "-R", str(every), "--regions-overlap", "0"]
                    + ["-i", f'INFO/{tag}=""', str(annotated)],
                    capture_output=True,
                    check=True,
                ).stdout
                clean = bool(untagged.strip())
                if clean:
                    manager.logger.info(
                        "Recipe left records at duplicate positions without "
                        f"INFO/{tag}; merging onto the input's own records there"
                    )
            manager.logger.info("Step 4/4: Overlay annotations; retain every record")

            def merge(number):
                """Add recipe annotations to one shard and drop any marker."""
                overlay = (
                    f"{regions[number]}--pair-logic exact -a {quoted(annotated)} -c INFO "
                    if missing
                    else ""
                )
                base = "overlay.bcf"
                if clean and duplicates[number]:
                    base = "base.bcf"
                    command(
                        f"{tool} view --no-version -T ^duplicates.tsv overlay.bcf "
                        "-o kept.bcf -Ob1 -W --threads 1 && "
                        f"{tool} concat --no-version -a kept.bcf duplicated.bcf "
                        "-o base.bcf -Ob1 -W --threads 1",
                        folders[number],
                    )
                if overlay or strip:
                    command(
                        f"{tool} annotate --no-version {overlay}{strip}{base} "
                        f"-o final.bcf -Ob{level} --threads 1",
                        folders[number],
                    )
                else:
                    command(
                        f"{tool} view --no-version overlay.bcf -o final.bcf "
                        f"-Ob{level} --threads 1",
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
                "duplicate_sites": sites,
                "duplicate_clean_base": clean,
                "cache_access": "windows" if windowed else "contigs",
                "output_level": level,
                "hit_marker": marked,
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
