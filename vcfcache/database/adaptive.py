# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Region-restricted, contig-parallel cache annotation using only bcftools.

The legacy path reads the whole cache on every run, because `bcftools annotate
-a` streams both files in lockstep unless regions are given. It also writes the
fully annotated input as a compressed intermediate and then reads it twice.

This engine keeps the same logical steps and the same bcftools commands, but
restricts cache access to the regions the input occupies, runs
coordinate-disjoint shards concurrently, and keeps the full-size intermediates
uncompressed so the payload is BGZF-compressed once rather than three times.

Region restriction only pays while the merged window list stays short: past a
few thousand windows the per-window index seeks cost more than a plain scan, so
larger inputs use whole-contig regions and keep the sharding and compression
benefits only.
"""

import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Optional, Tuple

from vcfcache.database.workflow_manager import BcftoolsCommand

WINDOW = 16384
UNBOUNDED = 2**31 - 1


def index_stats(bcftools: str, path: Path) -> List[Tuple[str, int, int]]:
    """Return (contig, length, records) for occupied contigs, from the index alone.

    Args:
        bcftools: bcftools executable from params.yaml
        path: Indexed BCF to inspect

    Returns:
        List of (contig, length, record count), omitting contigs without records.
        Contigs of undeclared length report the largest position bcftools accepts.
    """
    result = subprocess.run(
        [bcftools, "index", "--stats", str(path)],
        capture_output=True,
        text=True,
        check=True,
    )
    occupied = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3 and fields[2].isdigit() and int(fields[2]):
            length = int(fields[1]) if fields[1].isdigit() else UNBOUNDED
            occupied.append((fields[0], length, int(fields[2])))
    return occupied


def plan_shards(
    occupied: List[Tuple[str, int, int]], workers: int
) -> List[List[Tuple[str, int]]]:
    """Split contigs into consecutive shards of comparable record count.

    The shards stay in index order and each holds a contiguous run of contigs,
    so the gathered shard outputs concatenate into a correctly sorted file
    without decoding or re-sorting any record.

    Args:
        occupied: (contig, length, records) rows from index_stats, in index order
        workers: Maximum number of concurrent shards

    Returns:
        Lists of (contig, length) in index order; empty when nothing is occupied
    """
    if not occupied:
        return []
    count = min(max(1, workers), len(occupied))
    target = sum(records for _, _, records in occupied) / count
    shards: List[List[Tuple[str, int]]] = []
    current: List[Tuple[str, int]] = []
    load = 0
    for index, (contig, length, records) in enumerate(occupied):
        current.append((contig, length))
        load += records
        # Close once this shard carries its share, provided a contig is left
        # for the final shard. Uneven contigs make fewer shards than workers;
        # insisting on the full count instead would merge everything into one.
        if load >= target and index < len(occupied) - 1 and len(shards) + 1 < count:
            shards.append(current)
            current, load = [], 0
    if current:
        shards.append(current)
    return shards


def write_regions(
    bcftools: str,
    input_bcf: Path,
    shards: List[List[Tuple[str, int]]],
    root: Path,
    windowed: bool,
) -> List[Path]:
    """Write one regions file per shard, as merged windows or whole contigs.

    Args:
        bcftools: bcftools executable from params.yaml
        input_bcf: Normalized, indexed input whose positions define the windows
        shards: (contig, length) groups from plan_shards
        root: Directory to write the regions files into
        windowed: Emit merged position windows instead of whole-contig ranges

    Returns:
        Regions file paths in shard order
    """
    blocks: dict = {}
    if windowed:
        # Only reached for small inputs, where one extra pass is cheap.
        listing = subprocess.run(
            [bcftools, "query", "-f", "%CHROM\t%POS\n", str(input_bcf)],
            capture_output=True,
            text=True,
            check=True,
        )
        for line in listing.stdout.splitlines():
            contig, position = line.split("\t")
            blocks.setdefault(contig, set()).add((int(position) - 1) // WINDOW)
    paths = []
    for number, contigs in enumerate(shards):
        rows = []
        for contig, length in contigs:
            if windowed:
                for block in sorted(blocks.get(contig, ())):
                    rows.append(f"{contig}\t{block * WINDOW}\t{(block + 1) * WINDOW}")
            else:
                # Three BED columns throughout: bcftools cannot parse a bare
                # contig name from a regions file.
                rows.append(f"{contig}\t0\t{length}")
        path = root / f"regions-{number}.bed"
        path.write_text("".join(f"{row}\n" for row in rows))
        paths.append(path)
    return paths


def run_adaptive(
    manager,
    db_bcf: Path,
    preserve_unannotated: bool = False,
    skip_split_multiallelic: bool = False,
    workers: int = 0,
    max_windows: int = 5000,
) -> subprocess.CompletedProcess:
    """Annotate using the cache with region-restricted, sharded bcftools steps.

    Args:
        manager: WorkflowManager providing configuration, paths and logging
        db_bcf: Path to cache BCF (vcfcache_annotated.bcf)
        preserve_unannotated: Keep variants without annotation in the output
        skip_split_multiallelic: Skip splitting multiallelic variants
        workers: Concurrent shards; 0 uses the params.yaml thread budget
        max_windows: Largest merged window list worth restricting the cache to

    Returns:
        subprocess.CompletedProcess of the publishing step
    """
    manager.logger.info("Annotating sample using cache (adaptive engine)")

    work_task = manager.work_dir / "annotate"
    work_task.mkdir(parents=True, exist_ok=True)

    sample_name = manager.input_file.stem
    input_bcf = manager.input_file
    bcftools = manager.params_file_content["bcftools_cmd"]
    tag = manager.nfa_config_content["must_contain_info_tag"]
    threads = manager.params_file_content.get("threads", 1)

    # Step 0: identical to the legacy path, and the index it writes is what
    # lets every later step address shards without reading the records again.
    if skip_split_multiallelic:
        manager.logger.info("Step 0/4: Removing spanning deletion alleles (ALT=*)")
        filtered_input = work_task / f"{sample_name}_no_span_del.bcf"
        filter_cmd = f"{bcftools} view -e 'ALT=\"*\"' {input_bcf} -o {filtered_input} -Ob -W --threads {threads}"
        BcftoolsCommand(filter_cmd, manager.logger, work_task).run()
        input_bcf = filtered_input
    else:
        manager.logger.info(
            "Step 0/4: Splitting multiallelic variants and removing spanning deletions"
        )
        normalized_input = work_task / f"{sample_name}_normalized.bcf"
        norm_cmd = (
            f"{bcftools} norm -m- {input_bcf} -Ob | "
            f"{bcftools} view -e 'ALT=\"*\"' -o {normalized_input} -Ob -W --threads {threads}"
        )
        BcftoolsCommand(norm_cmd, manager.logger, work_task).run()
        input_bcf = normalized_input

    input_count = manager._count_variants(Path(input_bcf), bcftools)

    occupied = index_stats(bcftools, Path(input_bcf))
    shards = plan_shards(occupied, workers or threads or 1)
    windowed = bool(shards) and bool(input_count) and input_count <= max_windows
    # An empty input has no occupied contig to shard; one unrestricted pass
    # keeps the remaining steps and the recorded statistics well defined.
    regions: List[Optional[Path]] = list(
        write_regions(bcftools, Path(input_bcf), shards, work_task, windowed)
    ) or [None]
    shard_threads = max(1, threads // len(regions))
    manager.logger.info(
        f"Adaptive engine: {len(regions)} shard(s); "
        f"cache access={'merged windows' if windowed else 'whole contigs'}"
    )

    # Each concurrent step needs its own directory, because BcftoolsCommand
    # writes command.sh and the captured streams into the one it is given.
    dirs = []
    for number in range(len(regions)):
        path = work_task / f"shard-{number}"
        path.mkdir(parents=True, exist_ok=True)
        dirs.append(path)

    def restrict(number: int) -> str:
        """Region arguments for a shard, empty when nothing is restricted."""
        return f"-R {regions[number]} --regions-overlap 0 " if regions[number] else ""

    def run_shards(label: str, build) -> None:
        """Run one per-shard command concurrently and fail on the first error."""
        with ThreadPoolExecutor(max_workers=len(regions)) as pool:
            list(
                pool.map(
                    lambda number: BcftoolsCommand(
                        build(number), manager.logger, dirs[number]
                    ).run(),
                    range(len(regions)),
                )
            )
        manager.logger.debug(f"Completed {label} across {len(regions)} shard(s)")

    # Step 1: cache annotations, restricted to the occupied regions and written
    # uncompressed, because this full-size intermediate is only read again here.
    manager.logger.info("Step 1/4: Adding cache annotations")
    run_shards(
        "cache annotation",
        lambda number: (
            f"{bcftools} annotate {restrict(number)}-a {db_bcf} {input_bcf} -c INFO "
            f"-o {dirs[number]}/step1.bcf -Ob0 -W --threads {shard_threads}"
        ),
    )

    # Step 2: misses per shard, then a block copy into the single recipe input.
    manager.logger.info("Step 2/4: Identifying variants missing from cache")
    run_shards(
        "miss selection",
        lambda number: (
            f"{bcftools} filter -i 'INFO/{tag}==\"\"' {dirs[number]}/step1.bcf "
            f"-o {dirs[number]}/step2.bcf -Ob -W --threads {shard_threads}"
        ),
    )
    step2_bcf = work_task / f"{sample_name}_isecvst_miss.bcf"
    misses = " ".join(f"{dirs[number]}/step2.bcf" for number in range(len(regions)))
    BcftoolsCommand(
        f"{bcftools} concat --no-version --naive {misses} -Ob -o {step2_bcf} && "
        f"{bcftools} index -f {step2_bcf}",
        manager.logger,
        work_task,
    ).run()
    missing_count = manager._count_variants(step2_bcf, bcftools) or 0
    manager.logger.info(f"Found {missing_count} variants not in cache")

    # Step 3: one recipe invocation for every shard's misses, as on the legacy
    # path; annotator start-up is too expensive to pay per shard.
    step3_count = 0
    dropped_count = 0
    if missing_count > 0:
        manager.logger.info(f"Step 3/4: Annotating {missing_count} missing variants")
        step3_bcf = work_task / f"{sample_name}_missing_annotated.bcf"
        aux_dir = manager.output_dir / "auxiliary"
        aux_dir.mkdir(exist_ok=True)
        anno_cmd = manager._substitute_variables(
            manager.nfa_config_content["annotation_cmd"],
            extra_vars={
                "INPUT_BCF": str(step2_bcf),
                "OUTPUT_BCF": str(step3_bcf),
                "AUXILIARY_DIR": str(aux_dir),
            },
        )
        BcftoolsCommand(anno_cmd, manager.logger, work_task).run()

        stdout_file = work_task / "stdout.txt"
        stderr_file = work_task / "stderr.txt"
        annotation_log = manager.output_dir / "annotation_tool.log"
        if stdout_file.exists() and stdout_file.stat().st_size > 0:
            annotation_log.write_text(stdout_file.read_text())
        annotation_err_log = manager.output_dir / "annotation_tool_err.log"
        if stderr_file.exists() and stderr_file.stat().st_size > 0:
            annotation_err_log.write_text(stderr_file.read_text())
        if aux_dir.exists() and not any(aux_dir.iterdir()):
            aux_dir.rmdir()

        step3_count = manager._count_variants(step3_bcf, bcftools) or 0
        if step3_count < missing_count:
            dropped_count = missing_count - step3_count
            manager.logger.info(
                f"Annotation tool dropped {dropped_count} variants from input. "
                f"By default, these are also removed from cached output to match annotation tool behavior."
            )

        # Step 4: overlay the new annotations back onto each shard's step 1,
        # which already carries the input's order, so no reordering is needed.
        manager.logger.info("Step 4/4: Merging cache and newly annotated variants")
        run_shards(
            "annotation merge",
            lambda number: (
                f"{bcftools} annotate {restrict(number)}-a {step3_bcf} "
                f"{dirs[number]}/step1.bcf -c INFO "
                f"-o {dirs[number]}/step4.bcf -Ob0 -W --threads {shard_threads}"
            ),
        )
        merged = "step4.bcf"
    else:
        manager.logger.info("Step 3/4: Skipped (all variants found in cache)")
        manager.logger.info("Step 4/4: Merging cache and newly annotated variants")
        merged = "step1.bcf"

    # Gather by BGZF block copy: the shards are coordinate-disjoint and were
    # planned in index order, so no record is decoded or recompressed here.
    step4_bcf = work_task / f"{sample_name}_merged.bcf"
    parts = " ".join(f"{dirs[number]}/{merged}" for number in range(len(regions)))
    BcftoolsCommand(
        f"{bcftools} concat --no-version --naive {parts} -Ob0 -o {step4_bcf} && "
        f"{bcftools} index -f {step4_bcf}",
        manager.logger,
        work_task,
    ).run()

    # Publishing performs the single full-size compression of the run.
    output_bcf = manager._resolve_output_bcf(sample_name)
    output_target = output_bcf
    if output_bcf == "-":
        output_target = work_task / f"{sample_name}_vc.bcf"
    select = "" if preserve_unannotated else f"-i 'INFO/{tag}!=\"\"' "
    if preserve_unannotated:
        manager.logger.info(
            "Preserving unannotated variants in output (--preserve-unannotated flag set)"
        )
    if output_bcf == "-":
        publish = f"{bcftools} view {select}{step4_bcf} -o - -Ob --threads {threads}"
    else:
        publish = (
            f"{bcftools} view {select}{step4_bcf} -o {output_bcf} "
            f"-Ob -W --threads {threads}"
        )
    result = BcftoolsCommand(publish, manager.logger, work_task).run()

    manager.last_run_stats = {
        "input_variants": input_count,
        "output_variants": manager._count_variants(Path(output_target), bcftools),
        "missing_variants": missing_count,
        "missing_annotated": step3_count,
        "dropped_variants": dropped_count if missing_count > 0 else 0,
        "engine": "adaptive",
        "shards": len(regions),
        "cache_access": "windows" if windowed else "contigs",
    }
    manager.logger.info(f"Annotation complete: {output_bcf}")

    return result
