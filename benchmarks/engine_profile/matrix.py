"""Frozen experiment design, independent of Slurm and annotation tools."""

SAMPLES = ("LP6005519-DNA_A06", "LP6005592-DNA_C05")
RECIPES = ("fastvep", "veprs-core", "vep")
NODES = ("sl-w2", "sl-w4", "sl-w1", "sl-w3", "sl-w6", "sl-w2")


def design():
    """Return the complete bounded matrix; no adaptive scheduling decisions."""
    groups = []
    for sample_index, sample in enumerate(SAMPLES):
        for recipe in RECIPES:
            cells = []

            def add(
                assay,
                engine,
                *,
                purpose="primary",
                cpus=16,
                workers=0,
                memory=512,
                profile=True,
                form="bcf",
                skip_split=False,
                role="measure",
                build="profile",
                perf=False,
                sample=sample,
                recipe=recipe,
                cells=cells,
            ):
                cell = dict(
                    sample=sample,
                    recipe=recipe,
                    assay=assay,
                    engine=engine,
                    purpose=purpose,
                    cpus=cpus,
                    workers=workers,
                    memory_mb=memory,
                    profile=profile,
                    form=form,
                    skip_split=skip_split,
                    role=role,
                    build=build,
                    perf=perf,
                )
                cell["id"] = f"g{len(groups):02d}-c{len(cells):03d}"
                cells.append(cell)

            for assay in ("wes", "wgs"):
                add(
                    assay,
                    "legacy",
                    purpose="reference",
                    role="direct",
                    build="baseline",
                    profile=False,
                )
                add(
                    assay,
                    "adaptive",
                    purpose="reference",
                    role="cached",
                    build="baseline",
                    profile=False,
                )
                order = (
                    ("adaptive", "native")
                    if sample_index == 0
                    else ("native", "adaptive")
                )
                for engine in order:
                    add(assay, engine)
            if sample_index == 0 and recipe in ("fastvep", "veprs-core"):
                assay = "wes" if recipe == "fastvep" else "wgs"
                for cpus in (1, 8):
                    for engine in ("adaptive", "native"):
                        add(assay, engine, purpose="cpu-scaling", cpus=cpus)
                for workers in (2, 4, 8):
                    for engine in ("adaptive", "native"):
                        add(assay, engine, purpose="worker-scaling", workers=workers)
                if recipe == "veprs-core":
                    for memory in (2048, 8192):
                        add("wgs", "native", purpose="memory", memory=memory)
            if sample_index == 0 and recipe == "vep":
                for memory in (2048, 8192):
                    add("wes", "native", purpose="memory", memory=memory)
            if sample_index == 0 and recipe == "fastvep":
                for engine in ("adaptive", "native"):
                    add("wes", engine, purpose="input-format", form="vcf.gz")
                    add("wes", engine, purpose="input-preparation", skip_split=True)
                for purpose in ("first-use-after-copy", "repeated-access"):
                    add("wes", "native", purpose=purpose)
                for pair in range(3):
                    for profile in ((True, False) if pair % 2 == 0 else (False, True)):
                        add(
                            "wes",
                            "native",
                            purpose=f"instrumentation-pair-{pair+1}",
                            profile=profile,
                        )
                for engine in ("adaptive", "native"):
                    add("wes", engine, purpose="stack-sampling", perf=True)
            groups.append(
                dict(
                    id=f"g{len(groups):02d}",
                    sample=sample,
                    recipe=recipe,
                    node=NODES[len(groups)],
                    cells=cells,
                )
            )
    assert sum(len(g["cells"]) for g in groups) <= 96
    return groups
