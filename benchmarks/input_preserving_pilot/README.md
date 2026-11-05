# Input-preserving main / bcftools 1.25 pilot

The simple main workflow now uses the supplied sorted/indexed sample directly.
It never splits alleles or removes records. Cache/blueprint normalization remains
separate. Exact biallelic matches are reused; multiallelics are fresh misses.
Duplicate variant keys explicitly bypass cache reuse because bcftools VCF-to-VCF
pairing is one-to-one. The input is never modified. Recipe record identity and
FORMAT validation is included in elapsed time; it is not hidden preparation.

## Design

Use only Intel Broadwell workers **sl-w1/2/3/4/6**. **sl-w5 is excluded** from
all benchmark and review jobs. Each block has a fresh eight-CPU Slurm step per
condition, 28 GB memory, local scratch and a two-hour safety limit. The five
blocks start concurrently. Reader controls follow core WES on sl-w1; idle nodes
after their assigned work is finished are intentional.

One GRCh38 KPGP-00139 input, with the existing frozen WES target definition and
whole WGS, is tested on vep-rs core / fastVEP clinical (WES/WGS) and vep-rs
clinical (WES only). No full-cohort campaign and no technical best-picking.
Assay extraction/BCF encoding is untimed and does not normalize or remove sites.
Expected record counts must match the frozen definition or the block fails.

Each block times once: original main + 1.21 direct/cached, revised main + 1.25
direct/cached, and revised main + 1.25 cached with merged occupied 16-kb regions.
Mode order rotates across workers. Old/new comparisons include code and
dependency changes together; the separate reader control isolates access mode
and version. Zlib-based codec features are matched; no libdeflate A/B is added.
CPU allocation is enforced by Slurm, not just thread-count settings. Thread
wrappers cover supported commands inside the unchanged annotation recipes,
including output conversion/indexing and 1.25 sorting. Existing thread options
are not overridden. No redundant sorting is introduced for indexable outputs.

Semantic comparisons run only after every timed condition in a block. Complete
INFO/FORMAT and relevant headers are compared without ignored fields. Known
vep-rs clinical mismatches remain diagnostic, not validated speedups. Raw timings
and comparison results are retained even when validation fails.

The pilot-only `VCFCACHE_PILOT_REGION_MODE=windows` applies narrow `-R` regions to
both overlays; window construction is timed. It is off by default and adds no
public engine selection. Source/cache/recipe hashes, hardware metadata, commands,
resource counters, states and job IDs are durable. Large disposable outputs stay
in owned worker scratch for later audit; do not delete them automatically after
a failed or wall-limited block.

## Review

An after-any Slurm collector writes `REVIEW.json`, `conditions.tsv` and
`RESULTS.md`; it never submits another campaign. Review remains possible without
Slurm accounting days later. Preserve the original main source deployment and
the private 1.25 tools alongside the manifest/READY binding.

```sh
python -m benchmarks.input_preserving_pilot.run collect RESULT_ROOT
Rscript --vanilla benchmarks/input_preserving_pilot/render.R LOCAL_RESULT_ROOT figures-v1
```

The R/ggplot2 renderer refuses existing output folders, writes separate panels
and complete figures as PDF/SVG/600-dpi PNG, and saves the data as named tibbles.
Figure 1 shows single-run elapsed times per same-node condition; facet y scales
differ. These are engineering observations, not a cohort uncertainty estimate.
