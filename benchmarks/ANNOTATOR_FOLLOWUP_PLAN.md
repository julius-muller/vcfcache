# Follow-up benchmark: ANNOVAR, SnpEff, and fastVEP

> **Current priority:** Before this larger publication-quality matrix, run the
> exploratory [fastVEP pilot](fastvep_pilot/README.md). The pilot uses one WGS
> on ITCCcloud_dev to test output equality and current-engine speed at 80%, 90%,
> and 100% hits. If the current engine misses the product gate, use its stage
> profile to decide among the separately scoped native options in
> [FAST_ANNOTATOR_FOLLOWUP.md](fastvep_pilot/FAST_ANNOTATOR_FOLLOWUP.md). Do not
> schedule the matrix below until that decision is available.

## Question and primary comparison

The follow-up asks how much VCFcache reduces the runtime of each annotator, not
which annotator is fastest or which annotator produces the best annotations.
Every headline ratio must therefore compare the same pinned annotator,
configuration, database release, hardware allocation, and input in two modes:

1. annotate the complete input directly;
2. use VCFcache and invoke that exact annotator only for cache misses.

Cross-tool wall times may be shown as context, but they are not VCFcache
speedups because ANNOVAR, SnpEff, fastVEP, and Ensembl VEP do not implement
identical annotation semantics.

## Why a compact campaign is sufficient

Cache hit rate is determined by the input variants and the blueprint, not by
the annotator used to populate that blueprint. The completed real-world WGS
campaign therefore supplies the empirical hit-rate distributions. The new
campaign needs to measure only each pipeline's uncached per-variant cost and
VCFcache lookup/preprocessing overhead, then verify the fitted runtime model on
a small number of independent genomes.

Use the same bundled Zenodo blueprints as the current campaign. Build new,
annotator-specific caches locally from those blueprints; do not treat an
Ensembl-VEP-annotated Zenodo cache as interchangeable with another annotator.

## Integration gate

Before scheduling WGS jobs, prepare one immutable `annotation.yaml` and
`params.yaml` pair per configuration and test it on a small GRCh38 VCF:

- **SnpEff:** VCF input/output, pinned SnpEff and genome database, `ANN` as the
  required annotation tag, `-noLog`, and the same statistics setting in both
  arms.
- **ANNOVAR:** `table_annovar.pl -vcfinput`, pinned script distribution and
  `humandb` files, a stable emitted INFO tag selected from the real output, and
  conversion of the annotated VCF back to BCF.
- **fastVEP:** pinned Git commit and release-mode binary, GFF3, indexed FASTA,
  supplementary annotation files, VCF output with `CSQ`, and conversion back
  to BCF.

For every adapter, require exact per-variant annotation equality between direct
and cached execution after removing only volatile header provenance such as
command paths and timestamps. Record failures by INFO field; never accept only
record-count equality.

## Experimental matrix

### Controlled model

Use one real normalized WGS sample and deterministic blueprints giving observed
hit rates near 50%, 60%, 80%, 90%, and 95%. Run one direct baseline and one
cached job per hit rate. No full technical repeats are needed.

Configurations:

- SnpEff standard gene-effect annotation with HGVS;
- ANNOVAR standard multi-database annotation (gene, region, and filter
  operations);
- fastVEP consequence annotation with HGVS;
- fastVEP plus ClinVar and gnomAD supplementary annotations;
- fastVEP comprehensive configuration, using only databases that can be
  pinned and redistributed or reacquired reproducibly.

This is 30 primary executions: five configurations times one direct baseline
plus five cached hit-rate cells. Run one existing Ensembl-VEP configuration on
the same nodes as a hardware bridge if the follow-up occurs after the current
cluster has changed.

Fit, separately for every configuration:

`T_cached = T_lookup + beta * (1 - f) * T_uncached`

Report measured points, fitted lines, `T_lookup`, and the break-even hit rate.
The break-even hit rate is essential for fastVEP because its annotation stage
may already be shorter than VCFcache preprocessing.

### Real-world validation

Select three evaluation genomes from the existing external cohort before
looking at new runtimes: one near the low, median, and high observed cache-hit
rate. Run each of the five configurations directly and through VCFcache once.
This adds 30 paired executions and tests whether the controlled model predicts
real WGS wall time outside its construction sample.

If the three prediction errors are acceptably small, combine each fitted model
with all observed external-WGS hit rates to estimate a tool-specific speedup
distribution. If not, expand to three genomes per cohort before drawing the
figure.

## Resource and timing controls

- Use one Slurm worker type and identical CPU, RAM, and local-storage limits in
  both arms of every pair.
- Pin actual thread controls and confirm CPU use from cgroup accounting.
- Keep each tool's own indices or internal caches enabled in both arms. In
  particular, fastVEP's binary transcript cache is not a substitute for
  VCFcache and must not be disabled selectively.
- Prebuild databases and indices outside timed regions. Include input
  conversion, VCFcache lookup, miss extraction, annotation, and output merge in
  the end-to-end timed region.
- Balance direct/cached order across the three evaluation genomes. Record wall
  time, CPU time, peak RSS, bytes read/written, hit rate, miss count, and output
  equality.

## Recommended figure

Use one ggplot2 figure with tool/configuration facets:

- x-axis: observed cache-hit rate;
- y-axis: cached/direct wall time;
- points: controlled and real WGS observations with distinct symbols;
- line: fitted runtime model;
- annotation: measured lookup overhead and break-even hit rate.

A small adjacent strip should translate the fitted ratios into minutes saved
for a 10-minute, 1-hour, and 10-hour pipeline. For fastVEP, also show three
clearly separated comparisons: Ensembl VEP versus fastVEP as annotator choice,
fastVEP versus VCFcache plus fastVEP as incremental cache benefit, and the
combined Ensembl-VEP-to-VCFcache-plus-fastVEP ratio only as optional context.

## Primary references

- [fastVEP repository and benchmark protocol](https://github.com/Huang-lab/fastVEP)
- [fastVEP benchmark runner](https://github.com/Huang-lab/fastVEP/blob/master/benchmarks/run_benchmark.sh)
- [fastVEP validation against Ensembl VEP](https://github.com/Huang-lab/fastVEP/tree/master/validation)
- [SnpEff command-line documentation](https://pcingola.github.io/SnpEff/snpeff/commandline/)
- [ANNOVAR quick-start and `table_annovar.pl` examples](https://annovar.openbioinformatics.org/en/latest/user-guide/startup/)
- [ANNOVAR VCF-input documentation](https://annovar.openbioinformatics.org/en/latest/user-guide/input/)
