# Rust-adaptive figure captions

These figures use the frozen Rust-adaptive campaign, not the Python replay and
not the native/resident backend. Source/helper identities and complete timings
are retained in `source/original/`. Original campaign reference limitations and
annotation discrepancies are unchanged.

## Figure 1: Annotation reuse with VCFcache

A, Input variants are matched to a compatible recipe-specific cache; cached
annotations are reused, misses are annotated and results are combined into sorted
output. A recipe cache or blueprint can be obtained from Zenodo. B, Coverage,
annotation cost and amortization govern the benefit. This unchanged conceptual
workflow is copied from the archived 2026-08-27 v4 figure, not generated from new
runtime measurements. VCF, variant call format.

## Figure 2: Rust-adaptive runtime reduction across assays

A, Per-sample direct-to-cached elapsed-time ratios for 12 GRCh38 samples per recipe
and assay (six KPGP and six SGDP). Colours and symbols identify cohorts; black
diamonds show medians and bars show 95% intervals from 4,000 bootstrap resamples
within cohort. Dashed lines mark parity. B, Absolute WGS time saved; boxes show
quartiles/median and whiskers extend to 1.5 interquartile ranges, with every sample
visible. C, Serial processing capacity calculated as 24 hours divided by measured
elapsed time, not a measured 24-hour service workload. The diagonal, where visible,
denotes equal capacity. Each comparison uses its own sample/recipe's direct
denominator with the same eight-CPU allocation. Cache building/download is excluded.
Original paired direct measurements and approved historical references retain their
documented scope; generating these figures reran neither arm. Recipes are not
claimed to have identical biological information content across annotators.
WES, whole-exome sequencing; WGS, whole-genome sequencing; KPGP, Korean Personal
Genome Project; SGDP, Simons Genome Diversity Project.

## Figure 3: WGS cache-strategy coverage and speedup

A, Hit rates across 52 genomes (20 KPGP, 20 SGDP, 12 PGP), shown once using VEP.
B, Empirical cumulative distributions of the paired direct/cached ratios for the
same genomes, by recipe/cohort/cache. The strategies use gnomAD AF ≥ 10%, AF ≥ 1%
or three cohort genomes. The fastVEP recipe here is the basic campaign recipe,
not the clinical recipe of Figure 2. Twelve VEP plotting aliases are reused
measurement views, not additional genomes or runs. KPGP/SGDP use GRCh38 and PGP
uses GRCh37; neither assembly nor recipe information content is treated as
interchangeable. No 1000 Genomes cohort is plotted. Cache-seed separation does not
establish donor exclusion from gnomAD. AF, allele frequency; PGP, Personal Genome
Project; gnomAD, Genome Aggregation Database.

## Figure S1: Additional modeled annotation work

A, Direct and cached runtimes for one real WGS sample, two AF-threshold caches and
six configured wait levels in a synthetic VEP plugin. B, Corresponding runtime
ratios. The plugin waits per consequence and adds no annotations; this is not a
heavy fastVEP or vep-rs experiment. Each point is a single descriptive condition;
no replicate uncertainty is inferred. Direct historical times are retained and
per-delay output checks use a composed VEP semantic reference, not complete
retained per-delay outputs.

## Figure S2: Assembly-separated assay results

All accepted assay measurements are displayed by recipe, assay and assembly,
including PGP/GRCh37. Colours/symbols identify cohorts and the dashed line marks
parity. Twenty-four clinical vep-rs WES/WGS conditions with CSQ differences and
four SnpEff warning-mismatch conditions are excluded. Identical recipe names
across assemblies do not establish identical content or computational cost.

## Figure S3: Aggregate resources

Cached/direct ratios of aggregate CPU time and peak memory for accepted assay
conditions with paired direct resource measurements. CPU time uses the independent
Slurm step's cgroup CPU delta; peak memory is aggregate memory.peak from a fresh
step per condition. Boxes show median/quartiles with 1.5-interquartile-range
whiskers. The dashed line marks equal use; assemblies remain separate.

## Repository headliner: WGS time reduction

Normalized bars show the reciprocal of the median paired speedup for three
annotation-rich GRCh38 recipes; labels express the corresponding relative runtime
reduction, not the ratio of median times or a guarantee for every sample. Grey
represents direct annotation (100%); colour represents cached annotation. Setup
costs are excluded. Lightweight workflows need not benefit.

For the clinical vep-rs diagnostic and provisional Figure 4, hollow points indicate
CSQ differences and are timing evidence only. The full caption and hypothetical
Figure 5 model specification are in
`../../python_adaptive_results/interpretation-20261103/README.md`; apply the Rust
timings reported in `../../rust_adaptive_results/README.md`. These modeled panels
are identical across implementations and must not be presented as measured Rust
CPU scaling or a fitted performance predictor.
