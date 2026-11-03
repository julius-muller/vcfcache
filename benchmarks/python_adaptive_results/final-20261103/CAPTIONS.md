# Figure captions

## Figure 1: Annotation reuse with VCFcache

**A**, A compatible recipe-specific cache can be downloaded, or a blueprint can
be annotated with the user's recipe. Input variants are matched to cached variants;
cached annotations are reused, misses are annotated and both streams are combined
into sorted output. **B**, Cache coverage, annotation cost and amortization across
samples govern the benefit. This conceptual figure is copied unchanged from the
archived 2026-08-27 v4 workflow; it contains no replay measurements. VCF, variant
call format. Source: benchmarks/figures/output/final/2026-08-27T225600Z-story-campaign-final-v4/.

## Figure 2: Per-sample runtime reduction across assays and recipes

**A**, Paired direct-to-cached elapsed-time ratios for 12 GRCh38 samples (six KPGP,
six SGDP) per recipe and assay. Each point is an independent sample; colours and
symbols identify cohorts. Black diamonds show medians and bars show 95% bootstrap
intervals from 4,000 resamples within cohort. The dashed line marks parity; values
below one indicate slower cached execution. **B**, Absolute WGS time saved for the
same samples; boxes show quartiles and medians, whiskers extend to 1.5 interquartile
ranges and all samples remain visible. **C**, Serial-processing capacity calculated
as 24 hours divided by the measured elapsed time, not a measured 24-hour service
workload. The diagonal, where visible, denotes equal capacity. Every ratio uses
the same sample's own annotation recipe and matched eight-CPU allocation; recipes
are not compared as biologically equivalent annotation products. Direct baselines
were reused from compatible frozen measurements. Cache setup/build costs are
excluded. PGP uses GRCh37 and appears separately in Figure S2. Panel, targeted-gene
panel; WES, whole-exome sequencing; WGS, whole-genome sequencing; KPGP, Korean
Personal Genome Project; SGDP, Simons Genome Diversity Project.

## Figure 3: Cache strategy and individual WGS benefit

**A**, Cache-hit distributions across 52 genomes (20 KPGP, 20 SGDP and 12 PGP),
shown once using the VEP arm rather than counting annotator arms as extra genomes.
Cache strategies are gnomAD AF ≥ 10%, gnomAD AF ≥ 1% and a three-genome cohort
cache. **B**, Empirical cumulative distributions of per-sample runtime ratios for
each annotator/cohort/strategy combination. Each curve contains 20, 20 or 12 genomes
according to cohort, with 52 genomes per annotator/strategy overall. The fastVEP
recipe here is the original basic campaign recipe, **not** the clinical recipe
in Figure 2. Twelve aliased VEP conditions are views of already measured assay
conditions, not new runs or additional samples. GRCh38 is used for KPGP/SGDP and
GRCh37 for PGP. AF, allele frequency; gnomAD, Genome Aggregation Database; PGP,
Personal Genome Project. Cohort-cache seed separation is distinct from exclusion
of potential gnomAD donors, which was not independently audited here.

## Figure S1: Additional annotation work increases absolute time savings

**A**, Observed direct and cached runtime for one real WGS sample, two AF-threshold
caches and six configured delays in a synthetic VEP plugin. The plugin waits per
consequence and adds no annotation fields. **B**, Corresponding observed runtime
ratios. Points are single descriptive conditions, not technical replicate means;
there are no inferential error bars. Direct times are reused historical measurements;
output checks use a composed VEP semantic reference rather than retained complete
per-delay outputs. Delay models computational work, not ten real database plugins
or an experiment on fastVEP/vep-rs annotation complexity. AF, allele frequency.

## Figure S2: Validated assay results separated by reference assembly

All accepted assay conditions are shown by recipe, assay and assembly, retaining
PGP/GRCh37. Points identify samples and use cohort colours/symbols. The dashed line
marks parity. Missing clinical vep-rs WES/WGS points are excluded because CSQ differs;
four SnpEff conditions are excluded because of the documented ANN warning mismatch.
Identical recipe labels across assemblies do not establish identical information
content or workload. This panel must not be used to attribute assembly differences
solely to annotator efficiency. Definitions follow Figures 2–3.

## Figure S3: Aggregate CPU and memory use

Ratios of cached-to-direct aggregate CPU time and peak memory for accepted assay
conditions with retained paired direct resource measurements. CPU time is the
independent Slurm step's cgroup CPU-usage delta, not wall time or CPU quota. Memory
is aggregate memory.peak from a fresh step for each condition, not a sum of child
RSS values. Boxes show median/quartiles with 1.5-interquartile-range whiskers.
The dashed line denotes equal resource use. Reference assemblies remain separate.

## Diagnostic: Clinical vep-rs annotation differences

This diagnostic is **not a validated speedup figure**. Filled symbols are accepted
core/panel conditions; hollow symbols are the 24 clinical WES/WGS timing ratios
whose CSQ differs from direct annotation. Relevant headers and non-CSQ records
passed comparison. Ratios remain available for engineering assessment but must
not enter headline performance claims without resolving the annotation difference.
CSQ, consequence annotation field.

## Repository headliner: WGS time reduction

Bars show the reciprocal of the median paired speedup for the three annotation-rich
GRCh38 recipes in Figure 2. Labels express that normalized runtime reduction, not
the ratio of median runtimes or a claim about every sample. The grey bar represents
the direct baseline (100%); the coloured bar represents cached runtime. The plot
does not include setup costs or claim that lightweight annotators always benefit.
