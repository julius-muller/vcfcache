# Findings

## VCFcache reuses complete annotation recipes

Variant annotation repeatedly evaluates alleles shared among samples. Local
VEP or SnpEff data caches avoid repeated database downloads, but each invocation
still computes every requested consequence [@mclaren2016vep;
@cingolani2012snpeff]. Field-transfer tools can copy selected columns from a
donor VCF, but do not define how missing alleles traverse a composite pipeline
[@danecek2021bcftools; @vcfanno2016]. VCFcache instead stores the output of a
fixed recipe: annotator version, flags, plugins and supplementary databases.
For each new sample, cached alleles are retrieved, misses pass through that
same recipe, and both paths are assembled as a sorted, indexed BCF (Figure 1).

![](figures/main/figure1_workflow_v3.png){width=100%}

**Figure 1: Recipe-specific cache construction and per-sample annotation
reuse.** (A) Users may download a ready-to-use cache when build, annotator and
recipe match; annotate a public blueprint once with a custom recipe; or build a
blueprint from cohort variants. (B) A new sample is normalized and divided by
cache lookup. Stored annotations are returned for hits, only misses are
annotated, and both paths are merged into a coordinate-sorted, indexed BCF.
BCF, binary variant call format.

For $N$ variants, fixed annotator start-up $s_{\mathrm{ann}}$, per-variant cost
$t_{\mathrm{ann}}$, hit rate $f$ and VCFcache processing cost $C(N)$,

$$
T_{\mathrm{direct}}=s_{\mathrm{ann}}+Nt_{\mathrm{ann}}, \qquad
T_{\mathrm{cached}}\approx s_{\mathrm{ann}}+(1-f)Nt_{\mathrm{ann}}+C(N).
$$

Caching helps when $fNt_{\mathrm{ann}}>C(N)$. It can return per-variant work,
but not annotator start-up; consequently, high coverage alone does not
guarantee acceleration for very small inputs or minimal recipes. A custom
cache's one-time build cost is amortized as
$T_{\mathrm{effective}}\approx T_{\mathrm{cached}}+T_{\mathrm{build}}/S$
after $S$ uses.

## Public and cohort caches accelerate held-out genomes

We evaluated 20 Korean Personal Genome Project (KPGP), 20 Simons Genome
Diversity Project (SGDP) and 12 Harvard Personal Genome Project (PGP) genomes
[@kpgp_data; @mallick2016sgdp; @ball2012pgp]. KPGP and SGDP remained in
GRCh38, PGP in GRCh37. Three additional genomes from each provider formed
disjoint cohort caches; relatedness screening prevented a second-degree-or-
closer training/evaluation pair. The public alternatives were bundled,
assembly-matched gnomAD v4.1 caches containing alleles reaching AF ≥10% or
AF ≥1% in any reported frequency stratum [@gnomad2024release].

With VEP 115.2 `--everything`, median end-to-end speedups across 52 genomes
were 4.75-fold (AF ≥10%), 9.21-fold (AF ≥1%) and 4.59-fold (cohort cache).
With fastVEP 0.3.0 core consequence/HGVS recipes, they were 1.84-, 2.05- and
1.81-fold. Median AF ≥1% coverage was stable across KPGP, SGDP and PGP
(91.95–93.04%), whereas three-genome cohort coverage ranged from 77.45% to
82.88% (Supplementary Figure 1). The cohort caches repaid their measured build
cost after two subsequent genomes; bundled caches require no local build.

## Assay size and annotation load set the break-even point

Twelve independent GRCh38 genomes were represented as matched ACMG panel, WES
and WGS inputs. At similar median coverage (90.35–93.23%), VEP speedups were
0.71-, 2.29- and 9.50-fold, whereas core fastVEP speedups were 0.31-, 0.54-
and 2.02-fold (Figure 2). Panel and core-fastVEP exome cells were dominated by
fixed workflow or annotator costs; genome-scale work was large enough to
recover them.

![](figures/main/figure2_assay_annotator_v3.png){width=100%}

**Figure 2: Recoverable annotation work sets the cache break-even point.**
Twelve independent GRCh38 genomes (six KPGP and six SGDP) were evaluated as
matched ACMG SF v3.3 Panel, Twist Human Core Exome WES and WGS inputs with the
gnomAD AF ≥1% strategy. (A) Median direct runtime (open circle) and cached
runtime (arrowhead) on a logarithmic scale. Green arrows indicate time returned
and muted red arrows additional fixed workflow time. (B) Median speedup and
percentage wall-time change. Each sample/tool/assay cell contributed one
direct and one cached measurement. All outputs passed annotator-specific
semantic comparison. AF, allele frequency; KPGP, Korean Personal Genome
Project; SGDP, Simons Genome Diversity Project; VEP, Variant Effect Predictor;
WES, whole-exome sequencing; WGS, whole-genome sequencing.

A second campaign tested fastVEP and SnpEff 5.4c across 12 genomes per assembly
and the same three assay scales. Clinical SnpEff comprised consequence
annotation followed by ClinVar and gnomAD exome annotation. Its median speedups
were 2.66- and 3.62-fold for WES and 10.78- and 12.31-fold for WGS on GRCh38
and GRCh37, respectively. On the same GRCh37 genomes and cache membership,
moving from core to clinical SnpEff increased WES speedup from 1.16- to
3.62-fold and WGS speedup from 4.46- to 12.31-fold (Supplementary Table 1).
Clinical fastVEP reached 2.54-fold (GRCh38) and 1.70-fold (GRCh37) at WGS.
These are within-tool comparisons: native recipes differ in output scope and
do not rank annotators against one another.

An engineering experiment further tested whether additional CPUs eliminate
the benefit for fast native annotation. A dense ten-database fastVEP recipe
saturated between process-wide limits of 10 and 32 CPUs (706.2 versus 691.4 s
direct) but retained 3.07- and 3.37-fold speedups at 90% cache hits. At 10 CPUs,
caching returned 7.9 minutes per genome (Supplementary Figure 2). Parallelism
and reuse therefore reduce different parts of runtime.

## Cache coverage controls relative savings; pipeline cost controls hours saved

One held-out 4.80-million-variant KPGP genome was annotated with six controlled
VEP recipes spanning 4.2–22.7 hours. A no-output plugin added 0.5–10 ms per
transcript consequence without changing annotations. Bundled caches covered
80.23% and 90.26% of variants. The former left 19.7–21.7% of direct wall time
(4.62–5.07-fold speedup); the latter left 9.9–11.9% (8.39–10.10-fold). For the
longest recipe, 22.71 hours fell to 4.52 and 2.25 hours (Figure 3).

![](figures/main/figure3_pipeline_complexity_v3.png){width=100%}

**Figure 3: Absolute time returned grows while relative speedup approaches a
hit-rate ceiling.** One held-out KPGP WGS was processed with six VEP pipeline
loads and bundled AF ≥10% and AF ≥1% caches. (A) Direct runtime (open circles),
cached runtime (arrowheads) and wall time returned. (B) Speedup; dashed lines
show the asymptotic $1/(1-f)$ ceilings for observed hit rate $f$. The synthetic
plugin emitted no annotation field. All 12 cached outputs passed semantic
validation. AF, allele frequency; VEP, Variant Effect Predictor; WGS,
whole-genome sequencing.

Measured and predicted relative runtimes differed by at most 0.56 percentage
points. Thus relative speedup plateaus near $1/(1-f)$ even while absolute time
saved continues to increase—a distinction useful to both high-throughput
services and laboratories evaluating a particular recipe.

## Correctness was assessed outside timed cells

All 156 VEP and 156 fastVEP whole-genome pairs passed annotator-specific
semantic comparison. For VEP, 8,898 of 734,582,949 records (0.0012%) differed
only in `HGNC_ID`, a documented VEP 115.2 batching effect [@vep_issue_1959];
this field was reported but excluded from pass/fail. fastVEP comparison covered
all INFO, FORMAT and relevant header values, canonicalizing order only.

In the final fastVEP/SnpEff campaign, 176 of 180 pairs passed. Four SnpEff
WES/WGS cells from one GRCh37 genome were conservatively excluded because one
allele differed only by `WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS`; consequence,
transcript, HGVS, ClinVar, gnomAD, genotype and record counts agreed. A
two-record direct-SnpEff reproducer showed that the warning appears when an
adjacent insertion is processed first and disappears when record order is
reversed. VCFcache was absent from this reproducer, so no warning-specific
normalization was added. These results support equivalence under a fixed tool,
version, reference and recipe, not byte identity or equivalence between
annotators.
