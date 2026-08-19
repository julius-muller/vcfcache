# Findings

## Repeated variant annotation performs redundant work

Functional annotation is a standard step in sequencing analysis, but its cost
is repeatedly paid for variants shared among samples. Ensembl Variant Effect
Predictor (VEP), for example, locally caches reference and transcript data but
still evaluates every input allele for every invocation
[@mclaren2016vep]. This distinction matters in routine laboratories and
services, where the same annotation recipe may be applied to many incoming
samples. Faster native implementations such as fastVEP reduce the cost of an
individual annotation, but they do not remove computation already performed
for a recurring allele [@fastvep2026].

Reuse of precomputed annotation is not itself new, but existing mechanisms
reuse *inputs*, not the output of a pipeline. VEP's local cache stores reference
and transcript data, removing repeated database access while still evaluating
every input allele [@mclaren2016vep]. Field-transfer tools such as `bcftools
annotate` and vcfanno copy columns from a donor file onto matching positions
[@danecek2021bcftools; @vcfanno2016]. Both are limited to what a single
distributable resource already contains. Production annotation is instead a
*composite recipe*: a specific tool version with specific flags, plugins,
supplementary databases and custom tracks, often chained through several
processes. No public resource exists for the output of such a recipe, because
the recipe is site-specific, and a donor file also has no defined behaviour for
variants it lacks, so novel alleles are silently left unannotated.

VCFcache caches the result of exactly that composite recipe. The recipe is
stored with the cache and is immutable, so every cached record is reproducible
from a known command; variants absent from the cache are routed through the same
recipe rather than dropped; and cached and direct output are required to be
annotation-equivalent, which is verified here for every benchmark run. The unit
of reuse is therefore the pipeline's output rather than any one database, which
is what makes the approach applicable to arbitrary annotators.

VCFcache implements annotation reuse as a layer around the annotator. A
blueprint contains normalized, deduplicated alleles from a public resource or
a user cohort. The blueprint is annotated once with an immutable recipe, and
the resulting indexed BCF becomes an annotation cache. For a new sample,
VCFcache retrieves hits from that cache, sends only misses through the same
annotator recipe, and assembles a complete sorted and indexed output (Figure
1). The annotator remains authoritative for uncached variants. The command
recipe is therefore the integration boundary rather than a tool-specific
reimplementation of annotation semantics.

![](figures/main/figure1_workflow_v3.png){width=100%}

**Figure 1: Recipe-specific cache construction and per-sample annotation
reuse.** (A) Users can download a ready-to-use cache from Zenodo when genome
build, annotator, version and recipe match; annotate a downloaded blueprint
once with a custom fixed recipe; or construct a blueprint from existing cohort
variants before annotation. All paths produce a recipe-specific annotation
cache. (B) Each new sample is normalized and divided during cache lookup.
Cached annotations are reused for hits, only cache misses are submitted to the
same recipe, and both paths are merged into a coordinate-sorted, indexed BCF.
BCF, binary variant call format.

An annotator invocation has a fixed start-up cost $s_{\mathrm{ann}}$, paid once
to load transcript and reference data, and a per-variant cost
$t_{\mathrm{ann}}$. For a sample of $N$ variants,
$T_{\mathrm{direct}} = s_{\mathrm{ann}} + N t_{\mathrm{ann}}$. A cached run
still invokes the annotator once, for the misses, so with cache hit rate $f$

$$
T_{\mathrm{cached}} \approx s_{\mathrm{ann}} + (1-f) N t_{\mathrm{ann}}
+ C_{\mathrm{VCFcache}}(N),
$$

where $C_{\mathrm{VCFcache}}$ covers input processing, cache lookup, miss
selection and output construction. The start-up term is identical in both
expressions and therefore cancels: **caching can return only the per-variant
annotation work, never the cost of starting the annotator.** Reuse saves time
exactly when

$$
f N t_{\mathrm{ann}} > C_{\mathrm{VCFcache}}(N),
$$

that is, when the observed hit rate exceeds the break-even hit rate
$C_{\mathrm{VCFcache}} / (N t_{\mathrm{ann}})$. Relative speedup is therefore
governed by hit rate, while whether caching helps at all is governed by how much
per-variant work the annotation recipe performs. The one-time build cost is
amortized across $S$ uses:

$$
T_{\mathrm{effective}} \approx T_{\mathrm{cached}} + \frac{T_{\mathrm{build}}}{S}.
$$

## Public and cohort caches accelerate held-out whole genomes

We tested the model with 52 publicly accessible whole genomes that were not
used to construct the cohort caches: 20 Korean Personal Genome Project
(KPGP), 20 Simons Genome Diversity Project (SGDP), and 12 Harvard Personal
Genome Project (PGP) genomes [@kpgp_data; @mallick2016sgdp; @ball2012pgp]. KPGP and SGDP
inputs remained in GRCh38, whereas PGP inputs remained in their documented
GRCh37 coordinates. Three additional genomes from each provider cohort were
assigned to cache construction and were disjoint from evaluation. Relatedness
screening prevented a second-degree-or-closer relationship across that
training/evaluation boundary; one realistic related pair was retained within
the KPGP evaluation set.

Every evaluation genome was annotated directly and with three strategies: the
bundled, assembly-matched gnomAD v4.1 caches containing sites whose allele
frequency reached at least 10% or 1% in any reported frequency stratum, and the
corresponding three-genome cohort cache [@gnomad2024release]. The inclusion rule
was evaluated across all strata, not only the combined-population frequency. No
documented project overlap with gnomAD was found for these cohorts, although
undisclosed sample-level overlap cannot be excluded.

For VEP 115.2, median end-to-end speedups across the 52 genomes were
4.75-fold with the any-stratum AF ≥10% cache, 9.21-fold with AF ≥1%, and
4.59-fold with the cohort cache (Supplementary Figure 1). The deeper public
cache produced higher hit rates and left the smallest share of direct wall time
(Supplementary Figure 1).
The cohort cache nevertheless provided a similar median benefit to the
smaller public cache despite being constructed from only three genomes. This
supports both usage modes: a downloadable public cache gives an immediate
starting point, while local cohort history can be converted into a
recipe-specific cache.

Public-cache coverage was stable across cohorts differing substantially in
ancestry: median AF ≥1% hit rates were 91.95% (KPGP), 93.04% (SGDP, selected for
global diversity) and 92.66% (PGP). The three-genome cohort caches varied far
more (77.45–82.88%). A frequency-threshold public cache therefore delivers
comparable benefit across ancestries, whereas a small cohort cache reflects
whoever contributed to it.

Cache construction is a one-time cost that these data allow us to quantify
rather than assume. Building a three-genome cohort cache took 13,374, 14,229
and 6,666 seconds for KPGP, SGDP and PGP, against median per-sample savings of
7,503, 7,333 and 3,381 seconds. Each cohort cache therefore repays its
construction after **two** annotated genomes, and every genome thereafter is net
benefit. The bundled public caches remove this cost entirely.

## Caching remains useful around a fast native annotator

The same design was repeated with fastVEP 0.3.0. fastVEP is used here as an
example of a fast, natively compiled annotation pipeline rather than as a
recommended or widely adopted tool; it is a recent community reimplementation
with limited uptake. Its role is to represent the regime in which per-variant
annotation cost is very low, so that the break-even boundary of the caching
layer becomes observable at all. Throughout, fastVEP was run with a core
consequence and HGVS recipe carrying **no supplementary annotation databases**,
which is the minimum-work configuration rather than a production pipeline.

The annotations stored in each cache were rebuilt with the assembly-matched
fastVEP recipe; VEP-derived annotations were never presented as fastVEP cache
records. Median speedups were 1.84-fold, 2.05-fold and 1.81-fold for the
AF ≥10%, AF ≥1% and cohort strategies, respectively (Supplementary Figure 1).
Absolute direct and cached times were substantially shorter than for VEP, but
reuse still returned material time for every strategy.

The lower relative fastVEP speedups for PGP in Supplementary Figure 1 follow the
same rule rather than poorer coverage: median direct runtime was 207 seconds for
the GRCh37 PGP configuration against 402–427 seconds for the GRCh38 cohorts,
leaving less per-variant work to recover against an unchanged workflow cost.
This is an assembly and recipe interaction, not a cohort effect.

The matched assay extension isolated that fixed-cost boundary while retaining
real sample variation. With VEP, median direct and cached runtimes were 29.6
and 39.8 seconds for Panel, 407.5 and 177.6 seconds for WES, and 9,305 and 978
seconds for WGS, corresponding to 0.71-, 2.29- and 9.50-fold speedups. With
fastVEP, the corresponding speedups were 0.31-, 0.54- and 2.02-fold. Similar median cache coverage across
Panel, WES and WGS (90.35–93.23%) demonstrates why hit rate alone is
insufficient.

Decomposing each cell accounts for all six outcomes. Estimating
$s_{\mathrm{ann}}$ from the Panel condition, whose 255 median variants make
per-variant work negligible, gives 29.6 seconds for VEP and 9.0 seconds for
fastVEP. Comparing recoverable work with $C_{\mathrm{VCFcache}}$ then predicts
a gain in exactly the four cells that gained and a loss in the two that lost,
and leave-one-out refitting recovers all six outcomes without the held-out
cell.

For fastVEP the fitted cost is
$C_{\mathrm{VCFcache}}(N) = 21.7\ \mathrm{s} + 27.5\ \mathrm{\mu s} \times N$,
with residuals of at most 2.6 seconds across a four-order-of-magnitude range of
$N$. An independent diagnostic agrees: at an engineered 100% hit rate, where all
cached runtime is by definition VCFcache's own cost, a 4.33-million-variant
genome took 136.9 seconds, or 27 µs per variant, capping the core fastVEP
recipe at 2.84-fold genome-wide.

Two consequences follow. First, Panel cells cannot break even at any hit rate,
because a 255-variant input contains essentially no per-variant work to recover.
Second, the fastVEP WES cell is not a statement about fastVEP: its core recipe
performs only 10.6 seconds of per-variant work on an exome, so reuse could
return at most 9.9 seconds against a 26.5-second workflow cost.
This is a property of the recipe, not of the assay or the annotator, and the
same arithmetic explains why VEP — which performs 378.0 seconds of per-variant
work on the same exome — accelerates comfortably.

We tested that interpretation by enriching the fastVEP recipe with ClinVar,
REVEL and gnomAD v4.1 exome frequencies and repeating the measurement on the
same twelve genomes. Median direct runtime rose from 19.6 to 71.0 seconds and
the cached path moved from 0.54-fold to **1.10-fold**, with ten of twelve
genomes faster than direct and identical variant and consequence output in all
twelve. Enrichment does carry a fast annotator across the break-even boundary at
exome scale.

The size of that gain is the more informative result. Caching returned a median
of 6.3 seconds per exome, close to the 3.7 seconds of recoverable per-variant
work predicted from the observed miss count. Decomposing the enriched recipe
shows why: it costs 53.6 µs per variant on a 61.7-second fixed start-up, so
start-up is 94% of the run and is the one component caching never returns.
Database type therefore matters more than database count — allele-level
databases add per-variant work, whereas gnomAD exome frequencies are dominated
by load cost, and enrichment raises the recoverable and unrecoverable terms
together. A very fast annotator on a small assay is start-up-bound, so caching
is most valuable where per-variant work dominates: whole genomes for any
annotator, and exomes for annotators whose per-variant cost is substantial.

![](figures/main/figure2_assay_annotator_v3.png){width=100%}

**Figure 2: Per-variant annotation work, not assay label, sets the cache
break-even point.** Twelve independent GRCh38 genomes (six KPGP and six SGDP)
were evaluated as matched ACMG SF v3.3 Panel, Twist Human Core Exome WES and WGS
inputs with the gnomAD AF ≥1% strategy. fastVEP was run with a core consequence
and HGVS recipe carrying no supplementary annotation databases, which is its
minimum-work configuration. (A) Open circles denote median direct
runtime and arrowheads median cached runtime on a logarithmic scale. Green
arrows indicate time returned; muted red arrows indicate additional fixed
workflow time. (B) Each tile reports median speedup and the percentage by which
cached wall time was shorter or longer for one assay/annotator combination.
Exactly one direct and one
cached measurement were retained per sample, tool and assay after a separate
untimed calibration input established a consistent warm state. Every cached
output passed the annotator-specific semantic comparator. ACMG, American
College of Medical Genetics and Genomics; AF, allele frequency; KPGP, Korean
Personal Genome Project; SGDP, Simons Genome Diversity Project; VEP, Variant
Effect Predictor; WES, whole-exome sequencing; WGS, whole-genome sequencing.

The smaller fastVEP speedups expose the fixed-cost boundary of the current
implementation. Stage-resolved profiling of a 4.33-million-variant genome at
100% hits attributes the residual to cache lookup (24.5 s), miss selection
(20.1 s) and final output construction (21.6 s), with negligible subprocess
start-up. Compression dominates: writing the looked-up records and the final
output each cost over 180 CPU-seconds for roughly 522 MB. The workflow performs
three full compress-and-index passes over the annotated payload, so
$C_{\mathrm{VCFcache}}$ scales with stored annotation size, not variant count
alone, and a richer recipe raises both the work caching avoids and the cost of
the caching layer. These passes, rather than the recipe, are the target for a
future streaming implementation. A separate held-out-WGS experiment used a dense ten-database recipe under
process-wide CPU affinity. Direct runtime changed little between 10 and 32 CPUs
(706.2 and 691.4 seconds), indicating saturation, whereas a controlled 90% cache
gave 3.07-fold and 3.37-fold speedups. Additional CPUs did not make annotation
complexity free: at 10 CPUs the dense recipe remained 52.2% slower than the core
recipe, and caching returned 7.9 minutes per genome (Supplementary Figure 2).
Reuse and parallelism therefore address different components of runtime.

These data do not establish fastVEP as a semantic replacement for VEP. VCFcache
was evaluated within each annotator, so cross-annotator differences lie outside
the caching equivalence claim.

## Cache coverage and pipeline cost determine the user-visible benefit

To isolate pipeline cost from sample-to-sample variation, we used one held-out
GRCh38 KPGP genome containing 4,795,706 variants and two bundled Zenodo caches
[@vcfcache_af10_cache; @vcfcache_af1_cache]. The observed hit rates were
80.23% and 90.26%. Six VEP --everything recipes added 0.5–10 ms of work per
transcript consequence through a no-output plugin. The plugin changed runtime
but emitted no annotation fields, allowing the same annotation output to be
compared across controlled pipeline loads.

The AF ≥10% cache left 19.7–21.7% of direct wall time (4.62–5.07-fold) and the
AF ≥1% cache 9.9–11.9% (8.39–10.10-fold; Figure 3B). Absolute savings grew with
pipeline cost: the longest direct run took 22.71 hours, reduced to 4.52 and
2.25 hours respectively (Figure 3A). The same hit rate is therefore best read
two ways: relative speedup communicates throughput, hours returned communicates
practical value.

![](figures/main/figure3_pipeline_complexity_v3.png){width=100%}

**Figure 3: Absolute time returned grows while relative speedup approaches a
hit-rate ceiling.** One held-out 4.80-million-variant KPGP WGS was annotated
with six controlled VEP pipeline loads and two bundled caches. (A) Open circles
mark direct runtime and arrowheads mark cached runtime; labels report wall time
returned per genome. (B) Corresponding speedup over direct annotation. Dashed
lines show the asymptotic ceilings $1/(1-f)$ implied by the observed hit rates
$f$ of 80.23% and 90.26%. Synthetic delays emitted no annotation field and
affected only variants submitted to VEP. All 12 cached outputs passed semantic
validation. AF, allele frequency; VEP, Variant Effect Predictor; WGS,
whole-genome sequencing.

As direct annotation becomes more expensive, fixed overhead becomes a smaller
fraction of cached runtime and speedup approaches $1/(1-f)$, so it plateaus even
as absolute time returned keeps growing. Across all 12 cache-by-load
observations, predicted and measured relative runtimes differed by no more than
0.56 percentage points, supporting estimation for a user's own pipeline from
three locally measurable quantities: direct runtime, cache hit rate and fixed
VCFcache cost.

## Cached outputs preserve annotator-specific results

Correctness comparison was performed after, and excluded from, every timed cell.
All 156 VEP cached outputs passed a streaming comparator that
checked variant keys, selected input fields, genotype, the CSQ definition and
the complete set of CSQ entries, canonicalizing semantically irrelevant ordering
and excluding from pass/fail only the documented VEP 115.2 HGNC_ID batching
discrepancy. Across the 156 external-WGS comparisons, 8,898 of 734,582,949
compared records (0.0012%) differed only in HGNC_ID, at most 495 records in any
one output. A separate 89-output validation set showed the same pattern, with no
unexpected annotation mismatch.
All 156 fastVEP cached outputs passed a
stricter complete-record and relevant-header comparison, including
supplementary INFO and FORMAT fields, with only ordering canonicalized. All 12
outputs in the pipeline-complexity experiment also passed their specified
semantic gate.

The benchmarks have three boundaries. First, the claim is annotation
equivalence under a fixed tool, version, reference and recipe, not byte identity
or equivalence between annotators. Second, reported speedups include
normalization, lookup, annotation of misses, assembly, compression and indexing,
but exclude cache construction and the comparator. Third, construction is a
one-time cost whose value depends on reuse: for few samples or frequently
changing recipes, direct annotation may be simpler.

Two adoption routes have different economics. A **downloadable prebuilt cache**
costs only transfer and disk: the GRCh38 VEP caches are 4.03 GB and 1.74 GB as
distributed archives, covering 18.9 and 8.1 million alleles. There is no build
step and no amortization question, provided the local tool version matches.

A **custom cache** additionally costs one annotation pass over the blueprint.
Building each three-genome cohort cache took 1.85–3.95 hours, which is about
1.4 direct whole-genome VEP annotations, and returned about 2.1 hours per
subsequent genome. Each therefore repaid itself after **two** genomes and
reached a 2.9-fold effective per-sample speedup by ten. Neither route is
constrained by memory: peak resident memory during annotation was lower with
caching than without, because fewer records reach the annotator.
