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

If a sample contains $N$ variants and the average annotation cost is
$t_{\mathrm{ann}}$, direct runtime is approximated by
$T_{\mathrm{direct}} = N t_{\mathrm{ann}}$. With cache hit rate $f$, cached
runtime is approximated by

$$
T_{\mathrm{cached}} \approx T_{\mathrm{overhead}} + (1-f)T_{\mathrm{direct}},
$$

where $T_{\mathrm{overhead}}$ includes input processing, cache lookup and
output construction. The one-time build cost is amortized across $S$ uses:

$$
T_{\mathrm{effective}} \approx T_{\mathrm{cached}} + \frac{T_{\mathrm{build}}}{S}.
$$

The model predicts that hit rate primarily determines relative speedup,
whereas the absolute time returned to the user grows with annotator cost.

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
frequency reached at least 10% or 1% in any reported frequency stratum, and
the corresponding three-genome cohort cache [@gnomad2024release]. The public
cache wording is important: the inclusion rule was evaluated across all
frequency strata and was not limited to the combined-population frequency.
No documented project overlap with gnomAD was found for these provider
cohorts, although undisclosed sample-level overlap cannot be excluded.

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

## Caching remains useful around a fast native annotator

The same design was repeated with fastVEP 0.3.0. The annotations stored in
each cache were rebuilt with the assembly-matched fastVEP recipe; VEP-derived
annotations were never presented as fastVEP cache records. Median speedups
were 1.84-fold, 2.05-fold and 1.81-fold for the AF ≥10%, AF ≥1% and cohort
strategies, respectively (Supplementary Figure 1). Absolute direct and cached
times were substantially shorter than for VEP, but reuse still returned
material time for every strategy.

The lower relative fastVEP speedups visible for PGP in Supplementary Figure 1
do not track poorer public-cache coverage. Median direct fastVEP runtime was
207 seconds for the GRCh37 PGP configuration versus 402–427 seconds for the
two GRCh38 cohorts, while fixed parsing, lookup and output work did not halve.
The faster direct baseline therefore leaves a larger fixed-cost fraction; this
assembly/recipe interaction should not be interpreted as a cohort effect.

The matched assay extension isolated that fixed-cost boundary while retaining
real sample variation. With VEP, median direct and cached runtimes were 29.6
and 39.8 seconds for Panel, 407.5 and 177.6 seconds for WES, and 9,305 and 978
seconds for WGS, corresponding to 0.71-, 2.29- and 9.50-fold speedups. With
fastVEP, the corresponding speedups were 0.31-, 0.54- and 2.02-fold. Thus the
current workflow is counterproductive for very small panels and for an
ultrafast WES recipe, but materially accelerates WGS for both annotators
(Figure 2). Similar median cache coverage across Panel, WES and WGS
(90.35–93.23%) demonstrates why hit rate alone is insufficient: direct
annotation cost must exceed the fixed lookup and output cost before reuse
returns time.

![](figures/main/figure2_assay_annotator_v3.png){width=100%}

**Figure 2: Assay scale and annotator speed define the cache break-even
point.** Twelve independent GRCh38 genomes (six KPGP and six SGDP) were
evaluated as matched ACMG SF v3.3 Panel, Twist Human Core Exome WES and WGS
inputs with the gnomAD AF ≥1% strategy. (A) Open circles denote median direct
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
implementation. Identical cache hits remove a larger share of total work from
the slower VEP pipeline, while input parsing, lookup, miss selection and
output construction remain for both annotators. A separate held-out-WGS
experiment used a dense ten-database recipe under process-wide CPU affinity.
Direct runtime changed little between 10 and 32 CPUs (706.2 and 691.4
seconds), indicating saturation, whereas a controlled 90% cache reduced these
times to 230.0 and 205.4 seconds: 3.07-fold and 3.37-fold speedups,
respectively. Ten CPUs made the dense direct recipe 18.7% faster than the core
recipe restricted to one CPU, showing that additional hardware can offset
annotation complexity. However, at the same 10-CPU allocation the dense recipe
remained 52.2% slower than the core recipe, and caching returned 7.9 minutes
per genome (Supplementary Figure 2). This is treated as an engineering stress
test rather than an estimate of a typical annotation profile.

These data do not establish fastVEP as a drop-in semantic replacement for
VEP. VCFcache was evaluated within each annotator: cached VEP was compared
with direct VEP, and cached fastVEP with direct fastVEP. Cross-annotator
differences are outside the caching equivalence claim.

## Cache coverage and pipeline cost determine the user-visible benefit

To isolate pipeline cost from sample-to-sample variation, we used one held-out
GRCh38 KPGP genome containing 4,795,706 variants and two bundled Zenodo caches
[@vcfcache_af10_cache; @vcfcache_af1_cache]. The observed hit rates were
80.23% and 90.26%. Six VEP --everything recipes added 0.5–10 ms of work per
transcript consequence through a no-output plugin. The plugin changed runtime
but emitted no annotation fields, allowing the same annotation output to be
compared across controlled pipeline loads.

The AF ≥10% cache left 19.7–21.7% of direct wall time and produced
4.62–5.07-fold speedups. The AF ≥1% cache left 9.9–11.9% and produced
8.39–10.10-fold speedups (Figure 3B). Absolute savings grew with pipeline
cost: the longest direct run took 22.71 hours and was reduced to 4.52 hours
with the AF ≥10% cache and 2.25 hours with AF ≥1% (Figure 3A). Thus the same
cache hit rate can be presented in two complementary ways: relative speedup
communicates throughput, whereas hours returned communicates the practical
value for an expensive pipeline.

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
fraction of cached runtime and speedup approaches $1/(1-f)$. It therefore
plateaus even though the absolute time returned continues to grow. The simple
runtime model closely described these measurements. Across all 12
cache-by-load observations, predicted and measured relative runtimes differed
by no more than 0.56 percentage points. This supports estimation
for a user's own pipeline from three quantities that can be measured locally:
direct runtime, cache hit rate and fixed VCFcache overhead. Such estimates do
not replace matched measurement when storage, representation or annotation
recipes differ.

## Cached outputs preserve annotator-specific results

Correctness comparison was performed after, and excluded from, every timed
cell. All 156 VEP cached outputs passed a streaming semantic comparator that
checked variant keys, selected input fields, genotype, the CSQ definition and
the complete set of CSQ entries. It canonicalized semantically irrelevant
ordering and reported but excluded only the documented VEP 115.2 HGNC_ID
batching discrepancy from pass/fail. Across the 156 external-WGS comparisons,
8,898 of 734,582,949 compared records (0.0012%) differed only in HGNC_ID; the
maximum was 495 records in one 4.84-million-variant output. In a separate
89-output technical validation set, the 265 reported annotation differences
were likewise confined to HGNC_ID, with no unexpected annotation mismatch.
All 156 fastVEP cached outputs passed a
stricter complete-record and relevant-header comparison, including
supplementary INFO and FORMAT fields, with only ordering canonicalized. All 12
outputs in the pipeline-complexity experiment also passed their specified
semantic gate.

The benchmarks have three important boundaries. First, the primary claim is
annotation equivalence under a fixed tool, version, reference and recipe; it
is not byte identity or equivalence between annotators. Second, the observed
speedups include normalization, lookup, annotation of misses, assembly,
compression and indexing, but exclude cache construction and the external
comparator. Third, cache construction is a
one-time cost whose value depends on reuse. For few samples or frequently
changing recipes, direct annotation may be simpler. For recurring pipelines,
the bundled caches avoid local construction, and cohort-derived caches become
progressively less expensive per sample as their build cost is amortized.
