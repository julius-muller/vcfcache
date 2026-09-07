# Findings

## VCFcache reuses complete, recipe-specific annotations

Across separately annotated samples, the same alleles recur, yet conventional
annotation recomputes the requested annotations on every run. VEP and SnpEff
use local reference data to avoid repeated remote database access, but still
evaluate each input allele [@mclaren2016vep; @cingolani2012snpeff].
Lookup-based tools such as vcfanno and echtvar efficiently retrieve
pre-existing annotation fields from reference files or indexed archives, but
do not themselves invoke a user-defined annotation workflow for variants
absent from those resources [@vcfanno2016; @pedersen2023echtvar].
VCFcache combines reuse with automatic fallback: it retrieves complete
per-allele output produced by a fixed annotation recipe—defined by the
annotator version, command-line options, plugins and supplementary
resources—and sends cache misses through that same recipe. Cached and newly
annotated records are then merged into a sorted, indexed BCF (Figure 1).

![](figures/main_v4/figure1_workflow_and_value_v4.png){width=100%}

**Figure 1: VCFcache reuses annotation work without changing the completed
record set.** (A) Each input VCF is matched against a recipe-specific cache.
Existing annotations are retrieved, only cache misses are sent to the configured
annotator, and both streams are merged into a sorted, indexed output BCF.
Users may download a ready-to-use recipe cache or a blueprint from Zenodo. (B)
Benefit is governed by cache coverage, the cost of the underlying annotation
recipe and the number of samples over which any local cache build is amortized.
BCF, binary variant call format; VCF, variant call format.

For $N$ variants, annotator start-up $s_{\mathrm{ann}}$, per-variant cost
$t_{\mathrm{ann}}$, cache hit rate $f$ and VCFcache processing cost $C(N)$,

$$
T_{\mathrm{direct}}=s_{\mathrm{ann}}+Nt_{\mathrm{ann}}, \qquad
T_{\mathrm{cached}}\approx s_{\mathrm{ann}}+(1-f)Nt_{\mathrm{ann}}+C(N).
$$

Caching helps when $fNt_{\mathrm{ann}}>C(N)$. It can return per-variant work,
but not annotator start-up; high coverage therefore need not accelerate very
small inputs or minimal recipes. A custom cache's one-time build cost is
amortized as
$T_{\mathrm{effective}}\approx T_{\mathrm{cached}}+T_{\mathrm{build}}/S$
after $S$ uses.

## Assay scale and annotation cost determine practical benefit

Twelve independent GRCh38 genomes (six KPGP and six SGDP) were represented as
matched ACMG panel, WES and WGS inputs and evaluated with the bundled gnomAD
AF ≥1% cache. Each native annotator configuration was compared only with its
own direct baseline. VEP `--everything` produced median speedups of 0.71-,
2.29- and 9.50-fold for Panel, WES and WGS, respectively. The corresponding
values were 0.46-, 0.93- and 2.54-fold for clinical fastVEP and 1.06-, 2.66-
and 10.78-fold for clinical SnpEff (Figure 2). Thus fixed lookup and output
costs can outweigh saved work for small inputs, whereas WGS contained enough
recoverable work to accelerate all three configurations.

![](figures/main_v4/figure2_empirical_time_and_throughput_v4.png){width=100%}

**Figure 2: VCFcache returns most WGS annotation time across three native
annotator configurations.** Twelve independent GRCh38 genomes were tested as
ACMG SF v3.3 Panel, Twist Human Core Exome WES and WGS inputs with the bundled
gnomAD AF ≥1% strategy. VEP used `--everything`; fastVEP and SnpEff used their
clinical recipes. Configurations are compared only against their own direct
baseline. (A) Distribution of end-to-end runtime remaining after caching; the
dashed line marks equal runtime. (B) Median serial throughput calculated as
86,400 seconds divided by median per-sample wall time; open circles denote
direct annotation and arrowheads denote VCFcache. Cohort shapes identify
independent biological samples. AF, allele frequency; KPGP, Korean Personal
Genome Project; SGDP, Simons Genome Diversity Project; VEP, Variant Effect
Predictor; WES, whole-exome sequencing; WGS, whole-genome sequencing.

For WGS, median serial throughput increased from 9.3 to 88.4 genomes per day
for VEP, from 186.3 to 471.2 for fastVEP and from 9.9 to 105.4 for SnpEff under
the benchmark allocations. These values show the operational effect of reuse
within each recipe; they are not a ranking of annotators because the native
recipes intentionally differ in annotation scope.

## Cache coverage defines opportunity; pipeline cost determines time returned

Across the 108 matched GRCh38 sample-tool-assay pairs, runtime remaining
approached the zero-overhead expectation, $1-f$, as direct annotation became
more expensive (Figure 3A). The deviation for Panel and some WES cells shows
the fixed cost represented by $C(N)$; hit rate alone is insufficient to predict
benefit for short jobs.

![](figures/main_v4/figure3_runtime_rule_v4.png){width=100%}

**Figure 3: Cache coverage defines the theoretical opportunity, while
annotation cost determines how much time is returned.** (A) Observed runtime
remaining for 108 matched GRCh38 sample-tool-assay pairs. Dashed curves show
the zero-overhead expectation, $1-f$, for hit rate $f$; point area represents
direct annotation time. Fixed lookup, preprocessing and output costs dominate
very short Panel jobs but become small relative to WGS annotation. (B) One
4.80-million-variant KPGP genome was processed with six controlled VEP loads
and bundled gnomAD AF ≥10% and AF ≥1% caches. Solid lines show observed hours
returned as pipeline cost increases; dashed lines show the cache-coverage
opportunity, $f \times$ direct runtime. All 12 outputs passed semantic
comparison. AF, allele frequency; KPGP, Korean Personal Genome Project; VEP,
Variant Effect Predictor; WGS, whole-genome sequencing.

The controlled VEP recipes spanned 4.2--22.7 hours without changing their
annotation content. At 80.23% cache coverage, 19.7--21.7% of direct wall time
remained (4.62--5.07-fold speedup); at 90.26% coverage, 9.9--11.9% remained
(8.39--10.10-fold). For the longest recipe, VCFcache returned 18.19 and 20.45
hours per genome at the two cache depths. Measured and coverage-predicted
relative runtimes differed by at most 0.56 percentage points. Relative speedup
therefore approaches a hit-rate ceiling even while the absolute time returned
continues to grow with pipeline cost.

## Public and cohort caches support different deployment models

We next evaluated 20 KPGP, 20 SGDP and 12 Harvard PGP genomes
[@kpgp_data; @mallick2016sgdp; @ball2012pgp]. KPGP and SGDP remained in GRCh38
and PGP in GRCh37. Three additional, disjoint genomes from each provider formed
cohort caches; relatedness screening prevented a second-degree-or-closer
training/evaluation pair. The public alternatives were bundled,
assembly-matched gnomAD v4.1 caches containing alleles reaching AF ≥10% or
AF ≥1% in any reported frequency stratum [@gnomad2024release].

Across these 52 held-out genomes, VEP median end-to-end speedups were 4.75-fold
(AF ≥10%), 9.21-fold (AF ≥1%) and 4.59-fold (cohort cache). With the core
fastVEP recipe, they were 1.84-, 2.05- and 1.81-fold. Median AF ≥1% coverage
was 91.95--93.04% across providers, whereas three-genome cohort coverage was
77.45--82.88% (Supplementary Figure 1).

![](figures/main_v4/figure4_cohort_economics_v4.png){width=100%}

**Figure 4: A cohort-derived cache repays its build cost rapidly and returns
cumulative compute time.** Three held-out real-world cohorts were evaluated
with caches built from three disjoint genomes. (A) Effective speedup includes
the complete one-time cache-build duration divided across the indicated cohort
size; circles mark empirical break-even. (B) Compute days returned after build
cost for cohorts of 10, 100 and 1,000 samples, using measured cohort-specific
median direct and cached VEP WGS runtimes. KPGP, Korean Personal Genome Project;
PGP, Personal Genome Project; SGDP, Simons Genome Diversity Project; VEP,
Variant Effect Predictor; WGS, whole-genome sequencing.

All three cohort caches repaid their measured build cost after two subsequent
genomes; downloaded bundled caches require no local annotation build. At 100
genomes, the observed VEP rates corresponded to 8.5, 8.3 and 3.8 compute days
returned for KPGP, SGDP and PGP, respectively. At 1,000 genomes, these values
were 86.7, 84.7 and 39.1 days. The calculation includes the complete local
cache-build time and therefore distinguishes immediate public-cache use from
cohort-specific amortization.

## Caching complements efficient and parallel annotators

Increasing annotation load strengthened benefit without changing cache
membership. On the same GRCh37 genomes, extending SnpEff from core consequence
annotation to a clinical ClinVar and gnomAD recipe increased median speedup
from 1.16- to 3.62-fold for WES and from 4.46- to 12.31-fold for WGS
(Supplementary Figure 3 and Supplementary Table 1).

An independent engineering experiment asked whether more CPUs eliminate reuse
benefit for fast native annotation. A dense ten-database fastVEP recipe
saturated between process-wide limits of 10 and 32 CPUs (706.2 versus 691.4 s
direct) but retained 3.07- and 3.37-fold speedups at approximately 90% cache
hits. At 10 CPUs, caching returned 7.9 minutes per genome (Supplementary Figure
2). Parallelism and reuse therefore reduce different components of runtime.

## Correctness was assessed outside timed cells

All 156 VEP and 156 fastVEP whole-genome pairs passed annotator-specific
semantic comparison. For VEP, 8,898 of 734,582,949 records (0.0012%) differed
only in `HGNC_ID`, a documented VEP 115.2 batching effect [@vep_issue_1959];
the field was reported but excluded from pass/fail. fastVEP comparison covered
all INFO, FORMAT and relevant header values, canonicalizing order only.

In the final clinical campaign, all 176 eligible direct/cached pairs were
semantically concordant (Supplementary Figure 3). Four of 180 SnpEff pairs
were excluded before inference because one GRCh37 allele differed only by
`WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS`; consequence, transcript, HGVS,
ClinVar, gnomAD, genotype and record counts agreed. A direct-SnpEff reproducer
showed that the warning changed when two adjacent input records were reversed,
without VCFcache involvement. We therefore documented the annotator-order
effect rather than adding warning-specific normalization. These results
support equivalence under a fixed tool, version, reference and recipe, not byte
identity or equivalence between annotators.
