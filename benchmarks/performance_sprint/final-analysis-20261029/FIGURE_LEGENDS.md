# Figure legends

Final analysis, 29 October 2026. Preferred complete set: **`figures-v2/`**.
Earlier `figures/` exports are retained unchanged. Figures are vector PDF/SVG
and 600-dpi PNG. Panel letters identify subfigures; figure numbers, explanatory
captions and caveats belong below the graphics, not inside their plotting area.
These are performance/decision figures, not a replacement for the existing
workflow illustration.

## Figure 1: Assay-dependent acceleration and absolute WGS time savings

**A,** Per-sample speedup from VCFcache relative to the same sample's own direct
annotation baseline, using four frozen recipes and the adaptive engine. Dots
and triangles identify KPGP and SGDP genomes. Black diamonds indicate medians;
error bars are 95% percentile bootstrap intervals from 4,000 resamples within
cohort. The common GRCh38 subset has 12 donors, six per cohort, for each assay
and recipe. The dashed line denotes equal runtime; values below one indicate
slower execution with caching. **B,** Paired elapsed time saved per WGS sample.
Box centers are medians, boxes span the interquartile range and whiskers extend
to observations within 1.5 interquartile ranges; all donor observations remain
visible. Each recipe uses the gnomAD AF ≥1% blueprint with its compatible
annotation payload. Panel and WES inputs are interval-restricted WGS callsets,
not independent assay sequencing cohorts. Recipes differ in information content;
this is not an absolute speed ranking of annotators. VEP direct comparisons use
the predeclared HGNC_ID semantic handling and DOMAINS-order canonicalization;
other displayed recipes passed strict
record/header comparison. Clinical vep-rs timing-only results are excluded.

Abbreviations: WGS, whole-genome sequencing; WES, whole-exome sequencing;
VEP, Variant Effect Predictor; AF, allele frequency; gnomAD, Genome Aggregation
Database; KPGP, Korean Personal Genome Project; SGDP, Simons Genome Diversity
Project. Eight CPUs were allocated/bound per timed condition. Cache preparation
is excluded from these per-sample timing ratios.

## Figure 2: Population and three-genome cache strategies across 52 WGS genomes

**A,** Observed variant hit fractions for gnomAD AF ≥10%, gnomAD AF ≥1%, and
three-genome cohort caches, stratified by KPGP (20 genomes), SGDP (20) and PGP
(12). Each point represents a test genome; hit fractions are shown once rather
than duplicated for both annotators. Boxes and whiskers follow Figure 1.
**B,** Empirical cumulative distributions of per-sample speedup, with VEP
everything and basic fastVEP shown separately. Curves farther right indicate
larger acceleration. Dashed vertical lines denote equal runtime. Within each
recipe/donor the same direct baseline is used for all three strategies. Twelve
shared VEP measurements are represented as plotting aliases, not extra timings
or donors. The clinical fastVEP recipe from Figure 1 is not used in this family.
PGP uses GRCh37; KPGP and SGDP use GRCh38. Cross-cohort differences are therefore
not attributed solely to ancestry or relatedness.

Each cohort cache was built from three different seed genomes whose IDs do not
occur among that cohort's test genomes. “Held-out” refers only to these seeds;
absence from gnomAD's donor set has not been established. VEP population caches
are bundled Zenodo assets; fastVEP payloads were locally annotated from the
corresponding blueprints. PGP, Personal Genome Project. Other abbreviations and
comparison scopes follow Figure 1.

## Figure 3: Measured pipeline phases and conditional optimization headroom

**A,** Mean per-sample fractions of cached elapsed time spent annotating misses,
routing/lookup, preprocessing, output assembly, other engine work, and CLI/other
work in the common GRCh38 subset. Fractions sum to one per condition before
averaging; durations of parallel shards are not added as if sequential.
“Annotate misses” includes external recipe invocation, startup and I/O.
**B,** Hypothetical median direct/cached ratios when increasing fractions of
the measured non-annotation time are removed, holding each sample's measured
miss-annotation time fixed. The horizontal dashed line marks parity. These
curves are optimistic counterfactual calculations, **not additional benchmark
measurements, a guaranteed optimization path, or a model of changing hit rate**.
They identify where cache-side engineering could matter and where remaining
annotation dominates. CLI, command-line interface. Recipes, donors, comparisons
and allocation follow Figure 1.

## Figure S1: Pipeline complexity increases the absolute value of reuse

**A,** Direct and cached elapsed runtimes for the real WGS genome KPGP-00319,
using a configurable self-written VEP wait plugin and two gnomAD cache
strategies. The observed hit rates are 80.23% and 90.26%. Wait duration is
specified **per consequence**, not per input variant. **B,** Observed speedups
at each wait setting; lines connect descriptive observations and are not fitted
scaling laws. Each setting has one timing; no biological or technical uncertainty
interval is implied. The same historical direct measurement is shared across
the two cache strategies at each wait setting. Direct timings are reused, and
cached outputs were checked using the composed vanilla-VEP semantic reference;
retained full direct BCFs for each delay are not available. Waiting simulates
additional computational burden, not the information content or threading
behavior of real plugins. Eight allocated CPUs; cache construction excluded.

## Figure S2: Accepted assay results separated by genome assembly

Per-sample speedups for every accepted assay recipe, shown separately for
GRCh37 and GRCh38 and for Panel/WES/WGS callsets. Point color and shape identify
cohort. The vertical dashed line marks parity. GRCh38 has six KPGP and six SGDP
donors per available recipe/assay; GRCh37 has 12 PGP donors except SnpEff core
and clinical WES/WGS, where one donor is excluded after strict comparison
failure (n=11). SnpEff core was evaluated only for GRCh37; matched VEP and
vep-rs extensions only for GRCh38. Clinical vep-rs contributes accepted Panel
points only; its WES/WGS ratios are not validated evidence. Missing combinations
are not zero speedups. Full comparator outcomes, including excluded records,
are documented in the text and audit tables rather than plotted as performance
claims. Assay derivation, recipes and abbreviations follow Figures 1 and 2.

## Figure S3: Aggregate CPU cost and peak memory in fresh paired controls

Cached/direct ratios of aggregate CPU time and peak cgroup memory, separated
by assembly and assay. Boxes span the interquartile range, center lines are
medians and whiskers extend to observations within 1.5 interquartile ranges.
Only accepted conditions with newly measured paired direct resource controls
are included; historical direct resources are not reconstructed. Each timed
condition ran in its own Slurm step/cgroup, avoiding contamination of one
condition's memory peak by a previous condition. CPU time is the aggregate
cgroup counter difference. Memory includes charged file-backed/page-cache
memory, not only application heap or one child's resident set. Ratios below
one favor caching. Lower elapsed time can coexist with increased peak memory;
neither resource nor elapsed ratios prove cluster-wide throughput scalability.

## Figure S4: Historical cached runtime versus the current adaptive campaign

Distributions of historical/current cached-runtime ratios for accepted WGS
conditions, separated by cache strategy and colored by cohort. Ratios above
one indicate lower cached elapsed time in the current run. Boxes/whiskers follow
Figure S3. Compatible frozen inputs, recipes and historical baseline bindings
are retained, but dates, nodes, operating-system cache state and some storage
placements differ. **This is not a fresh randomized engine A/B trial**, and the
ratios cannot be attributed exclusively to one code change. Values from this
figure are not substituted for each recipe's direct-annotation speedup.

## Figure S5: Illustrative amortization under explicitly assumed build costs

Effective per-sample cost relative to direct annotation as the number of samples
increases. Three illustrative build costs—zero, ten and one hundred direct-sample
runtimes—are **assumptions**, not new build measurements. The steady-state
cached fraction is derived from the median observed WGS speedup for each
displayed recipe in Figure 1, and the assumed one-time build cost is divided
by sample count. The horizontal dashed line marks parity. Zero annotation-build
cost does not imply that cache download or other setup takes zero time. This
figure illustrates the runtime model rather than predicting a particular
laboratory's setup cost or cohort throughput.

## Figure S6: Conditional fastVEP amortization using payload-bound build provenance

Conditional cohort-level speedup combining current basic-fastVEP WGS **mean**
direct/cached sample runtimes with the recorded historical annotation-build
duration matching each tested cache payload's SHA256. Curves are separated by
cohort and cache strategy. The recorded builds start from available blueprints;
download and blueprint preparation are not included, and historical build
hardware/settings were not rematched to the current campaign. Thus, these are
conditional amortization scenarios, not new end-to-end preparation or cohort
benchmarks. The x-axis counts new genomes processed after preparation, excluding
the three seed genomes where applicable. The dashed line marks parity.

## Diagnostic D1: Clinical vep-rs timings retained outside validated results

Core and clinical per-sample timing ratios in the pinned vep-rs extension.
Filled points passed complete-record/header comparison; hollow points retain
ClinVar-derived CSQ discrepancies despite matching relevant headers and
non-CSQ content. Hollow points are **timing-only observations**, not validated
speedups. Their diagnostic separation prevents the promising clinical ratios
from being mistaken for output-equivalent manuscript evidence. CSQ, consequence
annotation INFO field. This graphic is for decision review, not a primary
performance figure.

## Repository graphic: One-glance WGS time reduction

Rounded gray bars represent each recipe's own direct runtime; green segments
represent its cached fraction derived from the median WGS speedup in Figure 1.
Displayed percentages summarize that same ratio. These are distinct recipes,
not an absolute runtime comparison between tools. The n=12 GRCh38 observations,
AF ≥1% caches, eight-CPU settings and comparator limitations follow Figure 1.
The graphic deliberately omits setup amortization and does not imply the same
benefit at WES, at every hit rate, or for every future sample.
