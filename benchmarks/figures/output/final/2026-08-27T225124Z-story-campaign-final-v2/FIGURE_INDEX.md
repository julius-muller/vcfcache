# Story-driven final figure set

This is an additive candidate set. It does not replace or overwrite earlier snapshots.
All plots are vector-native R/ggplot2 outputs and are also exported as PDF and 600 dpi PNG.

## Manuscript figures

**Figure 1: VCFcache reuses annotation work without changing the completed record set.** (A) Each input VCF is matched against a recipe-specific cache. Existing annotations are reused, only cache misses are sent to the configured annotator, and both streams are merged and sorted into a complete output VCF. Users may download a ready-to-use recipe cache or a blueprint from Zenodo. (B) Benefit is governed by cache coverage, the cost of the underlying annotation recipe and the number of samples over which any local cache build is amortized.

**Figure 2: VCFcache returns most WGS annotation time across three native annotator configurations.** Twelve independent GRCh38 genomes (six KPGP and six SGDP) were tested as ACMG SF v3.3 Panel, Twist Human Core Exome WES and WGS inputs with the bundled gnomAD AF ≥1% strategy. VEP used `--everything`; fastVEP and SnpEff used their clinical recipes. Configurations are compared only against their own direct baseline. (A) Distribution of end-to-end runtime remaining after caching; the dashed line marks equal runtime. (B) Median serial throughput calculated as 86,400 seconds divided by median per-sample wall time; open circles denote direct annotation and arrowheads denote VCFcache. Cohort shapes identify independent biological samples. KPGP, Korean Personal Genome Project; SGDP, Simons Genome Diversity Project; VEP, Variant Effect Predictor; WES, whole-exome sequencing; WGS, whole-genome sequencing.

**Figure 3: Cache coverage defines the theoretical opportunity, while annotation cost determines how fully it is realized.** (A) Observed runtime remaining for the 108 matched GRCh38 sample-tool-assay pairs. Dashed curves show the zero-overhead expectation, `1 − f`, for hit rate `f`; point area represents direct annotation time. Fixed lookup, preprocessing and output costs dominate very short panel jobs but become small relative to WGS annotation. (B) One 4.80-million-variant KPGP genome was processed with six controlled VEP loads and bundled gnomAD AF ≥10% and AF ≥1% caches. Realized returnable time is observed time saved divided by `f ×` direct runtime. Point area denotes absolute hours returned. All 12 outputs passed semantic comparison.

**Figure 4: A cohort-derived cache repays its build cost rapidly and returns cumulative compute time.** Three held-out real-world cohorts were evaluated with caches built from three disjoint genomes. (A) Effective speedup includes the complete one-time cache-build duration divided across the indicated cohort size; circles mark empirical break-even. (B) Compute days returned after build cost for cohorts of 10, 100 and 1,000 samples, using the measured cohort-specific median direct and cached WGS runtimes. KPGP, Korean Personal Genome Project; PGP, Personal Genome Project; SGDP, Simons Genome Diversity Project.

## Supplementary figure

**Supplementary Figure 1: Annotation load increases reuse benefit without a cache-attributable change in output.** (A) Paired SnpEff speedups for core and clinical recipes on the same GRCh37 PGP samples. Lines connect identical samples; the dashed line marks equal runtime. The documented warning-only sample was excluded before inference for WES and WGS. (B) Of 180 direct/cached cells in the final clinical campaign, 176 were eligible and semantically concordant. Four SnpEff cells were excluded because the direct annotator itself changed output with input order; these are not counted as evidence for or against VCFcache. No eligible pair contained a cache-attributable mismatch.

## Repository figures

- `headliner_one_day_of_compute`: intuitive serial WGS throughput under each native benchmark recipe.
- `headliner_wgs_clock`: direct-to-cached median WGS wall time with speedup and time returned.
- `cohort_100_time_returned`: presentation-friendly cumulative compute-time view for 100 WGS samples.

## Interpretation boundary

The three annotators intentionally retain their native recipes. These panels estimate within-tool VCFcache benefit and must not be interpreted as rankings of annotation completeness or cross-tool semantic equivalence. No 1000 Genomes sample appears in Figures 2 or 3A. Throughput assumes serial continuous processing under the benchmark allocation and excludes cache download/build time; Figure 4 explicitly includes local cohort-cache build time.

