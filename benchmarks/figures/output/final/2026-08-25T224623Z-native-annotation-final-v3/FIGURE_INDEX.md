# Final native-annotation figure set

This snapshot is additive and does not overwrite any prior figure set.
All plots were rendered in R with ggplot2; composite panels use patchwork.
PNG files are 600 dpi, and every manuscript panel is also exported as SVG and PDF.

## Proposed manuscript figures

**Proposed Figure 2: Assay scale and recoverable work determine VCFcache benefit across native annotators.** Twelve held-out GRCh38 genomes (six KPGP and six SGDP) were evaluated as matched ACMG SF v3.3 Panel, Twist Human Core Exome WES and WGS inputs with the bundled gnomAD AF ≥1% strategy. VEP used `--everything`; fastVEP and SnpEff used their clinical recipes. Each recipe is compared only with its own direct baseline because annotation content differs across tools. (A) Median direct runtime (open circle), cached runtime (arrowhead) and median paired speedup. (B) Sample-level speedups; colour identifies annotator and shape identifies cohort. The dashed line marks equal runtime. (C) Median paired speedup and wall-time change. All 108 pairs passed their annotator-specific semantic gate. AF, allele frequency; KPGP, Korean Personal Genome Project; SGDP, Simons Genome Diversity Project; VEP, Variant Effect Predictor; WES, whole-exome sequencing; WGS, whole-genome sequencing.

**Proposed Supplementary Figure: Additional annotation work increases the benefit of reuse within SnpEff.** The same 12 GRCh37 PGP genomes and the same gnomAD AF ≥1% cache membership were processed with core SnpEff or a clinical recipe adding ClinVar and gnomAD exome lookups. (A) Sample-level paired change in speedup. Lines connect the same sample across recipes; the dashed line marks equal runtime. The warning-only pair described in the manuscript was excluded before inference, leaving 11 paired genomes for WES and WGS. (B) Median direct runtime (open circle), cached runtime (arrowhead) and median paired speedup. Annotation load, rather than cache coverage, moved median WES speedup from 1.16- to 3.62-fold and WGS speedup from 4.46- to 12.31-fold. PGP, Personal Genome Project; WES, whole-exome sequencing; WGS, whole-genome sequencing.

**Proposed Figure 3: Absolute time returned grows while relative speedup approaches the cache-coverage ceiling.** One held-out 4.80-million-variant KPGP genome was processed with six controlled VEP `--everything` loads and bundled gnomAD AF ≥10% and AF ≥1% caches. (A) Direct runtime (open circle), cached runtime (arrowhead), and hours returned per genome. (B) Speedup as direct annotation becomes slower. Dashed lines are the theoretical `1/(1-f)` ceilings for observed hit rate `f`. The synthetic plugin emitted no annotation field, and all 12 cached outputs passed semantic comparison. AF, allele frequency; VEP, Variant Effect Predictor; WGS, whole-genome sequencing.

## Repository and alternative figures

- `headliner_native_wgs`: deliberately simple normalized-time view for the repository front page; bars compare each native recipe only with its own direct baseline.
- `native_wgs_clock_lines`: compact direct-to-cached runtime alternative; open circles are direct annotation and filled circles are VCFcache.

## Interpretation boundary

The combined annotator panels show within-tool speedups under native recipes. They must not be used to rank annotation completeness or claim cross-tool semantic equivalence. The 12 genomes are independent KPGP/SGDP evaluation samples and are not 1000 Genomes samples.

