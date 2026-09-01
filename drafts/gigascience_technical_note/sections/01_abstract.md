# Abstract

**Background:** Variant annotators repeatedly compute the same consequences
when alleles recur across samples. VCFcache stores the output of a fixed
annotation recipe and sends only cache misses to the original annotator.

**Findings:** In 12 independent GRCh38 genomes using a bundled gnomAD
allele-frequency (AF) ≥1% cache, median whole-genome sequencing speedups were
9.50-fold for Ensembl VEP `--everything`, 2.54-fold for clinical fastVEP and
10.78-fold for clinical SnpEff. Corresponding exome speedups were 2.29-, 0.93-
and 2.66-fold, showing that very short workflows may not recover fixed cache
costs. Across 52 held-out genomes, bundled AF ≥10%, bundled AF ≥1% and disjoint
three-genome cohort caches accelerated VEP 4.75-, 9.21- and 4.59-fold and
fastVEP 1.84-, 2.05- and 1.81-fold. Controlled VEP pipelines spanning
4.2–22.7 hours returned 3.3–20.5 hours per genome as annotation work increased.
Three-genome cohort caches broke even after two subsequent genomes; at 1,000
genomes, measured VEP rates corresponded to 39.1–86.7 compute days returned.
All 176 eligible comparisons in the final clinical campaign were semantically
concordant. Four SnpEff cells were excluded before inference because direct
SnpEff itself exhibited a processing-order-dependent warning.

**Conclusions:** Relative acceleration is governed by cache coverage and
recoverable per-variant work, whereas absolute time saved increases with
pipeline cost. VCFcache therefore benefits annotation-rich cohorts most while
retaining the configured annotator as the authority for every cache miss.
