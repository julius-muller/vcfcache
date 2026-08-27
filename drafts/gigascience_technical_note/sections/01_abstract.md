# Abstract

**Background:** Variant annotators repeatedly compute the same consequences
when alleles recur across samples. VCFcache stores the output of a fixed
annotation recipe and sends only cache misses to the original annotator.

**Findings:** In 52 held-out whole genomes, VCFcache accelerated Ensembl VEP
4.75-, 9.21- and 4.59-fold with bundled gnomAD allele-frequency (AF) ≥10%,
gnomAD AF ≥1% and disjoint three-genome cohort caches, respectively. The same
strategies accelerated fastVEP 1.84-, 2.05- and 1.81-fold. With clinical
SnpEff recipes, median speedups were 2.66- and
3.62-fold for exomes and 10.78- and 12.31-fold for genomes on GRCh38 and
GRCh37, respectively; a lower-cost core SnpEff recipe reached 1.16- and
4.46-fold on GRCh37. Controlled VEP pipelines spanning 4.2–22.7 hours were
reduced to 19.7–21.7% of direct runtime at 80.23% cache hits and 9.9–11.9% at
90.26% hits. All 312 VEP/fastVEP whole-genome outputs passed semantic
comparison. Of 180 SnpEff/fastVEP assay pairs, 176 passed; four excluded pairs
differed only by one processing-order-dependent upstream SnpEff warning.

**Conclusions:** Relative acceleration is governed by cache coverage and
recoverable per-variant work, whereas absolute time saved increases with
pipeline cost. VCFcache therefore provides the largest benefit for genomes and
annotation-rich recipes while preserving the configured annotator as the
authority for cache misses.
