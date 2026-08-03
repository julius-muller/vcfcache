# Publication figure blueprint and accelerated execution design

## Design decision

The biological sample is the independent unit. The default benchmark therefore
uses one validated paired or multi-condition execution per sample, with balanced
condition order. Whole-cohort warm-ups and three technical reruns were removed
before publication collection. The existing HG02079 repeat differed by 0.84%
uncached and 0.08% cached, while the 49 completed primary WGS samples show much
larger between-sample variation.

A smoke task validates staging, memory, commands, and semantic comparison. It is
not a second full-cohort phase. Any small repeatability control is shown
separately and never counted as an independent sample.

## Main figure

### A. What the cache avoids

For one representative sample, show total input variants split into cache hits
and misses, with only misses flowing into VEP. Put measured wall time under each
stage. This is the intuitive explanation of the mechanism.

### B. Real assays

Plot paired relative runtime (`uncached = 1`) for at least 12 panel, 12 WES, and
12 WGS samples using the bundled Zenodo gnomAD AF >= 1% VEP 115.2 `--everything`
cache. Show individual samples plus median and IQR. The already completed 49
primary WGS first-pass pairs are valid source measurements; the one cgroup-OOM
sample is reported as an infrastructure exclusion, not silently replaced.

### C. Real independent WGS cache strategies

For all 52 held-out KPGP, SGDP, and PGP genomes, show hit rate and relative
runtime for the bundled AF >= 10%, bundled AF >= 1%, and disjoint three-genome
cohort cache. Each sample has one common uncached baseline and all three cached
conditions. Facet or color by cohort and assembly rather than pooling GRCh37 and
GRCh38 as if they were one source population.

### D. Runtime scaling mechanism

On one representative sample, compare:

1. vanilla offline VEP;
2. VEP `--everything`;
3. vanilla VEP plus `SyntheticDelay,delay_us=N` at two calibrated delays.

For each scenario, measure controlled cache hit rates and overlay
`T_cached = T_overhead + (1-f) * T_uncached`. The synthetic plugin sleeps per
transcript consequence, emits no field, and is invoked only for cache misses.
This deliberately varies miss cost while holding lookup and output semantics
constant. Calibrate delay values in a short pilot to create clearly separated,
publication-readable baseline runtimes rather than choosing arbitrary long
delays.

## Cohort economics companion

Use the empirical external-WGS medians and measured custom-cache build cost in
`T_eff = T_cached + T_build/S`. Plot every integer `S=1..1000`, marking 5, 10,
and 100 samples. This is a model based on measured components, not a claim that
100 technical reruns were performed.

## Infrastructure decision

Do not resize workers in the middle of the primary/external comparison. The
current six identical nodes preserve comparability, and removing redundant
passes reduces the critical path more than a risky flavor migration. A future
controlled-mechanism campaign may use resized workers only if every scenario is
run on that same hardware and the hardware change is reported as a separate
campaign.
