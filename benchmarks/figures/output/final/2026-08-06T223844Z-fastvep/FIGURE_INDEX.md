# fastVEP and paired annotator figures

Source snapshot: `source_data/final/2026-08-06T223844Z-fastvep`

This versioned directory does not replace any earlier figure set. All plots
were produced in R with ggplot2. Every assembled figure and every component
panel is available as SVG, PDF, and 300 dpi PNG.

## Final fastVEP evidence

- `manuscript/fastvep_real_world_wgs.*` is the recommended standalone fastVEP
  figure. It contains all 52 KPGP, SGDP, and PGP evaluation genomes and all
  three cache strategies.
- `manuscript/panels/fastvep_real_world_wgs/` contains the separate speedup,
  hit-rate/overhead, and absolute-wall-time panels.
- All 156 cached fastVEP outputs passed strict complete-record and
  relevant-header comparison against the corresponding direct output.
- fastVEP used VCFcache `--statistics light`; each sample/condition was run
  once, so genomes rather than technical replicates are the inferential units.

Across all 52 genomes, the stratified median speedups are:

- gnomAD AF >= 10%: 1.84x (95% bootstrap interval 1.79-1.90x)
- gnomAD AF >= 1%: 2.05x (1.99-2.16x)
- three-genome cohort cache: 1.81x (1.79-1.85x)

The cohort-specific median fastVEP speedups span 1.46-2.16x. Median direct
wall time is about 3.5 minutes for PGP and 6.7-7.1 minutes for SGDP/KPGP; the
cache typically returns roughly 1.3-3.8 minutes per genome depending on cohort
and strategy.

## Provisional cross-annotator view

- `manuscript/vep_fastvep_impact_provisional.*` pairs the exact same genomes,
  variants, and cache blueprints between VEP and fastVEP.
- `manuscript/panels/vep_fastvep_impact_provisional/` contains its individual
  speedup, hit-rate/runtime, and absolute-wall-time panels.

This cross-annotator figure is deliberately labelled **provisional** because
the completed 52-genome VEP campaign used `--statistics full`, whereas fastVEP
used `--statistics light`. The standalone fastVEP evidence is unaffected. The
paired figure should be refreshed after the six-genome VEP light-statistics
calibration is complete.

Exact cohort medians and stratified bootstrap intervals are under `source/`.

