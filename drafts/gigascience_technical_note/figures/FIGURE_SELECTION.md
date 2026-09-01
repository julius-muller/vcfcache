# Figure selection

## Main figures

1. `main_v4/figure1_workflow_and_value_v4` connects ready Zenodo caches,
   Zenodo blueprints and cohort-derived blueprints to the per-sample cache-hit,
   miss-annotation and sorted-output paths.
2. `main_v4/figure2_empirical_time_and_throughput_v4` uses 12 matched GRCh38
   genomes to show runtime remaining and serial throughput for native VEP,
   clinical fastVEP and clinical SnpEff configurations across Panel, WES and
   WGS. It contains no 1000 Genomes-derived performance sample.
3. `main_v4/figure3_runtime_rule_v4` links empirical runtime remaining to cache
   coverage across 108 observations and uses controlled VEP loads to show that
   absolute time returned grows with pipeline cost.
4. `main_v4/figure4_cohort_economics_v4` includes measured cohort-cache build
   cost and reports effective speedup and compute days returned as cohort size
   increases.

## Supplementary figures

1. `supplementary_figure1_external_wgs_v2` reports cache coverage and wall
   time returned for all 52 independent evaluation genomes.
2. `supplementary_figure2_fastvep_cpu_complexity` tests whether process-wide
   CPU availability removes the benefit of reuse in core and deliberately
   annotation-dense fastVEP configurations on one held-out WGS.
3. `v4/supplementary_figure3_annotation_load_fidelity_v4` pairs core and
   clinical SnpEff results on the same genomes and accounts for every one of
   the 180 final-campaign direct/cached comparisons.

No figure contains an internal title, subtitle or provenance footer. Panel
letters remain in the artwork; numbered legends and abbreviation definitions
are placed immediately below figures in the manuscript.

Each composite is accompanied by separate 600-dpi PNG, vector PDF and SVG
exports for every panel under the corresponding `panels/` directory.

## Candidate headlines

The `candidates_v4/` directory contains three simplified alternatives from the
same frozen evidence snapshot: one-day throughput, WGS wall clocks and compute
time returned for 100 genomes. They are not cited by the Technical Note and
are retained for a bioRxiv landing page, repository README, graphical abstract
or later brief article format. They must not replace a manuscript panel without
updating the legend and evidence mapping.

## Provenance

`PROVENANCE.tsv` records every copied or rendered asset, its role, original
repository location, byte size and SHA-256 digest. The build check fails if a
listed asset changes without regenerating the provenance table.
