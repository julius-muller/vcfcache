# Figure selection

## Main figures

1. `figure1_workflow_v3` connects ready Zenodo caches, Zenodo blueprints and
   cohort-derived blueprints to the per-sample split, annotation-reuse,
   miss-annotation and merge paths.
2. `figure2_assay_annotator_v3` uses matched independent GRCh38 genomes to show
   where Panel, WES and WGS cross the fixed-overhead break-even point for VEP
   and fastVEP. It contains no 1000 Genomes-derived performance sample.
3. `figure3_pipeline_complexity_v3` uses time-return arrows and explicit
   hit-rate ceilings to reconcile increasing absolute savings with bounded
   relative speedup.

## Supplementary figures

1. `supplementary_figure1_external_wgs_v2` reports cache coverage and wall
   time returned for all 52 independent evaluation genomes.
2. `supplementary_figure2_fastvep_cpu_complexity` tests whether process-wide
   CPU availability removes the benefit of reuse in core and deliberately
   annotation-dense fastVEP configurations on one held-out WGS.
3. `supplementary_figure3_matched_assays_v3` shows every selected genome's
   cache coverage and speedup across the matched assay extension.

No figure contains an internal title, subtitle or provenance footer. Panel
letters remain in the artwork; numbered legends and abbreviation definitions
are placed immediately below figures in the manuscript.

Each composite is accompanied by separate 600-dpi PNG, vector PDF and SVG
exports for every panel under the corresponding `panels/` directory.

## Candidate headlines

The `candidates/` directory contains three existing simplified alternatives.
They are not cited by the Technical Note and are retained for a bioRxiv landing
page, repository README, graphical abstract, or a later brief article format.
They must not replace the measured manuscript panels without updating the
legend and evidence mapping.

## Provenance

`PROVENANCE.tsv` records every copied or rendered asset, its role, original
repository location, byte size and SHA-256 digest. The build check fails if a
listed asset changes without regenerating the provenance table.
