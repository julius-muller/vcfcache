# Complete VCFcache figure set

Rendered: 2026-08-07

This is the first complete figure directory containing the finalized VEP,
fastVEP, cross-annotator, statistics-calibration, and repository-headline
outputs. It does not replace or modify any earlier figure directory.

Every composite and component panel was produced in R with ggplot2 and is
available as SVG, PDF, and 300 dpi PNG. Exact derived values are under
`source/`; renderer paths and R session details are recorded in the four
`RENDERED_FROM*` and `R_SESSION_INFO*` files.

## Recommended manuscript figures

1. `manuscript/external_wgs.*`
   Primary real-world VEP evidence across 52 KPGP, SGDP, and PGP evaluation
   genomes, three cache strategies, and custom-cache amortization.
2. `manuscript/fastvep_real_world_wgs.*`
   The matching complete fastVEP campaign. All 156 cached outputs passed strict
   complete-record and relevant-header comparison.
3. `manuscript/vep_fastvep_impact_final.*`
   Exactly paired genomes, variants, and cache blueprints across VEP and
   fastVEP. Panel D reports the completed six-genome VEP light/full statistics
   sensitivity analysis; no original VEP measurement was rescaled.
4. `manuscript/assay_source_overlap_calibration.*`
   Panel, WES, and WGS runtime mechanics using 1000 Genomes-derived inputs.
   These data are source-overlap calibration and must not be presented as
   independent public-cache coverage estimates.

Each manuscript composite has separately exported panels under
`manuscript/panels/<figure-name>/`.

## Optional manuscript views

- `alternatives/external_wgs_decision.*`: observed wait, reuse, and hours
  returned per evaluation genome.
- `alternatives/runtime_model_diagnostics.*`: miss-fraction model, inferred
  overhead, and speedup efficiency.

## Repository headline alternatives

- `repository/headliner_before_after.*` — recommended quantitative GitHub
  headline. It shows the observed gnomAD AF >= 1% result for both annotators:
  8.6x median VEP speedup and 2.0x median fastVEP speedup across 52 genomes.
- `repository/headliner_workflow_doodle.*` — simplest package explainer. It
  combines bundled/cohort caches, VEP/fastVEP support, realistic reuse, observed
  speedup range, and validated output equivalence in one workflow-style image.
- `repository/headliner_scenario_map.*` — reader self-selection by cache hit
  rate and pipeline overhead. Values are explicitly modeled from median
  empirical overhead rather than measurements of those exact scenarios.
- `repository/user_impact_model_preview.*` — retained prior model graphic for
  comparison; the three new headline alternatives are preferable for the
  repository landing page.

## Statistical and validation notes

- Samples, not technical reruns, are the inferential units.
- Cross-cohort intervals use cohort-stratified bootstrap resampling.
- Historical VEP results used `--statistics full`; fastVEP used
  `--statistics light`.
- In six paired VEP calibration genomes, light/full median speedup ratios were
  +4.1% for gnomAD AF >= 10%, +7.2% for gnomAD AF >= 1%, and +3.2% for the
  three-genome cohort cache. This establishes the historical VEP results as
  slightly conservative without transforming the original observations.
- All 52 fastVEP tasks and all six VEP calibration tasks completed. All 156
  fastVEP cached outputs and all 18 calibration cached outputs passed their
  prespecified semantic comparators.

## Frozen source snapshots

- Original complete VEP/assay snapshot:
  `source_data/final/2026-08-06T095303Z`
- Final annotator and statistics-calibration snapshot:
  `source_data/final/2026-08-07T053247Z-complete`

