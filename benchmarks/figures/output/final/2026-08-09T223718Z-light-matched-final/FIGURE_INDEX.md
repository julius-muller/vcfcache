# Light-mode-matched final VCFcache figure set

Rendered: 2026-08-10

This is the recommended complete figure directory. VEP 115.2 and fastVEP 0.3.0
are both represented by `--statistics light` measurements across the same 52
real-world WGS genomes, variants, cache blueprints, and cache strategies. No
timing value in the primary comparison is calibrated, rescaled, or imputed.

Earlier snapshots and figure directories remain unchanged. Every composite and
component panel was produced in R with ggplot2 and is available as SVG, PDF,
and 300 dpi PNG.

## Recommended manuscript figures

1. `manuscript/external_wgs.*`
   Primary VEP evidence using the complete 52-genome light-statistics rerun.
   Median speedups across all genomes are 4.75x for gnomAD AF >= 10%, 9.21x
   for gnomAD AF >= 1%, and 4.59x for the three-genome cohort cache.
2. `manuscript/fastvep_real_world_wgs.*`
   Matching fastVEP evidence. Median speedups are 1.84x, 2.05x, and 1.81x for
   the same three strategies. All 156 cached outputs passed strict
   complete-record and relevant-header comparison.
3. `manuscript/vep_fastvep_light_matched_final.*`
   Direct annotator comparison with identical light-statistics timing mode.
   Panel D pairs the new VEP measurements with the historical full-rescan
   observations as a sensitivity analysis rather than altering either dataset.
4. `manuscript/assay_source_overlap_calibration.*`
   Panel, WES, and WGS runtime mechanics from 1000 Genomes-derived inputs.
   These remain source-overlap calibration data and must not be described as
   independent public-cache coverage estimates.

Every manuscript composite has separately exported panels under
`manuscript/panels/<figure-name>/`.

## Statistics-mode sensitivity

Across all 52 genomes, the median within-sample light/full speedup ratios were:

- gnomAD AF >= 10%: +3.6% (95% stratified bootstrap interval +2.8% to +4.5%)
- gnomAD AF >= 1%: +5.6% (+5.2% to +6.6%)
- three-genome cohort cache: +2.7% (+2.2% to +3.5%)

This confirms that the legacy full-output rescan made the historical VEP
speedups slightly conservative. The new light-mode observations are used
directly in all recommended primary figures.

## Repository headline alternatives

- `repository/headliner_before_after.*` — recommended measured headline:
  9.2x median VEP and 2.0x median fastVEP speedup with the bundled gnomAD
  AF >= 1% cache across 52 genomes.
- `repository/headliner_workflow_doodle.*` — simplest conceptual README image,
  combining bundled/cohort caches, VEP/fastVEP, realistic reuse, observed
  speedups, and semantic validation.
- `repository/headliner_scenario_map.*` — modeled user self-selection by hit
  rate and empirically inferred annotator overhead. Cells are estimates rather
  than measurements of those exact scenarios.
- `repository/user_impact_model_preview.*` — retained older model view; one of
  the three newer headlines is preferable for the repository landing page.

## Optional manuscript views

- `alternatives/external_wgs_decision.*`: wait, observed reuse, and hours
  returned per VEP-annotated genome using light-mode measurements.
- `alternatives/runtime_model_diagnostics.*`: miss-fraction behavior, inferred
  overhead, and speedup efficiency.

## Frozen source snapshots

- Complete VEP/assay source with light-mode external WGS:
  `source_data/final/2026-08-09T223718Z-light-complete-vep`
- Matched VEP/fastVEP source and full-versus-light sensitivity table:
  `source_data/final/2026-08-09T223718Z-light-matched-final`
- Historical VEP full-rescan source retained for sensitivity analysis:
  `source_data/final/2026-08-06T095303Z`

All 52 VEP light tasks, 52 fastVEP tasks, and their 312 cached outputs passed
the prespecified annotator-specific semantic comparators.

