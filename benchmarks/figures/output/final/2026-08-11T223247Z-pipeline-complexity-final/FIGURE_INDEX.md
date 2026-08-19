# Final real-WGS pipeline-complexity figures

Status: final. All six cached outputs passed semantic comparison against their
pipeline-specific direct annotation output.

## Recommended figures

- `manuscript/supplement_pipeline_complexity.*`: four-panel journal supplement
  showing relative runtime, absolute time saved, runtime-model agreement, and
  the dependency of speedup on direct pipeline runtime.
- `manuscript/panels/supplement_pipeline_complexity/`: the four manuscript
  panels as separate high-resolution files.
- `repository/headliner_pipeline_spectrum.*`: simplified measured spectrum from
  fastVEP through standard VEP to controlled all-day pipelines.
- `alternatives/pipeline_time_returned.*`: absolute cached and direct wall times
  on a log scale.

Every figure is available as SVG, PDF, and 300 dpi PNG.

## Frozen evidence

- Input: KPGP-00319, GRCh38, 4,795,706 variants.
- Cache: bundled Zenodo gnomAD AF >= 10% cache, DOI 10.5281/zenodo.18189447.
- Measured hit rate: 80.2347%.
- Six controlled VEP plugin costs: 0.5, 1, 2, 4, 7, and 10 ms per transcript
  consequence.
- Timed cells: 12, all with `--statistics light`; no technical repeats.
- Direct runtimes: 4.21 to 22.71 hours.
- Cached runtimes: 0.91 to 4.52 hours.
- Measured speedups: 4.62x to 5.07x.
- Time returned per genome: 3.30 to 18.19 hours.
- Runtime-model residuals: within 0.56 percentage points.

The synthetic plugin emits no annotation fields. It changes annotation cost
while keeping the input, cache, annotator configuration, and output semantics
fixed. fastVEP and standard VEP anchors come from the completed matched-light
52-genome annotator campaign and use the same KPGP genome/cache strategy shown
in the headline figure.

## Provenance

- Frozen inputs: `../../../source_data/final/2026-08-11T223247Z-pipeline-complexity/`
- Derived plotting tables: `source/`
- Renderer record: `RENDERED_FROM_PIPELINE.txt`
- R environment: `R_SESSION_INFO_PIPELINE.txt`
