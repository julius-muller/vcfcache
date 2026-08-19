# Final matched dual-cache pipeline-complexity figures

Status: final. The AF >= 10% campaign contributed six paired cached/direct
measurements. The AF >= 1% extension contributed six new cached measurements
against the same frozen direct baselines. All 12 cached outputs passed their
specified semantic validation.

## Recommended figures

- `manuscript/supplement_dual_cache_pipeline.*`: four-panel journal supplement
  comparing wall time remaining, absolute hours returned, measured speedup
  ceilings, and runtime-model agreement.
- `manuscript/panels/supplement_dual_cache_pipeline/`: all four manuscript
  panels as separate high-resolution files.
- `repository/headliner_cache_coverage.*`: simplified comparison of standard,
  heavier, and all-day pipelines using both bundled cache thresholds.
- `alternatives/af1_additional_hours.*`: additional time returned by the AF >=
  1% cache beyond the AF >= 10% cache.

Every figure is available as SVG, PDF, and 300 dpi PNG.

## Main measured result

| Bundled cache | Hit rate | Wall time remaining | Speedup | Maximum hours returned |
|---|---:|---:|---:|---:|
| gnomAD any-stratum AF >= 10% | 80.23% | 19.7-21.7% | 4.62-5.07x | 18.19 h |
| gnomAD any-stratum AF >= 1% | 90.26% | 9.9-11.9% | 8.39-10.10x | 20.45 h |

For the 22.7-hour direct pipeline, the measured remaining wall time was 4.52
hours with the AF >= 10% cache and 2.25 hours with the AF >= 1% cache. Across
all 12 measurements, runtime-model residuals were within 0.56 percentage
points.

## Evidence and provenance

- Input: KPGP-00319, GRCh38, 4,795,706 variants.
- AF >= 10% cache DOI: `10.5281/zenodo.18189447`.
- AF >= 1% cache DOI: `10.5281/zenodo.18190046`.
- Six controlled VEP loads: 0.5, 1, 2, 4, 7, and 10 ms per transcript
  consequence.
- Timed cells used `--statistics light`; no technical repeats.
- The synthetic plugin emits no annotation fields.
- AF >= 1% outputs share one canonical semantic and CSQ-header fingerprint and
  are anchored to the prior direct comparison of the same sample/public cache.
- Frozen source snapshot:
  `../../../../source_data/final/2026-08-12T223200Z-pipeline-complexity-dual-cache/`.
- Derived plotting tables: `source/`.
- Renderer record: `RENDERED_FROM_DUAL_CACHE.txt`.
- R environment: `R_SESSION_INFO_DUAL_CACHE.txt`.
