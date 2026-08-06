# Source-overlap-corrected VEP benchmark figures

Snapshot: `2026-08-06T095303Z` (`FINAL`)

This is the recommended non-FastVEP figure set. It was rendered into a new
directory so that the earlier output remains unchanged.

## Evidence roles

### Primary real-world evidence

- `manuscript/external_wgs.*` — measured cache reuse and speedup in 52 KPGP,
  SGDP, and PGP genomes with no documented project overlap with gnomAD.
- `alternatives/external_wgs_decision.*` — the same evidence expressed as
  wall-clock time and hours returned per genome.

The evaluation genomes were excluded from the three genomes used to build each
cohort-derived cache. Because gnomAD does not disclose all contributing
projects or identities, these cohorts are described as having no documented
project overlap, not as proven absent from gnomAD at the individual level.

### Runtime calibration only

- `manuscript/assay_source_overlap_calibration.*` — 20 panel, 20 WES, and 49
  WGS inputs derived from 1000 Genomes identities.
- `alternatives/runtime_model_diagnostics.*` — uses those calibration points
  together with the external cohorts to test runtime conditional on observed
  cache hit rate.

The gnomAD source universe includes 1000 Genomes samples. Consequently, the
1000 Genomes-derived cache hit rates and speedups are source-overlap upper
bounds. They must not be quoted as estimates for an independent user cohort.
They remain valid for demonstrating fixed lookup overhead and the relationship
between observed miss fraction and runtime.

### User-facing model

- `repository/user_impact_model_preview.*` — illustrative calculations over
  user-selected hit rates and pipeline durations. It remains explicitly marked
  as a model preview and does not derive its scenario hit rates from the
  overlapping 1000 Genomes cohort.

## Headline real-world WGS medians

- Bundled gnomAD any-stratum AF >= 1% cache: 8.2x KPGP, 9.0x SGDP, and 8.5x
  PGP speedup.
- Bundled gnomAD any-stratum AF >= 10% cache: 4.4x KPGP and 4.9x for both SGDP
  and PGP.
- Three-genome cohort cache: 4.7x KPGP, 4.5x SGDP, and 3.6x PGP; measured
  cache-build cost breaks even after 2-3 subsequent genomes.

All 52 external genomes and all 156 cached strategy outputs passed semantic
validation. Every composite and component panel is available as SVG, PDF, and
300 dpi PNG. Exact values are provided under `source/`.
