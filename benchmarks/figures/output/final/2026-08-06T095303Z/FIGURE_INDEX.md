# Complete VEP benchmark figure set

Snapshot: `2026-08-06T095303Z` (`FINAL`)

This directory contains the complete non-FastVEP benchmark dataset. Each
assembled figure and every component panel is available as vector SVG, print
PDF, and 300 dpi PNG. All plots were rendered in R with ggplot2.

## Primary figures

- `manuscript/assay_benchmark.*` — cache impact across 20 panels, 20 WES
  samples, and 49 1000 Genomes WGS samples.
- `manuscript/external_wgs.*` — 52 independent WGS genomes from KPGP, SGDP,
  and PGP, each tested with bundled gnomAD AF >= 10%, bundled gnomAD AF >= 1%,
  and a cohort-derived three-genome cache.
- `repository/user_impact_model_preview.*` — simplified communication graphic.
  This is intentionally marked as a model preview because its time examples
  still use the zero-overhead approximation.

## Alternative manuscript views

- `alternatives/external_wgs_decision.*` — converts relative runtime into
  median wall-clock hours and measured hours saved per genome. This is the most
  direct reader-facing complement to the main external-WGS figure.
- `alternatives/runtime_model_diagnostics.*` — tests the proposed runtime model
  over assays and cohorts, and shows why fixed lookup/preprocessing overhead
  dominates very short pipelines but becomes a smaller fraction of WGS runs.

The individual panels are in the corresponding `panels/` directories.

## Dataset status

- Primary WGS: 49/49 planned successful samples; HG02888 is the documented
  cgroup-OOM exclusion.
- Assay campaign: 40/40 samples (20 panel and 20 WES).
- External WGS: 52/52 genomes and 156/156 cache-strategy outputs semantically
  validated.
- Unexpected annotation mismatches: 0.

An external evaluation genome is held out only from construction of its
cohort-derived cache. It is then annotated against that cache and both bundled
Zenodo caches. This avoids evaluating a custom cache on the same three genomes
from which it was built.

## Headline medians

- Bundled gnomAD AF >= 1% cache: 8.2x KPGP, 9.0x SGDP, and 8.5x PGP WGS
  speedup.
- Bundled gnomAD AF >= 10% cache: 4.4x KPGP and 4.9x for both SGDP and PGP.
- Three-genome cohort cache: 4.7x KPGP, 4.5x SGDP, and 3.6x PGP; measured
  cache-build cost breaks even after 2-3 subsequent genomes.
- Assay-size series using the bundled AF >= 1% cache: 2.3x median WES and
  12.3x median WGS speedup. The tiny panel workload is slower with caching,
  clearly exposing the fixed-overhead boundary.

Exact values are in `source/assay_summary.tsv`,
`source/external_wgs_summary.tsv`, and the other derived TSV files.
