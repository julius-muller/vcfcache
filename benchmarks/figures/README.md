# Benchmark figures

This workspace separates frozen benchmark snapshots from plotting code and
rendered figures. All plots are produced in **R with ggplot2**. Python is used
only to collect compact metrics and provenance from completed campaign output.

## Layout

- `collect_snapshot.py`: extracts TSV/JSON source data without copying BCFs.
- `R/`: ggplot2 renderers and the package setup helper.
- `source_data/{preliminary,final}/<timestamp>/`: immutable metric snapshots.
- `output/<status>/<timestamp>/source/`: derived figure source tables.
- `output/<status>/<timestamp>/manuscript/`: evidence-heavy figures.
- `output/<status>/<timestamp>/alternatives/`: optional manuscript views.
- `output/<status>/<timestamp>/repository/`: simplified user-facing figure.

Each assembled figure is written as SVG, PDF, and 300 dpi PNG. Its component
plots are exported in the same three formats under `panels/<figure-name>/`, so
individual manuscript panels can be placed or revised without cropping the
assembled figure.

Every snapshot records campaign IDs, completion counts, hashes, and validation
status in `SNAPSHOT.json`. External rows awaiting corrected semantic comparison
may be retained in a preliminary snapshot but are explicitly marked
`comparator_pending_*`. A final snapshot requires all expected external genomes
and successful semantic validation for every cached output.

## Render the complete VEP snapshot

Install the small R dependency set once:

```bash
Rscript --vanilla benchmarks/figures/R/setup_packages.R
```

Render the immutable final snapshot:

```bash
Rscript --vanilla benchmarks/figures/R/render_all.R \
  benchmarks/figures/source_data/final/2026-08-06T095303Z \
  benchmarks/figures/output/final/2026-08-06T104928Z-source-overlap-corrected
```

Each figure is written as SVG, PDF, and PNG. Derived summary/model tables and
the R session information are written beside the figures.

The final set contains:

- `manuscript/assay_source_overlap_calibration.*`: runtime mechanics across
  1000 Genomes-derived panel, WES, and WGS inputs. These are calibration data,
  not real-world public-cache hit-rate or speedup estimates.
- `manuscript/external_wgs.*`: three independent WGS cohorts and three cache
  strategies, including custom-cache build-cost amortization.
- `alternatives/external_wgs_decision.*`: absolute wait, observed reuse, and
  hours returned per genome.
- `alternatives/runtime_model_diagnostics.*`: empirical evaluation of the
  miss-fraction scaling rule and lookup overhead.
- `repository/user_impact_model_preview.*`: the deliberately simple user-facing
  model graphic; it remains labelled as modeled until controlled-runtime data
  replace its zero-overhead assumption.

## Refresh after more jobs finish

Run the collector on `sl-head` against the three campaign roots:

```bash
python3 benchmarks/figures/collect_snapshot.py \
  --primary /mnt/data/slurm-results/campaigns/primary-wgs-32a760323fc5-bundled-af1 \
  --assay /mnt/data/slurm-results/campaigns/assay-singlepass-32a760323fc5 \
  --external /mnt/data/slurm-results/campaigns/external-wgs-6c422804f208-hg19 \
  --output /tmp/vcfcache-figure-snapshot-<UTC timestamp>
```

Copy that compact directory to a new timestamp under `source_data/preliminary/`
or `source_data/final/` and rerun `render_all.R`. Never overwrite a previously
published snapshot. The collector only assigns `FINAL` when all expected
external tasks and all three strategy rows per task are semantically validated.

## Interpretation of the complete VEP dataset

- The non-FastVEP campaigns are complete: 49 primary WGS samples, 20 panels,
  20 WES samples, and 52 independent real-world WGS genomes. HG02888 is the
  documented cgroup-OOM exclusion from the otherwise completed primary set.
- All 89 primary/assay cached-versus-direct pairs and all 156 external cached
  outputs pass semantic comparison; no unexpected annotation mismatch remains.
- The primary, WES, and panel inputs all derive from 1000 Genomes identities
  represented in the gnomAD source universe. They are retained only to examine
  runtime mechanics conditional on observed hit rate and must not support
  real-world cache-coverage or speedup claims.
- An external "evaluation genome" is a genome that was not one of the three
  genomes used to construct that cohort's custom cache. This prevents the
  custom-cache results from measuring reuse of the cache-building samples
  themselves. For bundled gnomAD caches, the plot reports ordinary lookup
  against the published cache rather than asserting individual-level exclusion.
- KPGP, SGDP, and PGP have no documented project overlap with gnomAD and carry
  the real-world WGS claims. Because gnomAD does not disclose every contributor,
  they are not described as proven absent at the individual level.
- The repository graphic currently uses the zero-overhead form of the runtime
  rule. The controlled-runtime campaign will replace that placeholder with
  measured lookup and preprocessing overhead.
