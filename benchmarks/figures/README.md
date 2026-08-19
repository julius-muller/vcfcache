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

## Render the fastVEP and paired annotator snapshot

The complete 52-genome fastVEP campaign is frozen separately so that it does
not overwrite the earlier VEP figure set:

```bash
Rscript --vanilla benchmarks/figures/R/render_annotators.R \
  benchmarks/figures/source_data/final/2026-08-06T223844Z-fastvep \
  benchmarks/figures/output/final/2026-08-06T223844Z-fastvep
```

This produces the final standalone `manuscript/fastvep_real_world_wgs.*` and a
clearly labelled provisional `manuscript/vep_fastvep_impact_provisional.*`.
The paired view remains provisional only because the historical VEP campaign
timed `--statistics full`; it can be refreshed after the small VEP light-mode
calibration. Individual high-resolution panels are stored below each figure's
`manuscript/panels/` subdirectory.

## Render the complete final figure round

The completed VEP light-statistics calibration, final fastVEP data, original
VEP/assay figures, and new repository headline alternatives are assembled in a
new directory with two R render passes:

```bash
Rscript --vanilla benchmarks/figures/R/render_all.R \
  benchmarks/figures/source_data/final/2026-08-06T095303Z \
  benchmarks/figures/output/final/2026-08-07T053247Z-complete

Rscript --vanilla benchmarks/figures/R/render_annotators.R \
  benchmarks/figures/source_data/final/2026-08-07T053247Z-complete \
  benchmarks/figures/output/final/2026-08-07T053247Z-complete
```

The recommended journal outputs are under `manuscript/`; their component
panels are under `manuscript/panels/`. Three contrasting GitHub-ready headline
styles are under `repository/`: measured before/after bars, an
empirical-overhead scenario map, and a workflow-style explainer. Read the
directory's `FIGURE_INDEX.md` before selecting figures because it distinguishes
measurements, sensitivity calibration, and modeled scenarios.

## Render the matched-light final figure round

After the complete 52-genome VEP light-mode rerun, the recommended final set is
rendered from measured light-mode data for both annotators:

```bash
Rscript --vanilla benchmarks/figures/R/render_all.R \
  benchmarks/figures/source_data/final/2026-08-09T223718Z-light-complete-vep \
  benchmarks/figures/output/final/2026-08-09T223718Z-light-matched-final

Rscript --vanilla benchmarks/figures/R/render_annotators.R \
  benchmarks/figures/source_data/final/2026-08-09T223718Z-light-matched-final \
  benchmarks/figures/output/final/2026-08-09T223718Z-light-matched-final
```

Use this directory in preference to the 2026-08-07 set. The VEP and fastVEP
comparison no longer has a timing-mode asymmetry: both use
`--statistics light`. The historical VEP full-rescan data appear only in the
paired sensitivity panel and are never used to rescale the new observations.

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

## Render the real-WGS pipeline-complexity spectrum

The controlled pipeline-cost campaign uses one 4.80-million-variant KPGP WGS,
the bundled Zenodo gnomAD AF >= 10% cache, and six no-output synthetic plugin
delays. All timed cells use light statistics and every cached output must pass
semantic comparison before the collector accepts the campaign.

```bash
python3 benchmarks/figures/collect_pipeline_complexity_snapshot.py \
  --wgs-spectrum-campaign-root <completed-campaign-root> \
  --annotator-snapshot \
    benchmarks/figures/source_data/final/2026-08-09T223718Z-light-matched-final \
  --output <new-source-snapshot-directory>

Rscript --vanilla benchmarks/figures/R/render_pipeline_complexity.R \
  <new-source-snapshot-directory> <new-output-directory>
```

The manuscript supplement is a four-panel view of relative runtime, absolute
hours returned, model agreement, and speedup saturation. Individual panels are
exported separately under `manuscript/panels/supplement_pipeline_complexity/`.
The same render also produces a repository headline spanning fastVEP, standard
VEP, and the controlled heavier pipelines, plus an absolute-time alternative.

## Render the matched dual-cache pipeline spectrum

After the cached-only AF >= 1% extension completes, freeze it together with the
existing AF >= 10% snapshot. The collector requires the same six direct
baselines, six successful AF >= 1% semantic fingerprints, the production
Zenodo alias/DOI, light statistics, and zero technical repeats.

```bash
python3 benchmarks/figures/collect_pipeline_complexity_dual_snapshot.py \
  --af10-snapshot <final-af10-pipeline-snapshot> \
  --af1-campaign-root <completed-af1-campaign-root> \
  --output <new-dual-cache-source-snapshot>

Rscript --vanilla benchmarks/figures/R/render_pipeline_complexity_dual.R \
  <new-dual-cache-source-snapshot> <new-output-directory>
```

This produces a four-panel journal supplement comparing runtime remaining,
absolute hours returned, speedup ceilings, and runtime-model agreement for both
bundled cache thresholds. Individual panels, a simplified cache-coverage
headline, and the incremental benefit of the AF >= 1% cache are also exported.

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
- The older `user_impact_model_preview` uses the zero-overhead form of the
  runtime rule. The newer `headliner_scenario_map` instead uses annotator-specific
  median overhead inferred from the 52 real-world genomes and labels its cells
  as modeled scenarios.
