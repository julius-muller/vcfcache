# Benchmark figures

This workspace separates frozen benchmark snapshots from plotting code and
rendered drafts. All plots are produced in **R with ggplot2**. Python is used
only to collect compact metrics and provenance from completed campaign output.

## Layout

- `collect_snapshot.py`: extracts TSV/JSON source data without copying BCFs.
- `R/`: ggplot2 renderers and the package setup helper.
- `source_data/preliminary/<timestamp>/`: immutable preliminary snapshots.
- `output/preliminary/<timestamp>/source/`: derived figure source tables.
- `output/preliminary/<timestamp>/manuscript/`: evidence-heavy figures.
- `output/preliminary/<timestamp>/repository/`: simplified user-facing figure.

Every snapshot records campaign IDs, completion counts, hashes, and validation
status in `SNAPSHOT.json`. External rows awaiting corrected semantic comparison
are retained for exploratory plotting but explicitly marked
`comparator_pending_*`. They must not be cited or promoted to final results.

## Render the current snapshot

Install the small R dependency set once:

```bash
Rscript --vanilla benchmarks/figures/R/setup_packages.R
```

Render a selected immutable snapshot:

```bash
Rscript --vanilla benchmarks/figures/R/render_all.R \
  benchmarks/figures/source_data/preliminary/2026-08-05T070826Z \
  benchmarks/figures/output/preliminary/2026-08-05T070826Z
```

Each figure is written as SVG, PDF, and PNG. Derived summary/model tables and
the R session information are written beside the figures.

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
and rerun `render_all.R`. Never overwrite a previous snapshot. The final data
refresh will use `source_data/final/`, require all expected samples, and reject
every external row that has not passed corrected semantic post-processing.

## Interpretation of the current draft

- The 49 primary WGS and 40 panel/WES pairs are semantically validated.
- Primary 1000 Genomes WGS illustrates assay size but is not the independent
  real-world WGS hit-rate estimate.
- The external figure is deliberately watermarked preliminary until all 52
  genomes complete and the comparator fixes have been applied.
- The repository graphic currently uses the zero-overhead form of the runtime
  rule. The controlled-runtime campaign will replace that placeholder with
  measured lookup and preprocessing overhead.
