# Exploratory fastVEP pilot

This pilot asks two deliberately narrow questions:

1. Does fastVEP produce credible core consequence annotations on a small human
   test input when compared with Ensembl VEP using the same GFF3 and FASTA?
2. Does the current VCFcache engine materially reduce total fastVEP wall time
   on one real WGS at controlled 80%, 90%, and 100% cache-hit rates?

It is an engineering diagnostic, not a publication benchmark. It runs once per
cell, uses Ensembl release 115 assets already available on ITCCcloud_dev, and
records enough provenance to interpret the outcome without attempting a fully
portable environment lock.

The completed exploratory outcome is documented in
[FASTVEP_PILOT_RESULTS.md](FASTVEP_PILOT_RESULTS.md). The rendered diagnostic
and its separate high-resolution panels are under `output/`.

## Frozen exploratory configuration

- host: `appuser@10.133.255.21` (`gvbrowse-preproc`, ITCCcloud_dev);
- runtime root: `/mnt/data/vcfcache_benchmarks/fastvep_pilot`;
- input: normalized, biallelic HG02079 GRCh38 WGS;
- fastVEP: 0.3.0, commit `e47216cebe3abcd8dff722b7fb0ab1b19d4fcc80`;
- transcript data: Ensembl release 115 GFF3 and matching GRCh38 FASTA;
- threads: 16 for Rayon and bcftools;
- timed WGS cells: direct, 100%, 90%, and 80% cache hits.

Both VCFcache arms use `--skip-split-multiallelic`. The frozen input contains
no multiallelic or spanning-deletion records, so this removes redundant
normalization symmetrically without changing the input contract.

## Commands

Run from the repository checkout on the pilot host:

```bash
.venv/bin/python benchmarks/fastvep_pilot/run.py setup
.venv/bin/python benchmarks/fastvep_pilot/run.py smoke
.venv/bin/python benchmarks/fastvep_pilot/run.py wgs
.venv/bin/python benchmarks/fastvep_pilot/run.py collect
```

`all` runs those four phases in order. Completed phases are resumable. Existing
run directories are never silently overwritten.

The runner writes raw logs, resource measurements, equality reports,
`summary.tsv`, `summary.json`, and `FASTVEP_PILOT_RESULTS.md` below the runtime
root. Plot after copying `summary.tsv` back to this checkout:

```bash
Rscript --vanilla benchmarks/fastvep_pilot/plot.R \
  /path/to/summary.tsv \
  benchmarks/fastvep_pilot/output
```

Every cached result must match the direct result after canonicalizing only INFO
tag order and CSQ entry order. Variant order, all ordinary VCF fields, FORMAT,
and sample values remain exact.
