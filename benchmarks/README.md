# VCFcache publication benchmarks

This directory contains the reproducible preparation and execution plan for the
VCFcache publication benchmarks. Large datasets and results are kept outside Git
under `/mnt/data/vcfcache_benchmarks`.

## Documents

- [PLAN.md](PLAN.md): implementation phases and publication figure specification.
- [DATA_SETUP.md](DATA_SETUP.md): exact public-data layout, transformations, and QC.
- [INFRASTRUCTURE.md](INFRASTRUCTURE.md): measured VM capacity and scaling decision.
- [PILOT_RESULTS.md](PILOT_RESULTS.md): frozen paired-pilot results and scale-up gate.

## Prepare the public VCF cohort

```bash
.venv/bin/python benchmarks/prepare_data.py preflight
.venv/bin/python benchmarks/prepare_data.py select
.venv/bin/python benchmarks/prepare_data.py smoke
.venv/bin/python benchmarks/prepare_data.py all
```

The default cohort contains 50 unrelated 1000 Genomes samples (10 from each
superpopulation, sex-balanced) and the seven GIAB v4.2.1 GRCh38 benchmark
genomes. Final inputs are sorted, BGZF-compressed, tabix-indexed, single-sample
VCFs covering autosomes 1–22.

All stages are resumable. Successfully verified downloads and completed
chromosome/sample outputs are skipped on subsequent invocations.

After QC, all files are also available through the flat, non-duplicating
directory `/mnt/data/vcfcache_benchmarks/samples/GRCh38/all`.

## Tests

```bash
.venv/bin/python -m pytest tests/test_benchmark_prepare.py -q
.venv/bin/python -m pytest -q
VCFCACHE_BENCHMARK_NETWORK=1 \
  .venv/bin/python -m pytest tests/test_benchmark_prepare.py -m benchmark_network -q
```

## Paired publication pilot

The runner uses the local gnomAD GRCh38 AF≥0.001 cache and its immutable VEP
115.2 `--everything` recipe. Cached and uncached commands are identical except
for VCFcache's `--uncached` switch.

```bash
.venv/bin/python benchmarks/run_pilot.py preflight
.venv/bin/python benchmarks/run_pilot.py all
```

Each run records the exact command, Git commit, environment, wall/CPU time,
maximum RSS, filesystem counters, VCFcache variant counts, and cache hit rate.
The final comparison streams both BCFs and requires identical variant keys,
input AF/AC/AN and GT values, and equivalent CSQ sets. It canonicalizes only
split-allele order within a locus and CSQ item order.

## Slurm cohort

Generate the frozen 50-sample, three-replicate task matrix on the controller:

```bash
.venv/bin/python benchmarks/run_cohort.py prepare
.venv/bin/python benchmarks/run_cohort.py submit
.venv/bin/python benchmarks/run_cohort.py status
.venv/bin/python benchmarks/run_cohort.py collect
```

For the publication infrastructure smoke test, prepare only the confirmed
HG02079 pair:

```bash
.venv/bin/python benchmarks/run_cohort.py prepare \
  --sample HG02079 --replicates 1 \
  --tasks /mnt/data/vcfcache_benchmarks/manifests/slurm_smoke.tsv
.venv/bin/python benchmarks/run_cohort.py submit \
  --tasks /mnt/data/vcfcache_benchmarks/manifests/slurm_smoke.tsv \
  --concurrency 1
```

Each array task runs cached and uncached modes in a deterministic randomized
order on worker-local storage, requires semantic equivalence, and archives the
pair atomically below `/results/tasks/`.
