# VCFcache publication benchmarks

This directory contains the reproducible preparation and execution plan for the
VCFcache publication benchmarks. Large datasets and results are kept outside Git
under `/mnt/data/vcfcache_benchmarks`.

## Documents

- [MATERIALS_AND_METHODS.md](MATERIALS_AND_METHODS.md): publication-ready
  methods draft, with explicit placeholders for the final full-cohort run.
- [SOURCE_PROVENANCE.md](SOURCE_PROVENANCE.md): source releases, citations,
  checksum chain, data-use notes, and submission-time provenance gate.
- [manifests/](manifests/): tracked machine-readable source, sample-selection,
  and interval snapshots for peer review and deposition.
- [PLAN.md](PLAN.md): implementation phases and publication figure specification.
- [DATA_SETUP.md](DATA_SETUP.md): exact public-data layout, transformations, and QC.
- [INFRASTRUCTURE.md](INFRASTRUCTURE.md): measured VM capacity and scaling decision.
- [PILOT_RESULTS.md](PILOT_RESULTS.md): frozen paired-pilot results and scale-up gate.
- [ASSAY_DATA_RESULTS.md](ASSAY_DATA_RESULTS.md): completed HPRC/WES/panel counts,
  checksums, storage, and QC gate.
- [CACHE_STRATEGY_COMPARISON.md](CACHE_STRATEGY_COMPARISON.md): small,
  leakage-free comparison of gnomAD AF thresholds and a three-genome custom
  cache, including the vector-figure workflow.

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

## Prepare the assay-size and HPRC R2 cohorts

The second, independent pipeline leaves the original 57-sample cohort
immutable and adds three publication cohorts:

- 20 population-stratified HPRC Release 2 graph-derived GRCh38 VCFs;
- 50 matched capture-like WES VCFs using the official Twist Human Core Exome
  hg38 targets padded by 125 bp, calibrated to approximately 80,000 variants;
- 50 matched strict-target WES controls using the unpadded Twist BED; and
- 50 matched small-panel VCFs covering the 84 ACMG SF v3.3 genes, defined as
  Ensembl 115 MANE Select coding sequence plus 20 bp on each side.

The matched WES and panel cohorts use chromosomes 1–22 and X. The official
1000 Genomes chrX callset is prepared as assay-only shards, so the established
autosomal WGS cohort does not change.

```bash
.venv/bin/python benchmarks/prepare_assay_data.py preflight
.venv/bin/python benchmarks/prepare_assay_data.py download --workers 3
.venv/bin/python benchmarks/prepare_assay_data.py select-hprc
.venv/bin/python benchmarks/prepare_assay_data.py regions
.venv/bin/python benchmarks/prepare_assay_data.py prepare-x
.venv/bin/python benchmarks/prepare_assay_data.py prepare-hprc
.venv/bin/python benchmarks/prepare_assay_data.py prepare-wes --workers 4
.venv/bin/python benchmarks/prepare_assay_data.py prepare-wes-targets --workers 4
.venv/bin/python benchmarks/prepare_assay_data.py prepare-panel --workers 4
.venv/bin/python benchmarks/prepare_assay_data.py qc
```

`prepare_assay_data.py all --workers 4` runs the same resumable sequence.
All final files remain ordinary sorted, BGZF-compressed, tabix-indexed,
single-sample `.vcf.gz` files; no BCF dataset is created.

The primary `wes_twist_core` cohort uses a transparent ±125-bp capture
footprint: 73,266–95,429 records per sample (median 77,056; mean 80,211). The
unpadded 33-Mb files are retained as `wes_twist_core_targets`, a mechanistic
target-only control rather than a claim about typical raw WES output.

## Tests

```bash
.venv/bin/python -m pytest tests/test_benchmark_prepare.py -q
.venv/bin/python -m pytest tests/test_benchmark_assay_prepare.py -q
.venv/bin/python -m pytest -q
VCFCACHE_BENCHMARK_NETWORK=1 \
  .venv/bin/python -m pytest tests/test_benchmark_prepare.py -m benchmark_network -q
```

## Paired publication pilot

The runner uses the local gnomAD GRCh38 AF≥0.01 cache and its immutable VEP
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
split-allele order within a locus and CSQ item order. Differences confined to
the VEP `HGNC_ID` CSQ field are counted and reported but do not fail the
comparison because VEP 115.2 assigns that field non-deterministically according
to buffering and batch composition (Ensembl VEP issue
[#1959](https://github.com/Ensembl/ensembl-vep/issues/1959)). Any difference in
another CSQ field still fails the comparison.

The single-sample pilot is feasibility evidence only. Before manuscript
submission, replace all square-bracketed fields in
`MATERIALS_AND_METHODS.md` from the full Slurm archive and complete the
submission-time gate in `SOURCE_PROVENANCE.md`.

## Slurm cohort

Prepare and submit one campaign from the controller. This creates a one-pair
HG02079 smoke job, a 50-pair warm-up array, and a 150-pair measured array. The
arrays are linked with `afterok` dependencies, so measured jobs cannot start
unless the smoke and all warm-ups succeed.

```bash
.venv/bin/python benchmarks/run_cohort.py prepare \
  --campaign-id primary-wgs-<commit>
.venv/bin/python benchmarks/run_cohort.py submit-chain \
  --campaign-id primary-wgs-<commit> --concurrency 6
.venv/bin/python benchmarks/run_cohort.py status \
  --campaign-id primary-wgs-<commit>
.venv/bin/python benchmarks/run_cohort.py collect \
  --campaign-id primary-wgs-<commit>
```

Each task runs cached and uncached modes in separate Slurm steps on worker-local
storage, requires semantic equivalence, and atomically archives the full pair
below `/results/campaigns/<campaign-id>/`. Manifests live on the same shared
export. The controller sees it at `/mnt/data/slurm-results`; workers mount it as
`/results`. Failed attempts are retained separately and can be resubmitted with
the `submit` subcommand and a sparse `--task-ids` array specification.
