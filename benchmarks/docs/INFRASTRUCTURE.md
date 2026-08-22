# Benchmark infrastructure assessment

Assessment date: 2026-07-30 UTC.

## Current preparation VM

| Resource | Observed capacity | Assessment |
|---|---:|---|
| CPU | 32 vCPUs, Broadwell/KVM | Sufficient for concurrent download and VCF preparation |
| RAM | 62 GiB, no swap | Sufficient; preparation peaked near 6.3 GiB active RAM |
| Root filesystem | 19 GiB total, 3.4 GiB free | Not suitable for benchmark data or temporary files |
| `/mnt/data` | 2.0 TiB total, 592 GiB free before setup | Sufficient with a wide safety margin |
| Source transfer | approximately 9 MiB/s in the initial single-stream probe | Multiple resumable streams are appropriate |

All large inputs, shards, sort files, and final VCFs live under
`/mnt/data/vcfcache_benchmarks`. The preparation program rejects a data root
outside `/mnt/data`, requires at least 100 GiB free, and places `TMPDIR` below
the benchmark root.

The 1000 Genomes autosomal sources are approximately 28 GiB. Seven prepared
GIAB VCFs occupy approximately 953 MiB. Even allowing for source files,
per-chromosome shards, final VCFs, partial downloads, and sort amplification,
the preparation peak is expected to remain below 75 GiB.

The completed setup used 37 GiB: 29 GiB of immutable sources, 3.3 GiB of
1000 Genomes chromosome shards, 3.3 GiB of final 1000 Genomes VCFs, and
920 MiB of final GIAB VCFs. `/mnt/data` retained 556 GiB free. From freezing
the sample manifest at 14:42 UTC to the final QC report at 15:46 UTC, setup
took 64 minutes on this VM, including the long-tail public-server transfer.

The assay expansion needs approximately 5.0 GiB of additional large sources
(HPRC R2 wave VCF and 1000 Genomes chrX), 105 MiB for Ensembl 115, and less
than 10 MiB for interval definitions. Derived WES, panel, chrX-shard, and HPRC
files are expected to remain below 10 GiB. With 552 GiB free at the expansion
preflight, preparation remains comfortably within the current VM; no Slurm
cluster is needed for this data-construction step.

Actual expansion setup completed in approximately 25 minutes and left 544 GiB
free. The complete benchmark root is 48 GiB (34 GiB sources, 5.0 GiB work,
5.8 GiB final samples, and 3.8 GiB retained pilots). The current VM was more
than sufficient; Slurm remains necessary only for the repeated annotation
matrix.

## Decision

Keep public-data preparation on the current VM. Moving this streaming workload
to a cluster would add setup and transfer overhead without solving a current
resource constraint.

The paired HG02079 pilot measured 2 h 07 min 47 s uncached and 15 min 25 s
cached. Each retained run directory occupies 969 MiB. The workflow used eight
VEP forks and the current VM completed both paths safely.

Both original modes reported approximately 44.1 GiB peak RSS. Inspection
showed that this peak occurred after annotation in final statistics accounting,
which captured a complete `bcftools view` stream in Python.

Commit `88c018086b21` replaced that operation with a newline-only,
incrementally counted `bcftools query`. A direct 4,329,621-record test used
31 MiB peak RSS and returned the exact count. In the repeated end-to-end pair,
GNU time's largest-process peak fell to 5.22 GiB uncached and 749 MiB cached.
Manual snapshots during uncached VEP observed approximately 20–23 GiB aggregate
worker RSS; GNU time does not sum concurrent child RSS.

Keep sequential pilots on the current VM. Use Slurm for the 50-sample,
single-pass annotation matrix. A conservative initial request is eight CPUs
and 32 GiB RAM per job, allowing 15 concurrent jobs within the project's
500-GiB aggregate memory. Confirm cgroup peak memory on the first Slurm jobs
before increasing concurrency.
