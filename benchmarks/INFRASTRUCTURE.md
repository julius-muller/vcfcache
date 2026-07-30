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

## Decision

Keep public-data preparation on the current VM. Moving this streaming workload
to a cluster would add setup and transfer overhead without solving a current
resource constraint.

The paired HG02079 pilot measured 2 h 07 min 47 s uncached and 15 min 25 s
cached. Each retained run directory occupies 969 MiB. The workflow used eight
VEP forks and the current VM completed both paths safely.

Both modes reported approximately 44.1 GiB peak RSS. Inspection showed that
this peak occurs after annotation in the final statistics accounting, which
captures a complete `bcftools view` stream in Python. It is not the working set
required by VEP. This should be converted to streaming record counting and the
paired pilot repeated before sizing the full matrix.

Keep sequential pilots on the current VM. Use Slurm for the 50-sample,
three-replicate annotation matrix after the memory fix. Until then, do not
schedule more than one benchmark job per 62-GiB worker; the project's 500 GiB
aggregate memory would otherwise limit safe concurrency to roughly ten jobs.
