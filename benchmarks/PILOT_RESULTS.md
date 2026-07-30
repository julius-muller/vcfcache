# Confirmed paired publication pilot

Pilot date: 2026-07-30 UTC.

## Configuration

| Item | Frozen value |
|---|---|
| Sample | 1000 Genomes HG02079 (EAS), GRCh38 |
| Input variants | 4,329,621 |
| Initial VCFcache code | commit `0d1d794242ff` |
| Confirmed VCFcache code | commit `88c018086b21` |
| Annotation | VEP 115.2, `--everything`, eight forks |
| Public cache | gnomAD v4.1 GRCh38, AF ≥ 0.001 |
| Host | 32-vCPU, 62-GiB VM |
| Replicates | one paired pilot; not an inferential result |

The cached and uncached commands were identical except for the `--uncached`
flag. GNU time covered the complete VCFcache command, including normalization,
annotation, output writing, indexing, and final accounting.

## Confirmed result

| Metric | Uncached | Cached | Change |
|---|---:|---:|---:|
| Wall time | 7,613.4 s (2:06:53.4) | 796.1 s (0:13:16.1) | 9.56×; 89.54% saved |
| User CPU time | 57,705.5 s (16.03 h) | 4,414.5 s (1.23 h) | 14.80 CPU-h saved |
| Variants sent to VEP | 4,329,621 | 172,755 | 4,156,866 avoided |
| Effective cache-hit rate | — | 96.01% | — |
| Largest-process peak RSS | 5.22 GiB | 749 MiB | statistics spike removed |
| Retained run directory | 969 MiB | 969 MiB | no material difference |

This single-sample pilot supports the mechanism and establishes feasibility. It
does not support population-level confidence intervals or a final headline
effect size; those require the prespecified sample-level replicates.

## Correctness gate

The semantic comparator streamed all 4,329,621 output records and checked
CHROM/POS/REF/ALT, input AF/AC/AN, genotype, the CSQ header, and the complete
CSQ item set. Results:

- zero variant or input-field mismatches;
- zero annotation mismatches;
- identical canonical semantic SHA-256:
  `c70dfac9c066aeac238b7089ed799e6a5af8bdf55f8e92afa8b2bd32fc6821a9`;
- 1,884 multi-allelic loci differed only in split-record order.

The first validator pass compared split records positionally and exposed the
within-locus ordering difference. A regression test was added before rerunning
the complete comparison; the revised comparator canonicalizes record order
only within the same CHROM/POS locus.

## Streaming-accounting confirmation

The initial pair used a statistics implementation that captured every matching
VCF text record in memory and produced a 44.1-GiB peak in both modes. The
streaming implementation emits only one newline per matching record and counts
the iterator incrementally.

- A direct full-output count returned 4,329,621 in 38.3 s at 31 MiB peak RSS.
- Uncached end-to-end peak fell from 44.1 GiB to 5.22 GiB.
- Cached end-to-end peak fell from 44.1 GiB to 749 MiB.
- Cached wall time fell by 128.4 s; uncached wall time fell by 53.2 s.
- The repeated pair retained exact semantic equivalence.

GNU time reports the maximum RSS of the largest process, not the sum of
concurrent VEP workers. Read-only snapshots observed approximately 20–23 GiB
aggregate worker RSS during uncached annotation. Slurm jobs should initially
request 32 GiB and use cgroup memory metrics for final reporting.

## Scale-up

The implementation is ready for the prespecified 50-sample matrix: one warm-up
and three randomized measured cached/uncached pairs per sample. Start with 15
concurrent eight-CPU, 32-GiB Slurm jobs and adjust only after observing the
first batch's cgroup memory.

Machine-readable results are outside Git at
`/mnt/data/vcfcache_benchmarks/pilot/HG02079/88c018086b21/summary.json`.
