# Frozen paired publication pilot

Pilot date: 2026-07-30 UTC.

## Configuration

| Item | Frozen value |
|---|---|
| Sample | 1000 Genomes HG02079 (EAS), GRCh38 |
| Input variants | 4,329,621 |
| VCFcache code | commit `0d1d794242ff` |
| Annotation | VEP 115.2, `--everything`, eight forks |
| Public cache | gnomAD v4.1 GRCh38, AF ≥ 0.001 |
| Host | 32-vCPU, 62-GiB VM |
| Replicates | one paired pilot; not an inferential result |

The cached and uncached commands were identical except for the `--uncached`
flag. GNU time covered the complete VCFcache command, including normalization,
annotation, output writing, indexing, and final accounting.

## Result

| Metric | Uncached | Cached | Change |
|---|---:|---:|---:|
| Wall time | 7,666.6 s (2:07:46.6) | 924.6 s (0:15:24.6) | 8.29×; 87.94% saved |
| User CPU time | 56,614.6 s (15.73 h) | 4,624.3 s (1.28 h) | 14.44 CPU-h saved |
| Variants sent to VEP | 4,329,621 | 172,755 | 4,156,866 avoided |
| Effective cache-hit rate | — | 96.01% | — |
| Peak RSS | 44.1 GiB | 44.1 GiB | accounting-limited |
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

## Scale-up gate

Before launching the 50-sample matrix:

1. change final annotated-variant accounting to stream/count records rather
   than retaining full `bcftools view` output;
2. rerun this exact pair and semantic validation;
3. freeze that commit and measure its true per-job peak RSS;
4. size Slurm concurrency from the new peak and execute one warm-up plus three
   randomized measured pairs per sample.

Machine-readable results are outside Git at
`/mnt/data/vcfcache_benchmarks/pilot/HG02079/0d1d794242ff/summary.json`.
