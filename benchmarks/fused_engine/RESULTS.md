# Fused-routing prototype: local engineering results

Measured on 22 October 2026. These compare **two cached VCFcache workflows**,
not VCFcache against a direct annotator and not vep-rs. The native engine is
opt-in; the production/default engine remains legacy.

## Outcome

Across twelve matched conditions, the fused prototype reduced cached workflow
time by **24.6–60.8%**, or **1.33–2.55-fold**. Every complete-record and relevant
header comparison passed, with no annotation-field exceptions. The native
engine did not fall back to legacy in any measured condition.

| Synthetic CSQ payload | Cache hits | CPU affinity limit | Legacy (s) | Fused (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| 32 bytes | 80% | 1 | 11.20 | 7.22 | 1.55x |
| 32 bytes | 80% | 6 | 10.06 | 4.83 | 2.08x |
| 32 bytes | 80% | 12 | 9.67 | 4.22 | 2.29x |
| 32 bytes | 90% | 1 | 10.94 | 6.69 | 1.64x |
| 32 bytes | 90% | 6 | 9.11 | 4.13 | 2.21x |
| 32 bytes | 90% | 12 | 8.50 | 3.48 | 2.44x |
| 1,024 bytes | 80% | 1 | 16.86 | 10.97 | 1.54x |
| 1,024 bytes | 80% | 6 | 9.79 | 7.38 | 1.33x |
| 1,024 bytes | 80% | 12 | 11.16 | 6.64 | 1.68x |
| 1,024 bytes | 90% | 1 | 18.89 | 10.78 | 1.75x |
| 1,024 bytes | 90% | 6 | 9.66 | 7.02 | 1.37x |
| 1,024 bytes | 90% | 12 | 10.61 | 4.16 | 2.55x |

In the annotation-rich, 90%-hit condition, increasing the enforced CPU limit
from one to twelve reduced fused time from 10.78 to 4.16 seconds (**2.59-fold**).
This is useful scaling, **not linear scaling**. At twelve CPUs, measured child
CPU time also fell from 36.81 to 20.67 CPU-seconds compared with legacy.

The speed improvement trades compression work for temporary space: that cell
retained about 1.40 GB of workflow/output artifacts versus 0.30 GB for legacy
(4.7-fold). Successful ordinary CLI runs clean the work directory unless debug
mode is enabled. The recorded RSS is the largest child-process RSS, not an
aggregate simultaneous memory peak; it must not be interpreted as the latter.

## Design and interpretation

- One million generated records across 22 contigs, including biallelic SNPs
  and indels, one sample with GT/DP, and ten deterministic supplementary INFO
  fields. The payload column is the extra generated CSQ body, not total bytes
  per record. The repeated synthetic strings are compressible and are not a
  clinical annotation recipe.
- Both engines used the same indexed input, cache, bcftools lookup recipe,
  ordered output contract and compression level. The mock recipe performed a
  real INFO overlay from a complete annotation table; it did not sleep.
- Each pair ran serially, legacy then fused, with identical process-wide CPU
  affinity inherited by all children. CPUs were not exclusively reserved on
  the workstation. There was one run per condition and no inferential intervals.
- Timing covers `WorkflowManager.run`, including all routing, recipe execution,
  compression and final indexing. CLI construction, data/cache generation and
  the external complete-record comparison were outside timing. Neither arm
  performed a full-statistics rescan.
- Inputs and intermediate files were on `/tmp` (tmpfs), so these are primarily
  CPU/memory-path measurements, not evidence about SSD, network filesystem or
  cold-storage performance. This is an engineering diagnostic, not a manuscript
  cohort benchmark or a measured vep-rs break-even result.

The final implementation uses grouped, region-restricted contig routing,
uncompressed hit spools, one invocation of the annotation recipe, a miss-only
second INFO overlay, ordinal reconstruction, direct native compressed shard
writing, and header-checked BGZF concatenation without recompression.
libdeflate and normal portable release optimization are enabled; no CPU-native
compiler flags or bespoke SIMD record parser were used.

## Provenance and next validation

Authoritative source files are in [output/20261022-native-writer](output/20261022-native-writer):
`summary.tsv`, manifest with native binary and source hashes, per-cell resource
and stage reports, execution logs, and twelve complete-record equality reports.
The measured native binary/core are unchanged. Exact measured Python sources
are archived under `measured_sources/` and match the manifest hashes. Subsequent
maintenance only avoids an unnecessary missing-index diagnostic during
header-only reads and emits LF rather than CRLF in generated TSV summaries;
neither changes routing, annotations or the recorded timing values.
Earlier local development trials, including an
invalid synthetic input preparation attempt, are not part of this result set
and remain only in ignored `scratch/`.

Raw inputs, caches and workflow files are retained locally under
`/tmp/vcfcache-fused-bench-ub3duiu0`; the manifest records the path. These are
temporary, reproducible artifacts rather than permanent manuscript source data.

Reproduce with a new output directory:

```sh
python -m benchmarks.fused_engine --output /scratch/new-fused-run \
  --records 1000000 --payloads 32 1024 --hit-rates 0.8 0.9 --cpus 1 6 12
```

Next validate on existing real WGS caches and the independently prepared
vep-rs recipe, including native annotation-rich output. Do not promote the
prototype to the default or infer a vep-rs speedup from these synthetic results.
