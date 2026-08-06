# fastVEP 90%-hit CPU-scaling pilot

## Outcome

VCFcache remained materially faster than direct fastVEP at every tested CPU
limit. On the same 4,329,621-variant HG02079 WGS, a 90%-hit cache produced
**2.17x, 2.68x, and 2.56x speedups** with 1, 10, and 32 logical CPUs,
respectively.

| Logical CPUs | Direct fastVEP | VCFcache, 90% hits | Speedup |
|---:|---:|---:|---:|
| 1 | 1,202.3 s (20.0 min) | 553.3 s (9.2 min) | 2.17x |
| 10 | 360.1 s (6.0 min) | 134.5 s (2.2 min) | 2.68x |
| 32 (host maximum) | 326.0 s (5.4 min) | 127.4 s (2.1 min) | 2.56x |

Both paths benefited from parallel execution. From 1 to 10 CPUs, direct
fastVEP became 3.34x faster and VCFcache became 4.11x faster. Increasing the
allocation from 10 to 32 CPUs yielded only another 1.10x and 1.06x,
respectively, so this configuration is already close to its practical scaling
limit at ten CPUs. The slight drop in relative VCFcache benefit at 32 CPUs is
therefore not a loss of caching utility: it is the result of small, unequal
marginal gains once both paths reach their serial and I/O-bound work.

Every CPU-matched cached output was compared with its direct output across all
4,329,621 ordered records after canonicalizing only INFO-tag and CSQ-entry
order. All three comparisons had zero mismatches and the same canonical
SHA-256:
`15c5d795296ec1598237b907f1f3f1af0bf6bb70075ce00bfce66fb90fdacbf6`.

## Design and interpretation boundary

- fastVEP Rayon workers and all VCFcache/bcftools steps received the same CPU
  limit in each pair;
- the cells ran serially on the otherwise idle 32-logical-CPU
  `gvbrowse-preproc` VM in ITCCcloud_dev;
- the cache contained an observed 90.02% of this WGS, leaving 432,201 variants
  for fastVEP;
- VCFcache used lightweight streaming statistics;
- there was one technical run per cell, because this is an engineering scaling
  diagnostic rather than an estimate of biological or timing variance;
- cache construction and untimed output-equality scans were excluded.

The result supports a simple user-facing claim: even around a fast native
annotator, VCFcache cut end-to-end time by about 54-63% at a realistic 90% hit
rate across small and large CPU allocations. It does not establish portable
hardware-independent scaling and should not be presented with error bars.

Machine-readable source data are in `source_data/core_scaling.tsv`. Raw logs,
BCFs, resource reports, the manifest, and equality reports are stored under
`/mnt/data/vcfcache_benchmarks/fastvep_pilot/core_scaling` on ITCCcloud_dev.
