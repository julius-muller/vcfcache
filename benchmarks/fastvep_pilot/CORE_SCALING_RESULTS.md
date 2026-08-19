# fastVEP CPU availability and annotation-complexity result

## Outcome

VCFcache remained faster than direct fastVEP at every enforced CPU limit. On
the held-out 4,795,706-variant KPGP-00319 WGS, a controlled 90%-hit cache
produced **1.51x, 2.57x and 2.55x speedups** for the core consequence recipe at
1, 10 and 32 CPUs, respectively.

| CPUs available to process | Direct fastVEP | VCFcache, 90% hits | Speedup |
|---:|---:|---:|---:|
| 1 | 868.4 s (14.5 min) | 576.5 s (9.6 min) | 1.51x |
| 10 | 463.8 s (7.7 min) | 180.3 s (3.0 min) | 2.57x |
| 32 | 379.1 s (6.3 min) | 148.8 s (2.5 min) | 2.55x |

Both paths benefited from parallel execution, but marginal gains diminished.
From 10 to 32 CPUs, direct runtime improved only 1.22-fold and cached runtime
1.21-fold. The slight change in relative cache benefit at 32 CPUs is therefore
not a loss of caching utility; both paths are approaching their serial and
I/O-bound work.

Every CPU-matched cached output was compared with its direct counterpart using
the strict fastVEP complete-record comparator. All INFO and FORMAT values and
relevant header definitions were checked after canonicalizing only INFO-tag
and CSQ-entry ordering. All comparisons passed.

## Annotation-complexity interaction

The same genome and cache design were tested with a deliberately dense recipe
containing ten supplementary databases. At 10 CPUs, direct and cached runtime
were 706.2 and 230.0 seconds (3.07x); at 32 CPUs, they were 691.4 and 205.4
seconds (3.37x). The 10-CPU dense direct run was 18.7% faster than the one-CPU
core run, confirming that extra hardware can offset annotation complexity.
However, it remained 52.2% slower than the core run given the same 10 CPUs,
and VCFcache returned 7.9 minutes per genome. Increasing the dense allocation
from 10 to 32 CPUs changed direct runtime by only 2.1%, while the cache benefit
persisted.

## Design and interpretation boundary

- fastVEP Rayon workers and VCFcache/bcftools received the same CPU allocation
  in each pair, enforced process-wide using `taskset`;
- the cells ran serially on the otherwise idle 32-logical-CPU
  `gvbrowse-preproc` VM in ITCCcloud_dev;
- the controlled cache contained 4,316,659 records, or 90.01% of this WGS;
- there was one technical run per cell, because this is an engineering scaling
  diagnostic rather than an estimate of biological or timing variance;
- cache construction and untimed output-equality scans were excluded.

The result is a controlled interaction experiment, not portable
hardware-independent scaling, and it is presented without inferential error
bars. The biological hit-rate estimates remain those of the independent
52-genome campaign; the controlled 90% cache here isolates computing effects.

Machine-readable source data and frozen manifests are in
`source_data/core_scaling_kpgp00319*` and
`source_data/heavy_core_scaling_kpgp00319*`. Raw logs, BCFs, resource reports
and equality reports are stored under the corresponding run directories below
`/mnt/data/vcfcache_benchmarks/fastvep_pilot/` on ITCCcloud_dev.
