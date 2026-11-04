# Python publication-pass fusion: controlled engineering pilot

Completed 5 November 2026 (CET). **Do not promote this experiment to the default.**
Eliminating a separate final decode/recompress pass did not yield a convincing
end-to-end improvement. All 80 timed invocations succeeded, and all 64 cached
outputs passed strict complete-record/relevant-header comparison against their
matched fresh direct output, with no ignored fields. No manuscript measurement
or existing figure has been replaced.

## What changed

The opt-in `VCFCACHE_EXPERIMENT_FUSED_PUBLISH=1` pipes each shard's final overlay
into its post-overlay required-tag filter and final compression, then uses
checked naive BCF concatenation and a final index. It avoids a separately
materialized merge and a whole-file final filter/recompression pass. It does not
remove the overlay, filtering, compression or final index work themselves.
The input/cache formats, normalization and recipe are unchanged; misses are
annotated once per sample, not once per shard. No Rust dependency is introduced.

Production code remains opt-in on `feature/python-pass-fusion-20261104`, based
on `feature/python-adaptive`. The timed source is frozen at
`e72d6fe11fc2df302bea6dd6f7987d0b0fbb08d0`; later commits add analysis, audit
recovery, tests and documentation, not replacement timed measurements.

## Design and measured results

vep-rs **core**, GRCh38, eight-CPU Slurm affinity per condition. The typical WGS
input, KPGP-00139, was repeated on all six workers. The lower-hit input,
LP6005592-DNA_C05, was additionally tested on sl-w2 and sl-w5. Each block has
two technical repeats of five modes: fresh direct annotation, original default
shards (seven here), original four shards, fused four shards and fused two
shards. Condition order rotates across workers and reverses in repeat two.
Strict comparisons are untimed and follow all five timings in each repeat.
OS/filesystem caches were not flushed; the recipe/cache assets remain on shared
NFS and disposable intermediates use worker-local scratch.

Median elapsed seconds across the two technical repeats:

| Input | Worker | Direct | Original default | Original 4 | Fused 4 | Fused 2 |
|---|---|---:|---:|---:|---:|---:|
| Typical | sl-w1 | 84.41 | 98.32 | 101.23 | 106.33 | 97.51 |
| Typical | sl-w2 | 95.70 | 106.24 | 107.19 | 105.42 | 100.35 |
| Typical | sl-w3 | 98.03 | 110.69 | 113.07 | 120.75 | 118.14 |
| Typical | sl-w4 | 83.71 | 98.89 | 103.36 | 110.83 | 97.38 |
| Typical | sl-w5 | 51.34 | 58.09 | 60.54 | 63.50 | 59.64 |
| Typical | sl-w6 | 85.98 | 99.28 | 101.12 | 102.41 | 97.67 |
| Lower-hit | sl-w2 | 128.10 | 137.85 | 144.76 | 143.11 | 137.39 |
| Lower-hit | sl-w5 | 67.61 | 80.63 | 83.09 | 85.54 | 80.76 |

Typical input: 4,634,658 records; 418,884 misses (**90.96% hits**).
Lower-hit input: 6,169,350 records; 1,296,702 misses (**78.98% hits**).
Every cached condition retained the expected complete record count.

For the typical input, first form each repeat's original-default/candidate
ratio, take its two-repeat median within each node, then summarize across the
six nodes:

| Candidate | Median gain over original default | Node range | Faster nodes |
|---|---:|---:|---:|
| Original 4 shards | 0.976× | 0.957–0.991× | 0/6 |
| Fused 4 shards | 0.921× | 0.898–1.008× | 1/6 |
| Fused 2 shards | 1.012× | 0.939–1.060× | 4/6 |

For the lower-hit input, fused two-shard gains over default are 1.004× on sl-w2
and 0.998× on sl-w5. All node/mode median speedups against fresh direct
annotation remain below 1.0×. This is an engineering test of two genomes and
two technical repeats, not a new biological cohort estimate or a significance
test; a roughly 1% typical median improvement is not a reliable product gain.

## Interpretation and recommendation

The earlier 1.06×/1.19× proposals were optimistic **counterfactual stage
subtractions**, not measured forecasts. This implementation relocates work;
filtering and final compression still cost time, and the extra main process per
pipe changes how a fixed CPU budget is shared. Four publication pipes receive
no additional compression threads at eight CPUs; two pipes allow two codec
threads each. Fewer shards also change cache-overlay and miss-selection costs.
Consequently, removing a materialization does not remove its entire old wall
time. The measurements do not isolate one universal bottleneck or rule out a
more substantial streaming design.

Keep the experimental flag off and **do not rerun the full manuscript cohort
on this candidate**. Retain the existing Python-adaptive baseline and separate
Rust results. If further performance work is approved, investigate end-to-end
merge/index/compression costs and a true streaming overlay under the same
fixed-CPU controls; do not invest first in loading the whole genome cache into
RAM on the strength of the old stage-subtraction model. Cohort miss batching
remains a separate proposal, not a feature tested or implemented here.

## Infrastructure lesson for the symposium

The same fresh-direct workload took 51.34 s on sl-w5 versus 83.71–98.03 s on the
other workers: up to **1.91×** difference at the same allocated CPU count.
sl-w5 exposes AMD EPYC-Rome; the others expose Intel Broadwell. Both cached
and direct annotation benefit on sl-w5. This supports a VM/node-performance
warning, but **does not isolate disks as the cause**, prove physical-host CPU
specifications, or measure cBioPortal performance. Use the slide-ready wording
in [SYMPOSIUM_SLIDE.md](../SYMPOSIUM_SLIDE.md).

## Evidence, recovery and figures

- `conditions.tsv`: all 80 timings, counts, CPU models and resources.
- `paired_ratios.tsv` and `summary_by_node.tsv`: within-repeat ratios and
  within-node summaries; no fastest-repeat selection.
- `ACCEPTANCE.json`: eight accepted blocks, including one explicit audit-only
  recovery; `attempts.tsv` retains the distinction.
- `terminal-evidence.tar.gz`: frozen compact raw evidence, including commands,
  recipe/asset bindings, hardware metadata, comparisons and logs. `SNAPSHOT.json`
  binds its SHA-256 and the original manifest; `EVIDENCE_INVENTORY.json` provides
  file hashes. Extract into a new directory to reproduce the tables.
- `SUPERSEDED_LAUNCH.json`: the initial runner failure, retained separately and
  excluded from accepted results. Its `adaptive --uncached` argument mistake
  was fixed before the v2 campaign; these are not annotation-engine failures.
- Jobs 2782–2789 produced all timings. Job 2788 reached its one-hour limit after
  all ten timings, during untimed audits. Recovery job 2792 compared the two
  remaining unchanged retained outputs, preserving the original state and all
  timings. Its receipt is inside the evidence archive. Final analysis job 2793
  accepted all eight blocks. No further jobs were queued at final review.
- Safety verification: **472 tests passed, 18 skipped**; targeted Ruff and Black
  checks passed. Tests cover post-overlay filtering, empty inputs, all hits,
  preservation, incompatible headers, pipeline failure and binary stdout.

[Figure 1: complete pilot](figures-v1/Figure1_same_node_pass_fusion.pdf).
Elapsed times by condition, with separate points/lines for the two technical
repeats and separate facets for each node/input. Facet y scales differ; use the
table above for absolute cross-node comparisons. Individual panels, SVG and
600-dpi PNG exports, the named tibble-list RDS and R session information are in
`figures-v1/`. `figures-preliminary-v1/` is the earlier seven-block snapshot and
is preserved, not overwritten; it must not be used as the complete result.
