# Python publication-pass pilot

An isolated experiment based on `feature/python-adaptive`; no Rust routing
helper, new cache format, assay-specific dispatch or product-default change.

## Implementation

Opt in with `VCFCACHE_EXPERIMENT_FUSED_PUBLISH=1`. The default remains the original
Python-adaptive algorithm. With at least two configured CPUs the experiment:

1. Keeps normalization, cache overlay, miss selection, and one annotator
   invocation under the existing contracts.
2. Pipes each shard's final annotation overlay as uncompressed BCF into a
   **post-overlay** filter and final-level compressor. `annotate -i` would filter
   before overlay and lose newly annotated misses, so it is not used.
3. Concatenates final compressed shards with checked `--naive`, never
   `--naive-force`, and writes the final CSI once. No full-size `step4.bcf`,
   gathered intermediate followed by a final decode/recompress pass, or new
   annotation-store conversion is required.

At most `CPUs // 2` piped shards run concurrently. Each pipe budgets two main
processes plus its output-compression threads; non-piped shard phases budget
one main process plus compression threads. Eight CPUs therefore allow four
two-process publication pipes without extra codec threads, or two pipes each
with two extra codec threads. The one-CPU path retains original publication.

Zero misses, empty inputs, `--preserve-unannotated`, output indexing and binary
stdout are covered. Compatible BCF dictionaries are checked by bcftools naive
concatenation; a mismatch or producer failure must not publish the output.
Publication is staged on the destination filesystem and count-checked before
replacing the requested BCF/index. Original intermediate files and the original
algorithm remain available when the experiment is disabled.

## Controlled design

- Frozen original recipes/assets, GRCh38, **vep-rs core WGS**, eight CPUs per
  independently timed Slurm step. No normalization advantage is granted to either
  cached or direct annotation. Strict complete-record and relevant-header
  comparison runs outside timing, without field exceptions.
- Typical input **KPGP-00139** on **all six workers**, explicitly pinned. This
  gives the same input on each guest CPU/storage configuration rather than
  inferring a node effect from different genomes.
- Lower-hit **LP6005592-DNA_C05** follows on `sl-w2` and `sl-w5`, conditional on
  the corresponding typical block passing. It remains a real stress case, not
  an excluded outlier. The available data give approximately 91% versus 79% WGS
  hits, respectively; current measurements record their actual counts.
- Every block runs two repeats of: **fresh direct**, **original default shard
  count**, **original four shards**, **fused four shards**, **fused two shards**.
  Mode order rotates across blocks and reverses for repeat two. These few
  technical repeats estimate engineering variability; they are separate from the
  publication's independent-sample/no-technical-repeat campaign.
- Input preparation, asset verification, hardware capture and semantic audits
  are untimed. Processes start fresh; OS/filesystem caches are **not flushed**.
  Comparisons occur after all timed conditions in a repeat, avoiding comparison
  reads between arms. Large disposable files are confined to owned local scratch.
- Capture guest `lscpu`, kernel, local mount and free-space metadata per block.
  Guest CPU labels are not a complete physical-hypervisor inventory. Resource
  samples contain cgroup-v2 CPU, memory and local block-I/O counters; memory
  includes page cache, and local block-I/O counters do not cover all NFS traffic.

There are eight blocks / **80 planned timed invocations**, not a new full-cohort
campaign. Only one AMD-labelled worker is available in this design; the two
lower-hit blocks serialize behind their same-node typical blocks. Idle nodes
after their assigned work finishes are intentional, not missing scheduling.

## Execution and later review

Worker paths use `/results`; the controller exposes the same shared tree under
`/mnt/data/slurm-results`.

The accepted retry definition is:

- Root: `/results/python-pass-fusion-20261104-v2`.
- Source: `e72d6fe11fc2df302bea6dd6f7987d0b0fbb08d0`.
- Typical jobs: **2782–2787**; lower-hit jobs: **2788–2789**.
- Automatic after-any review: **2790**.
- Each block has a one-hour Slurm limit, a fresh eight-CPU step per condition,
  two same-node replicates, and durable JSON checkpoints.

The first launch, `/results/python-pass-fusion-20261104`, is explicitly
superseded. Its runner incorrectly combined `--engine adaptive` with `--uncached`,
which this Python branch rejects. Partial/failed evidence is retained and
`ABORTED.json` records the reason; jobs 2773–2781 are not accepted benchmark
results. The retry uses a separate immutable deployment and output root. The
CLI builder and review collector now have explicit regression tests.

Review never submits another campaign:

```sh
python -m benchmarks.pass_fusion_pilot.run collect /results/python-pass-fusion-20261104-v2
```

`manifest.json` plus `READY` bind the definition; `submissions.json` preserves
commands and IDs. Each case's `state.json` contains its node, completed timings
and terminal status. Each condition retains resources, command/task/parameter
bindings, statistics and CLI logs; completed repeats retain strict comparisons.
Failed/partial conditions are never promoted or replaced with their fastest
retry. Slurm accounting may be disabled; the checkpoints remain reviewable.

No publication figure or manuscript number is changed by this pilot. Promotion
requires consistent same-node benefit and complete-record equivalence, not the
optimistic subtraction of entire old stages from elapsed time.
