# First completed block — 7 November 2026

This is a **partial engineering result**, not the final pilot report. Remaining
blocks continue automatically; see the controller location in [LAUNCH.md](LAUNCH.md).

KPGP-00139, GRCh38, vep-rs core WES; job 2801, Intel sl-w1, eight-CPU affinity.
Input: 75,946 records; 7,364 misses in every cached condition (90.30% hits).
One timing per condition, no best-of-repeat selection. All six cached outputs
pass strict complete-record/relevant-header comparison against fresh direct
annotation, with **zero mismatches and zero ignored fields**.

| Condition | Elapsed seconds |
|---|---:|
| Original Python-adaptive + 1.21 | 19.25 |
| Same source + 1.25 | 20.04 |
| Prepared input + 1.25, original auto shards | 19.07 |
| Prepared input, contigs, 2 shards | 18.65 |
| Prepared input, merged windows, 2 shards | 15.61 |
| Prepared input, merged windows, 4 shards | 17.50 |
| Prepared direct annotation + 1.25 | 8.99 |

Interpretation within this block:

- The toolchain upgrade alone did not help this original whole-contig route.
- Two-shard windows reduced cached elapsed time by approximately 19% versus
  the original 1.21 control, or 18% versus prepared auto on 1.25.
- At fixed two-shard budgeting, windows took 16% less time than whole contigs.
- Four shards did not beat two; more subprocesses/codec workers are not
  automatically faster under a fixed CPU allocation.
- The best cached observation is still only 0.58-fold direct throughput
  (15.61 versus 8.99 seconds): **no cached vep-rs core WES win**.

These observations do not establish statistical confidence or apply to other
recipes/samples. All conditions share a worker but not an isolated storage
backend; order and shared I/O can affect these single measurements. No code is
promoted to main based on this partial result. vep-rs clinical remains a
diagnostic block because of the existing cache/fresh annotation discrepancy.

The immutable campaign definition and this completed block's commands,
resources, logs, statistics, and comparisons are copied into
`results-20261107/`. Final collection/figures remain automatic on the controller.
