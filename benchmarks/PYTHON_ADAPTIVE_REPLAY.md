# Python-adaptive manuscript replay

Decision: ship the Python/bcftools adaptive CLI default, retain legacy as an
explicit compatibility option, and keep historical Rust experiments separate.
No compiled helper is included in this deployment. Historical Rust measurements
are not Python measurements and must not be relabelled.

## Safety gate

The execution wrapper now actually runs its recorded `set -euo pipefail` script.
Both engines, direct annotation and cache builds reject unexpected recipe-boundary
record loss/expansion and unreadable counts. Tests exercise failing pipelines and
zero-exit, indexed but truncated output. Recipe counts are taken **after** input
normalization; they supplement, not replace, untimed semantic comparison.

Job 2043 must not be described as an accepted truncated campaign cell: its outer
runner failed on WGS direct annotation. The application-level shell gap is real,
but the outer failure and subsequent replacement measurements remain distinct.

## Replay and urgency

`python_adaptive_replay.py` freezes the original 600-condition manuscript matrix,
including repaired task definitions, under a new result directory. It binds the
252 retained paired direct controls by command, CPU allocation, resources and
artifact SHA256. Other historical direct references retain their documented
validation scope. No direct annotation is scheduled implicitly.

Six first conditions gate automatic dispatch; failed checks stop release of the
remainder. Each condition has an exclusive Slurm job and independent measured
step. All six eligible workers are used. Main data are queued before the 12
synthetic-delay supplement conditions, with longer conditions first within each
priority. Jobs checkpoint results and run an afterany collector without a live
assistant session. Known clinical vep-rs differences and SnpEff exclusions are
not silently accepted as validated speedups.

The prior cached measurements total roughly 96 worker-hours (including about
20 hours of artificial-delay conditions): six workers imply approximately 16
hours of annotation, **before** semantic auditing and Python-engine differences.
This is an estimate, not a promised completion time; early primary data should
arrive much sooner. Do not sacrifice matched CPU allocations or validation to
shorten the estimate. Old results and printed figures are untouched.

Benchmark-only runner helpers were restored from performance-sprint commit
`2bee82a4d`; they do not add Rust code or dependencies to the product. The new
freeze records the exact deployed product and runner files. Review with:

```bash
PYTHONPATH=DEPLOYMENT PYTHON -m benchmarks.performance_sprint.manuscript status CAMPAIGN
```

## Accepted launch (2026-11-01, 22:00 CET)

- Product freeze: `ee2162473b50bd57eb83c0495ad63bfc22deaf83`.
- Controller deployment: `/mnt/data/slurm-results/deployment/python-adaptive-20261101-ee2162473`.
- Controller campaign: `/mnt/data/slurm-results/python-adaptive-manuscript-20261101-ee2162473/campaign`.
- Worker paths use `/results/` in place of `/mnt/data/slurm-results/`.
- Launch cells: **2070–2075**, one on each of sl-w1 through sl-w6.
- First collector: **2076**; automatic continuation: **2077**.
- Local immutable submission receipts: [python_adaptive_replay_20261101](python_adaptive_replay_20261101/).

Five launch conditions passed by 22:02 CET; the sixth (fastVEP WGS) was still
running normally. The initial idle workers are intentional while that last
gate completes; the automatic continuation then fills all eligible nodes.

**ETA correction:** 16 hours above is annotation time alone, not end-to-end
completion. The previous campaign occupied about **191 worker-hours after
subtracting direct annotation**, including preparation and semantic audits:
approximately **32 hours on six workers**. Python performance may change this.
Use roughly **30–40 hours for the entire audited 600-cell matrix**, not a promise
of overnight completion. Partial validated results are already arriving;
artificial-delay supplement cells remain last. This distinction must accompany
any timeline reported to the user.
