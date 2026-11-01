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

Submission receipts and exact deployment/campaign paths will be recorded below
after Slurm accepts the launch jobs.
