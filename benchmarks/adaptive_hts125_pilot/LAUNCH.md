# Launch record — 7 November 2026

Frozen implementation: `b9633e72dba1a40dd6909200a58f7430a16a1fcf` on
`feature/python-adaptive-hts125-pilot-20261107`. Main is unchanged.

| Job | Block | Worker |
|---|---|---|
| 2801 | vep-rs core WES | sl-w1 |
| 2802 | vep-rs core WGS | sl-w2 |
| 2803 | fastVEP clinical WES | sl-w3 |
| 2804 | fastVEP clinical WGS | sl-w4 |
| 2805 | vep-rs clinical WES, diagnostic | sl-w6 |
| 2806 | After-any collection, no new annotation jobs | Intel workers only |

All five compute jobs were verified running at 18:25 CET. The AMD worker
sl-w5 is excluded. Each timed condition uses its own Slurm step with eight-CPU
affinity. This is a bounded engineering pilot, not the full manuscript cohort.

## Evidence and automatic review

- Controller: `ubuntu@10.133.241.128:/mnt/data/slurm-results/adaptive-hts125-20261107/`
- Worker spelling of the same share: `/results/adaptive-hts125-20261107/`
- Candidate deployment: `/results/deployment/adaptive-hts125-20261107/`
- Unmodified control: `/results/deployment/python-adaptive-reference-8ec84894f/`
- Private matched-codec toolchain: `/results/input-preserving-20261105-tools/`

The frozen manifest/READY checksum, submissions, per-block `state.json`,
commands, resources, complete-record/header comparisons, and logs are retained.
Collection writes `RESULTS.md`, `conditions.tsv`, `REVIEW.json` and an evidence
archive. The bounded controller renderer (initial PID 276694) waits for terminal
review and writes `figures-head-v1/`; `render_state.json` records its outcome.
The renderer expires after three hours and never submits jobs.

The initial estimate is 10 minutes to the first WES block and 30–60 minutes for
all timings and audits. The job wall limit is two hours. If a job fails, inspect
its state and logs; accepted measurements must not be silently overwritten.

## Checks before submission

- Local suite: 456 passed, 18 skipped.
- Focused tests against deployed bcftools/HTSlib 1.25: 11 passed.
- Black, Ruff, whitespace checks and R syntax checks passed.
- Source, input, cache and annotation-asset hashes are checked before timing.
- No performance or output-equivalence result is claimed by this launch record.

See [README.md](README.md) for the seven conditions and interpretation limits.
