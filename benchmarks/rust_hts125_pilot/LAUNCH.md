# Launch — 7 November 2026, 23:15 CET

Frozen source: `533dcd36022001e9b364b42fca9f2474464cb9ed`, branch
`feature/rust-hts125-pilot-20261107`. Main and previous measurements are unchanged.

| Job | Tool / assay / input | Worker |
|---|---|---|
| 2807 | vep-rs core WES, KPGP-00139 | sl-w1 |
| 2808 | fastVEP clinical WES, KPGP-00139 | sl-w2 |
| 2809 | vep-rs core WES, LP6005592-DNA_C05 | sl-w3 |
| 2810 | fastVEP clinical WES, LP6005592-DNA_C05 | sl-w4 |
| 2811 | vep-rs core WGS, KPGP-00139 | sl-w6 |
| 2812 | fastVEP clinical WGS, KPGP-00139 | sl-w1, after job 2807 |
| 2813 | After-any collection; never submits annotation | Intel workers only |

At 23:15 CET the first five blocks were verified running. sl-w5 is explicitly
excluded. Every timed CLI has a separate Slurm step with eight-CPU affinity.
The queued sixth block intentionally reuses sl-w1 after the first block ends.

## Review locations

- Controller: `ubuntu@10.133.241.128:/mnt/data/slurm-results/rust-hts125-20261107/`
- Worker spelling of the same share: `/results/rust-hts125-20261107/`
- Source deployment: `/results/deployment/rust-hts125-20261107/`
- Private helper builds/provenance: `/results/rust-hts125-20261107-tools/`
- External 1.25 tools: `/results/input-preserving-20261105-tools/`

`manifest.json`/`READY`, submissions, per-block state, commands, resources,
statistics, hashes and strict comparisons are durable. Collection writes
`RESULTS.md`, `conditions.tsv`, `REVIEW.json` and `terminal-evidence.tar.gz`.
It labels interrupted checkpoints honestly if a scheduler timeout occurs.
The bounded controller R watcher (initial PID 281639) renders new PDF/SVG/PNG
figures and individual panels into `figures-head-v1/`; its result appears in
`render_state.json`. It expires after five hours and submits no jobs.

Expected full pilot/audit duration: roughly 30–60 minutes from launch, with a
two-hour safety wall limit per block. This is an estimate, not a success claim.
No 12-WES validation or full manuscript cohort is queued automatically.

## Build / preflight evidence

Both helpers have the same Rust sources/dependency versions and reported codec
features, including libdeflate. Runtime `build-info` confirms **HTSlib 1.19.1**
versus **1.25**; both execute on the controller. Source/ABI bindings are rebuilt,
not swapped at runtime. Binary SHA-256 values are frozen in the manifest.

The final full local suite passed **724 tests (21 skipped)**; all 13 focused
tests also passed on the controller with each helper/toolchain.
The old helper also passed 123 native/adaptive
regression tests locally (3 skipped). Rust unit tests, formatting, lint,
whitespace and synthetic R rendering checks passed. No benchmark performance
or output-equivalence claim is made by this launch record.

See [README.md](README.md) for conditions, scope and interpretation limits.
