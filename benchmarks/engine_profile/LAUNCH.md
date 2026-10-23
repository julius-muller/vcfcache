# Engineering profiling campaign — 23 October 2026

Branch: `feature/engine-overhead-profile`.
Frozen profiling snapshot: `97545d739`; baseline: `296855d72`.
The optional helpers both report 0.4.0; their distinct binary SHA256 hashes and
source revisions are recorded in [launch-20261023.json](launch-20261023.json).

Campaign on workers:
`/results/engine-profile/20261023T190736Z-97545d739`.
Controller (`ubuntu@10.133.241.128`):
`/mnt/data/slurm-results/engine-profile/20261023T190736Z-97545d739`.

Preparation job: **1333**. This job stages and hashes immutable fixtures before
the 86-cell matrix is submitted. Submission acceptance and live progress are
recorded separately; preparation is not evidence that timed cells succeeded.
Launcher job: **1334**, dependency `afterok:1333`, with invalid dependencies
cancelled. It submits the six groups and their automatic `afterany` collector
only after successful preparation, without requiring this session to survive.

## Check status or collect later

On the controller:

```sh
cd /tmp
sudo -u appuser env \
  PYTHONPATH=/mnt/data/slurm-results/deployment/engine-profile-20261023T190736Z-97545d739-profile \
  /mnt/data/slurm-results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python \
  -m benchmarks.engine_profile.run status \
  /mnt/data/slurm-results/engine-profile/20261023T190736Z-97545d739
```

Replace `status` with `collect` to regenerate the content-addressed report.
`submissions.json` stores accepted jobs; `latest-review.json` points to the latest
report. Every attempt keeps its commands, resource counters, raw comparisons,
fallback details and traces. Inspect failures rather than treating an empty queue
as proof of success. See [README.md](README.md) for explicit retry semantics.

## Verification before launch

- Host full suite: **556 passed, 21 skipped**.
- Container full suite with the deployment helper: **556 passed, 21 skipped**.
- Rust: **2 tests passed**; clippy with warnings denied passed.
- Python lint and formatting checks passed.
- [GitHub CI for the frozen profiling commit](https://github.com/julius-muller/vcfcache/actions/runs/40009158498): all required jobs passed.
- Real-BCF tests confirm profiling-on/off record equivalence for both backends
  and preserve traces of a native memory-budget fallback.
- Harness tests cover the 86-cell bound, interval unions, partial submission with
  an `afterany` collector, offline/read-only status, idempotent collection, and
  manifest-drift rejection.
- R/ggplot2 smoke test generated both PDF and 400-dpi PNG outputs. On this
  workstation run R with `--vanilla` and
  `R_LIBS_USER=/home/j380r/.local/lib/R/4.5`; its default startup profile points
  to an unavailable RStudio directory. Worker collection does not require R.

No existing publication results, figures, cache files, or default-engine settings
were changed. This is a profiling campaign, not a new validated speedup claim.
