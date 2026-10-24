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

Preparation and launcher completed successfully. Slurm accepted these jobs at
21:18 CEST; all six groups have a six-hour allocation limit:

| Job | Group | Sample | Recipe | Node |
| --- | --- | --- | --- | --- |
| 1335 | g00 | LP6005519-DNA_A06 | fastVEP clinical | sl-w2 |
| 1336 | g01 | LP6005519-DNA_A06 | vep-rs core | sl-w4 |
| 1337 | g02 | LP6005519-DNA_A06 | VEP everything | sl-w1 |
| 1338 | g03 | LP6005592-DNA_C05 | fastVEP clinical | sl-w3 |
| 1339 | g04 | LP6005592-DNA_C05 | vep-rs core | sl-w6 |
| 1340 | g05 | LP6005592-DNA_C05 | VEP everything | sl-w2 |

Collector **1341** depends on all six with `afterany`. After submission its
`ExcNodeList` was explicitly set to `sl-w5`, matching the worker exclusion in the
design. This scheduler-only correction does not change the frozen executable
snapshot; later source includes the same exclusion when submitting collectors.
At initial acceptance five groups were running and g05 was waiting for sl-w2.

Frozen fixture counts: A06 **82,714 WES / 4,858,092 WGS** records; C05
**108,539 WES / 6,169,350 WGS** records. All four are verified biallelic.
The [manifest and accepted submissions](campaign-20261023/) are archived locally;
the manifest SHA in `READY` protects the source/resource fingerprints.

## Initial progress snapshot (not final results)

The final pre-handoff collection recorded **12 complete, 5 running, 69 not yet
started**, with no failed cells. Four instrumented vep-rs core WES cells
(`g01-c002`, `g01-c003`, `g04-c002`, `g04-c003`) passed both strict direct and
uninstrumented-cache comparisons. Their Python/Rust spans and independent Slurm
cgroup counters were verified. No duplicate cgroup paths or trace parse errors
were reported. The latest archived report is
[6ce7d5758814fcfc](campaign-20261023/reports/6ce7d5758814fcfc/REPORT.md).

The local campaign directory is an explicitly partial metadata snapshot: live
jobs continue remotely. Large inputs, output BCFs and perf binaries are excluded
from git; reference/failed outputs remain in the shared campaign directory.
Recollect remotely for current progress or the final optimization assessment.

## Check status or collect later

On the controller:

```sh
cd /tmp
sudo -u appuser env \
  PYTHONPATH=/mnt/data/slurm-results/deployment/engine-profile-recovery-20261025-183514646 \
  /mnt/data/slurm-results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python \
  -m benchmarks.engine_profile.run status \
  /mnt/data/slurm-results/engine-profile/20261023T190736Z-97545d739
```

Replace `status` with `collect` to regenerate the content-addressed report.
The review command above uses the corrected 25 October harness; the original
engine deployments and manifest remain frozen. See [RECOVERY-20261025.md](RECOVERY-20261025.md)
for the separately submitted VEP WGS memory recovery.
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
