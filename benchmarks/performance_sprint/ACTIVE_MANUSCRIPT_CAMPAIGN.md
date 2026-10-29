# Active manuscript campaign

This is a durable launch record, **not live status**. Refresh Slurm/checkpoints
before reporting progress.

## Final collection and analysis — 29 October

The main and repair campaigns and final collectors **2000, 2039 and 2040
completed**. Last condition finished at **12:57 UTC / 13:57 CET**. The later
live queue check found no jobs and six idle workers; no automatic annotation
stage remains. The launch/recovery sections below are preserved history.

- **600** canonical terminal conditions; **604** execution attempts retained.
- **572** accepted under their declared comparator scopes, **24** clinical
  vep-rs timing-only, **four** known SnpEff strict warning differences excluded.
- All **36** repaired clinical SnpEff GRCh37 slots finished: 34 accepted, two
  known content discrepancies. This is not an outstanding runtime/header failure.
- **12** aliases restore 612 plotting views without increasing independent
  measurements or donor counts. Fifty-four distinct genome IDs overall.
- Commands, comparisons, states, recipes, resource metrics and source hashes
  were copied read-only into the new
  [final analysis folder](final-analysis-20261029/README.md); large BCFs remain
  at their original remote paths. Existing figures were not overwritten.
- Read the [decision brief](final-analysis-20261029/DECISIONS.md) before changing
  defaults or authorizing another benchmark. No new jobs were submitted by the
  final consolidation, and no further work advances automatically after it.

## Targeted header repair, 28 October

The user approved the INFO-header compatibility repair after three clinical
SnpEff GRCh37 WGS cached conditions failed (`Integer/.` versus `Integer/A`).
The original controls, caches, candidate freeze and failed attempts are intact.
See [repair details](HEADER_COMPATIBILITY_REPAIR.md).

- Corrected build: **787a3081a**, helper **0.7.1**. Both host and Docker suites
  pass **666 tests, 21 skipped**. All other conditions retain their frozen build.
- Repair campaign: `/results/manuscript-final/20261028T210624Z-787a3081a/campaign`.
  Receipt: `header-repair-launch-20261028.json`.
- Gate **2001 passed**: 165 records, zero full-record mismatches, matching relevant
  headers, no ignored fields. This was an unstarted final condition, not a repeat
  of a successful sample/tool/configuration.
- **35 remaining replacements accepted, jobs 2004–2038**; collector **2039**.
  Final parent collector **2040** depends on both original collector 2000 and
  repair collector 2039. Metadata-only gate jobs 2002/2003 were superseded by
  controller-side collection/dispatch after the gate passed, avoiding exclusive
  worker allocation just for queue bookkeeping. `CONTROLLER_DISPATCH.json`
  records this; it is not a failed benchmark.
- Original pending counterparts were held and then cancelled only after their
  replacements were accepted. No running original was killed and no successful
  condition was repeated. Only repair jobs' queue priority was raised; CPU,
  memory, annotation configuration and node-exclusivity settings are unchanged.
- Three retained WGS direct controls passed exact reuse preflight, preserving
  **4.8 node-hours** of completed annotation. A fourth original control (1568)
  is allowed to finish; replacement 2007 has `afterany:1568` and freezes/reuses
  its successful direct artifacts before cached timing.
- Replacement attempts are **a003/a004 in the original checkpoint tree** and
  their accepted job IDs are in the original submission history. The original
  600 logical conditions therefore remain the authoritative inventory. Historical
  failed a002 attempts stay archived; until a retry starts, a raw latest-state
  scan may still show its old failure. Check its replacement job ID in the
  submission history rather than mistaking an accepted pending retry for no action.
- At 21:20 UTC / 22:20 CET, all six workers were active and first WGS retry 2004
  was running. Repair execution and final collection need no further user input.
  This snapshot does not claim that all 36 replacements have already passed.

The main status command below still reads all 600 slots, including replacement
attempts. For repair-only checkpoints, use the repair deployment and campaign.

- Candidate: `fe197f34d` (same product files as the previously validated checks).
- Worker campaign: `/results/manuscript-final/20261027T191111Z-fe197f34d/campaign`.
- Deployment: `/results/deployment/manuscript-final-20261027T191111Z-fe197f34d`.
- Controller mapping: replace `/results/` with `/mnt/data/slurm-results/`.
- Receipt: `manuscript-launch-resumed-20261027.json`.
- Preparation job **1398 completed**. It froze all 600 conditions and submitted
  the continuation. Primary jobs **1400–1403** started on four workers; collector
  **1404** and dispatcher **1405** were registered with dependencies. The VEP and
  vep-rs checks were carried forward with identical product/configuration hashes.
  **All six primary conditions passed**. At **19:33 UTC on 27 October**, automatic
  bulk dispatch was verified: **594 jobs (1406–1999)** accepted, collector **2000**
  registered. All six workers were running jobs; 588 annotation jobs plus the
  collector were pending. Checkpoints showed six completed primary conditions,
  six running bulk conditions and **zero failures**. This is a launch snapshot,
  not a claim that the full campaign has completed.
- Matrix: **600 logical cached conditions**, preserving the original cohorts and
  experiment families, with 12 duplicate historical rows represented as aliases.
  There are 72 pinned vep-rs extension pairs; the 180 clinical projection-only
  references require fresh direct controls. Compatible historical direct results
  are reused. No technical replicate averaging or best-run selection is planned.
- Primary engine: adaptive, 8 allocated/bound CPUs, 28 GB, exclusive nodes.
  All six healthy workers are eligible. Native remains opt-in.
- Deadline: Monday **2 November 2026, 18:00 CET**. Initial provisional runtime
  estimate is **30–48 hours after bulk release**, including new direct controls
  and audits. Re-estimate from measured progress, not from this file's age.
- Verification: **646 passed, 21 skipped** in both the host and Docker suites;
  Ruff and whitespace checks pass. GitHub CI passed deployed revision `fe197f34d`.

## Refresh later

Run on the controller as appuser, with mapped paths:

```sh
cd /tmp
sudo -u appuser env \
  PYTHONPATH=/mnt/data/slurm-results/deployment/manuscript-final-20261027T191111Z-fe197f34d \
  /mnt/data/slurm-results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python \
  -m benchmarks.performance_sprint.manuscript status \
  /mnt/data/slurm-results/manuscript-final/20261027T191111Z-fe197f34d/campaign
squeue -o '%.12i %.10T %.10M %.20N'
```

Before `campaign/READY` exists, inspect the preparation log:
`/mnt/data/slurm-results/manuscript-final-prepare-20261027T191111Z-fe197f34d-1398.out`.
After preparation, inspect `START.json`, `CARRIED_FORWARD.json`, `submissions.json`
and `AUTO_DISPATCH.json`. `DISPATCH_BLOCKED.json` means the primary checks failed
and active repair is required; it must never be described as progressing.

The primary gate uses actual final conditions, not extra replicates. Unchanged
successful checks can be carried forward only after exact product/configuration
matching. Once the gate passes, all remaining ready conditions are submitted
automatically, with an afterany collector. Slurm and collectors survive disconnect;
there is no autonomous agent repairing unexpected failures after the turn ends.

## Validation and recovery boundaries

Full records/headers are compared outside timing. Known clinical vep-rs CSQ
differences remain timing-only, with both outputs retained. Unexpected differences
fail; no broad HGNC_ID exception applies to non-VEP tools. Synthetic delay cells
retain an explicitly composed vanilla-VEP reference, not invented per-delay BCFs.

Cache payloads, indexes, recipes, parameters and provenance are frozen together.
Drifted worker replicas are repaired in private node-local storage, untimed, while
original caches are preserved. Node-local repair was executed successfully on
sl-w2 in setup probe 1399 (7,967,504 expected records).

Historical clinical SnpEff payloads have a non-ASCII text-representation limitation.
Both new arms explicitly use compatible ASCII text encoding and `-noId`. The
corrected 242-record direct control passed with no ignored fields. This does not
establish Unicode fidelity to source databases; see `launch-diagnostics-20261027/`.

Preparations 1376/1397 failed safely before bulk release; 1395 was deliberately
cancelled during preparation. The prior first-check campaign also remains intact.
Verified payloads are reused after rehashing, and delayed filesystem metadata
changes trigger bounded fresh-hash retries, never a content-hash waiver.

## Handoff

No new scientific decision is currently required. The next useful check-in is
approximately **21:35 UTC / 22:35 CET on 27 October** (two hours after bulk release).
Inspect failures and six-worker
utilization then; do not resume optimization experiments or repeat successes.
