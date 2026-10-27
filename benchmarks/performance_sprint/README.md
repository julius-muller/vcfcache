# Deadline-bounded performance sprint

See [the sprint handoff](../../PERFORMANCE_SPRINT.md) for the deadline, gates and
scope. These candidates do not change default annotation behavior.

## Reproducible first screening

Build the compatible helper, run tests and commit the candidate before deployment:

```sh
docker build -f docker/Dockerfile.vcfcache --target native-build -t vcfcache:sprint-helper .
.venv/bin/python -m benchmarks.performance_sprint.deploy \
  --helper COMPATIBLE_HELPER --receipt UNIQUE_RECEIPT.json
```

The deployment freezes source, helper SHA/version and 14 cached-only cells;
submitters write job IDs immediately and register an afterany collector.
Original direct/cached output references are reused and checksum-verified.
Baseline and candidate profile the same coarse stages, on pinned nodes with
16-CPU affinity. Sample B reverses the sample A order. No benchmark annotation
resource is patched or downloaded.

On the controller use mapped deployment/campaign paths and the shared Python:

```sh
PYTHONPATH=DEPLOYMENT PYTHON -m benchmarks.engine_profile.run status CAMPAIGN
PYTHONPATH=DEPLOYMENT PYTHON -m benchmarks.performance_sprint.pilot report CAMPAIGN
PYTHONPATH=DEPLOYMENT PYTHON -m benchmarks.engine_profile.run collect CAMPAIGN
```

`pilot report` reports single-run changes, **not automatic promotion**. Candidate
winners need focused paired confirmation within the 36-cell budget. Inherited
clinical discrepancies never justify new engine-output differences.

Do not submit the full publication rerun before baseline-reuse decisions and
the release-candidate freeze are recorded. Existing cached result figures are
historical comparisons, not replacement evidence for a new engine.

## Review days later

Launch receipts bind each campaign, source/helper identity and accepted Slurm IDs.
`deploy` locks the canonical receipt directory and counts uncertain submissions
conservatively against the **36-cell total**; selecting another receipt path
does not reset the budget.

```sh
python -m benchmarks.performance_sprint.review NEW_LOCAL_SNAPSHOT
python -m benchmarks.performance_sprint.references prepare NEW_REFERENCE_REQUEST.json
```

The first command reads remote small evidence only; no annotation calls, retries
or source changes. A snapshot may contain running states and is not live status.
Large BCFs stay on the cluster. Each snapshot contains checksums and `cells.tsv`.

`evidence-final-engineering-20261027/` contains the completed 36-cell screen.
Historical reference auditing is separate. The corrected request reads annotator
campaign IDs from the accepted metrics snapshot, and synthetic conditions with
missing direct BCFs remain explicitly metadata-only. To retry only unresolved
reference lookups without repeating resolved hashes or annotation:

```sh
python -m benchmarks.performance_sprint.references audit NEW_AUDIT_DIRECTORY \
  --request CORRECTED_REQUEST.json \
  --only-unresolved-from ORIGINAL_AUDIT/references.jsonl
```

Supplemental summaries report their subset size and the full request size; they
do not pretend to independently re-audit all 540 historical conditions.

## Cached-only manuscript runner

`manuscript.py` is implemented, but **the full matrix is not submitted or approved**.
Tasks require explicit input/cache/annotation-asset fingerprints, frozen runtime
parameters, matched CPU count, historical cached duration and a baseline object.
Baseline approval binds its audit evidence files and direct duration. Its kind is
either `complete_records` (SHA-bound direct and historical cached BCFs) or
`historical_annotation_digest` (the original direct/cached digest, query format
and SHA-bound blocksort implementation). The latter never claims full-record
direct equivalence. New cached-reference differences always fail; inherited
annotator differences must remain explicitly declared.

The freeze document must bind the exact task-document SHA256, complete required
ID list, approved/tested candidate and immutable build file hashes. No implicit
direct measurements or automatic retries exist. After a reviewed deployment:

```sh
python -m benchmarks.performance_sprint.manuscript prepare NEW_CAMPAIGN \
  --tasks REVIEWED_TASKS.json --freeze APPROVED_FREEZE.json
python -m benchmarks.performance_sprint.manuscript submit CAMPAIGN
python -m benchmarks.performance_sprint.manuscript status CAMPAIGN
python -m benchmarks.performance_sprint.manuscript resume CAMPAIGN --cell EXPLICIT_FAILED_ID
```

Use worker-visible `/results` paths in documents, and mapped
`/mnt/data/slurm-results` paths for controller invocation. All conditions receive
individual exclusive jobs and an `srun` cgroup; historical cached duration sets a
30-minute to 24-hour limit including audit allowance. Longest conditions submit
first. Records, metrics, commands, comparisons and failure outputs persist;
successful large BCFs are not duplicated to shared storage. Slurm jobs/collectors
persist after disconnect, but candidate decisions and figure preparation are
not an unattended background agent.
