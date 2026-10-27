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
