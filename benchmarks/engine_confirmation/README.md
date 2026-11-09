# Bounded same-worker engine confirmation

The subsequent prepared-input/window/resident comparison is documented
separately in [CANDIDATES.md](CANDIDATES.md); it does not replace this first run.

This is the small decision pilot authorized on 9 November 2026, **not a rerun
of the manuscript cohort**. It tests whether the apparent Python/Rust gap
survives fresh controls on the upgraded cluster. It does not change defaults,
retire branches, optimize engines or automatically start another campaign.

## Preregistered design

Six exclusive worker blocks, using every worker once:

| Worker | Recipe | Sample | Assay |
|---|---|---|---|
| sl-w1 | vep-rs core | KPGP-00139 | WGS |
| sl-w2 | fastVEP clinical | KPGP-00139 | WGS |
| sl-w3 | vep-rs core | LP6005519-DNA_A06 | WGS |
| sl-w4 | fastVEP clinical | LP6005519-DNA_A06 | WGS |
| sl-w5 | vep-rs core | KPGP-00139 | WES |
| sl-w6 | fastVEP clinical | KPGP-00139 | WES |

Each block runs three repeats of **fresh direct**, **Python-adaptive**, and
**Rust-adaptive**: 54 complete CLI invocations altogether. The three-condition
order rotates so each mode occupies each position once. Repeats are only for
engineering timing variability; two WGS genomes/one WES genome do not establish
a population-level confidence interval. Different recipes on different nodes
are not a controlled cross-annotator ranking.

The cached engines are the frozen candidates from the preceding toolchain
pilots: Python source `b9633e72dba1a40dd6909200a58f7430a16a1fcf` and Rust
source `533dcd360` (full source-file hashes bound in the manifest). Both use
the same bcftools/HTSlib **1.25**, and Rust links HTSlib 1.25 with libdeflate.
No best-per-assay tuning, publication-pass fusion, prepared-input bypass,
experimental memory store or resident/native WES is enabled. Rust means the
adaptive engine, not an implicit native fallback. Branch-default preprocessing
is retained in all timed arms; this is not a normalization-removal experiment.
Assay extraction/BCF encoding before timing does not alter input variants.

CPU affinity is enforced with a fresh `srun --exact -c8 --cpu-bind=cores`
step per timed CLI, not merely configured thread counts. Memory is 28,000 MiB;
the three-hour limit is a safety cap, not an expected runtime. All annotator
parameters and cache/input/asset payloads are frozen and hash-checked. Cache
payloads are copied **byte-identically to private local scratch** before timing,
and outputs/intermediates also stay local. This removes concurrent NFS cache
lookup as an engine confound, but changes storage relative to historical shared
cache campaigns. It is explicitly a confirmation under local-cache conditions.
Annotator databases retain their frozen recipe locations: for example,
fastVEP's supplementary annotation directory is still on `/results` (NFS).
This does not eliminate every source of shared-storage contention; it isolates
the cache payload location and records the remaining recipe paths. Repeats and
counterbalanced order reduce, but cannot rule out, time-varying contention.
No OS page-cache flushing or cold-cache claim is made. Full output hashes and
semantic audits run only after all timings; timed statistics are light.

The replacement sl-w5 is eligible only after verifying Broadwell, the matching
kernel/container versions, helper binary and Python runtime. Its historical
AMD timings remain excluded from the earlier v2 analysis. Preflight found the
old benchmark Python 3.12.3 interpreter missing on the replacement worker; the
existing interpreter distribution was copied there without modifying system
packages or any frozen benchmark environment.

## Acceptance and review days later

Each cached output is compared against its **own repetition's direct output**.
The later direct outputs are compared to the first direct, checking repeat
stability as well. Comparison includes complete records, INFO/FORMAT and
relevant headers, canonicalizing only the existing comparator's ordering
rules. No HGNC_ID/ClinVar or other fields are ignored. Every output must be
indexed and retain the expected record count; Rust may not silently fall back.
A failed audit preserves timings but does not produce accepted speedup claims.

Large outputs stay in owned worker scratch; do not purge them before review.
Shared evidence includes manifest/READY hashes, submitted job IDs and commands,
hardware/mount/boot identity, resource/cgroup samples, CLI logs, output hashes,
and atomic state checkpoints after every condition/audit. An **after-any review
job** writes `REVIEW.json`, `conditions.tsv`, `summary.tsv`, `RESULTS.md`, and a
hash-bound evidence archive even if a job fails or is interrupted. It never
submits a further benchmark. Review works without historical Slurm accounting.

Worker-visible results use `/results/engine-confirmation-20261109-*`; the
controller exposes the same tree at `/mnt/data/slurm-results/`. Exact launch
receipts are saved alongside this README after submission.

Accepted launch: **jobs 2842–2847**, one block on each worker in the table,
with automatic after-any review **2848**. The result root is
`/results/engine-confirmation-20261109-62a9fdb98` (controller:
`/mnt/data/slurm-results/engine-confirmation-20261109-62a9fdb98`). Frozen runner
source is `62a9fdb98`; later documentation commits do not alter timed sources.
[launch-20261109.json](launch-20261109.json) binds the manifest SHA-256,
submission commands/IDs, tool binaries, backend commits and successful imports
on all six nodes. Initial launch observation showed all six jobs running in
asset preflight; that alone does not establish successful timings or audits.
The full node is reserved exclusively (Slurm may display 16 allocated CPUs),
while every timed condition must report exactly eight allowed CPUs.

Initial health check subsequently confirmed **vep-rs core WES completed all
nine timings and all eight strict comparisons passed**. The other five blocks
were still running, with no recorded failures. See
[initial-health-20261109.json](initial-health-20261109.json). This is a launch
health receipt, not the final pilot conclusion; use the shared `REVIEW.json`
and `RESULTS.md` after review job 2848 finishes.

```bash
# On the controller, as appuser, using the frozen deployed runner:
PYTHONPATH=DEPLOYMENT PYTHON -m benchmarks.engine_confirmation.run collect RESULT_ROOT
```

Decision rule: inspect strict acceptance first, then per-block median gains,
individual repeats/ranges and consistency across both WGS genomes. A roughly
5% gain alone does not justify another backend; reproducible gains around
20% or a reversal of a caching loss are decision-relevant. These are engineering
criteria, not statistical significance thresholds. Packaging is checked on the
deployed Linux workers only; this does not validate other operating systems.
The subsequent implementation/publication decision remains with the user.
