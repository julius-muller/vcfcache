# Strong-candidate follow-up: prepared Python versus Rust/native

This second bounded pilot addresses the omissions in the first default-adaptive
confirmation. It **does not overwrite that campaign**, change vcfcache defaults,
or start the full publication cohort. The runner's `--profile candidates`
freezes this design before submission, reusing the existing source deployments,
strict comparator, independent timed Slurm steps and durable review mechanism.

## Fixed conditions and candidates

All measurements use eight allowed CPUs, bcftools 1.25 and, for the Rust helper,
linked HTSlib 1.25. Both receive the **same sorted/indexed BCF**, verified to be
biallelic, star-free, duplicate-free and without the required annotation tag.
No timed arm normalizes/splits/removes records. The common input-domain check
is untimed and its duration is saved in `prepared_preflight.json`.
Thus these are **prepared-BCF annotation timings**, not a claim that arbitrary
raw VCF requires no preparation or validation. The eligibility cost must be
considered before promoting a fast path into a general product default.

| Assay | Prespecified cached candidates |
|---|---|
| WES | Python prepared windows, 2 workers; Rust resident/native 4 workers + 0 decoders, 4 + 4, and 8 + 0; resident 4 + 4 with CSI chunk coalescing |
| WGS | Python prepared auto routing; Python prepared whole-contig routing with 2 workers; Rust prepared adaptive routing |

Each block also includes **fresh direct annotation**. Native WGS is deliberately
absent because its known adaptive fallback would not measure the native engine.
Native memory budget is fixed at 2048 MiB and fallback is a failure, not an
accepted alternate implementation. Requested worker/decoder settings and
prepared-input activation are verified from actual run statistics.

The candidate list is based on earlier measurements: windows2 beat the other
Python WES settings; resident WES substantially beat Rust-adaptive; worker and
decoder allocation mattered; coalescing helped some core WES cells. WGS Python
auto/contigs2 are included because their ordering can depend on CPU scheduling.
The unsuccessful Python publication-pass fusion experiment is not repeated.
This is **not an exhaustive search or proof of a universal fastest configuration**:
older library versions, every worker count, other recipes and resident memory
stores are outside this small pilot. Every measured candidate remains visible;
we will not publish a per-sample minimum assembled across configurations.

## Six blocks, run concurrently

| Worker | Recipe | Sample | Assay |
|---|---|---|---|
| sl-w1 | vep-rs core | KPGP-00139 | WES |
| sl-w2 | fastVEP clinical | KPGP-00139 | WES |
| sl-w3 | vep-rs core | LP6005592-DNA_C05 | WES |
| sl-w4 | fastVEP clinical | LP6005592-DNA_C05 | WES |
| sl-w5 | vep-rs core | KPGP-00139 | WGS |
| sl-w6 | fastVEP clinical | KPGP-00139 | WGS |

The two real WES inputs include the typical and lower-hit cases from the earlier
pilot; actual hit counts are recorded anew. Three repetitions per condition give
**96 timed invocations**, with **90 strict comparisons** including direct-repeat
stability. Most extra conditions are short WES runs; no cohort expansion occurs.
Order rotates across repeats and blocks, but three repetitions cannot perfectly
balance every position for four or six modes. Report ranges and paired results;
do not equate technical repeats with biological replication or confidence limits.

Identical cache payloads and working files use private node-local scratch.
Annotator database paths remain as frozen, including fastVEP's NFS supplementary
annotation directory. This does not remove every shared-host/storage effect.
No page caches are flushed. Asset hashes, staging, input extraction and shared
eligibility verification are outside timing; complete output hashing and strict
audits follow all timed runs. Statistics are light. INFO/FORMAT and relevant
headers are compared with **zero ignored fields**; a mismatch blocks acceptance.

The first campaign's historical direct times are **not reused** or pooled here.
The follow-up has its own immutable deployment/result directory and submission
receipts. A three-hour per-block cap is only a safety limit; expected completion
is approximately 1–2 hours including audits. Early-finished WES workers may be
idle because this is a bounded pilot, not a request to generate additional work.

## Execution and later review

Controller paths start `/mnt/data/slurm-results/engine-candidates-20261109-...`;
worker paths start `/results/engine-candidates-20261109-...`. Launch IDs and
manifest/source/tool hashes are saved in a separate candidate launch receipt.
`state.json` checkpoints show phase, completed timings and audits. Large outputs
stay in owned local scratch for review. Do not purge them before acceptance.

An after-any review job automatically creates `RESULTS.md`, `summary.tsv`,
`conditions.tsv`, `REVIEW.json` and a hash-bound evidence archive. It does not
change defaults, delete branches or submit a publication campaign. Those remain
separate user decisions after inspecting both correctness and paired performance.

```bash
# Explicit controller-only freeze before submission (as appuser):
PYTHONPATH=DEPLOYMENT PYTHON -m benchmarks.engine_confirmation.run freeze RESULT_ROOT \
  --deployment WORKER_DEPLOYMENT --commit COMMIT --profile candidates
# Review can be repeated days later without Slurm accounting:
PYTHONPATH=DEPLOYMENT PYTHON -m benchmarks.engine_confirmation.run collect RESULT_ROOT
```
