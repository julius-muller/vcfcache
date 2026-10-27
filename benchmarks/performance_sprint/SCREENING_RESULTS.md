# Screening interpretation — 27 October 2026

These are engineering controls on two pinned genomes/nodes, not additional
biological replicates or publication estimates. All durations below are timed
cached CLI steps, in seconds. Direct annotation was not rerun. Acceptance needs
reproducible end-to-end benefit, not a favorable single stage.

| Candidate / workload | A06 control → candidate | C05 control → candidate | Interpretation |
|---|---:|---:|---|
| Indexed miss spools / fastVEP WES | 47.05 → 48.95 | 55.92 → 55.47 | No consistent benefit |
| Indexed miss spools / core vep-rs WGS | 45.75 → 43.59 | 60.62 → 60.89 | Below gate; inconsistent |
| 4-kb windows / native fastVEP WES | 44.08 → 46.30 | 43.02 → 49.61 | Regresses; retain 16-kb default |
| 1-kb windows / native fastVEP WES | 44.08 → 42.34 | 43.02 → 55.02 | No consistent benefit |
| Trusted hit assembly / fastVEP WES | 47.75 → 55.10 | 55.13 → 62.08 | No demonstrated end-to-end benefit |
| Trusted hit assembly / core vep-rs WGS | 44.08 → 42.69 | 61.11 → 61.78 | Below gate; inconsistent |
| Early memory fallback / core vep-rs WGS, 8-GiB native budget | 62.68 → 45.90 | 79.32 → 62.22 | 26.8% / 21.6% reduction; confirmation submitted |
| Bounded native route / fastVEP WES | 51.91 → 47.81 | 52.51 → 49.67 | 7.9% / 5.4% reduction, sparse workload only |
| Bounded native route / core vep-rs WGS | 43.90 → 54.11 | 60.44 → 79.17 | 23.3% / 31.0% slower; not a universal default |
| Early memory fallback / reversed-order A06 confirmation | 65.33 → 45.60 | Not repeated | 30.2% reduction; strict parity passes |

All **36/36** engineering cells completed and passed strict complete-record/header checks
against both retained direct output and the frozen cached reference. No ignored
fields were needed in these core/fastVEP pilot cells. Clinical vep-rs inherited
discrepancies are not being re-investigated or excused for new engine changes.

The trusted-hit fastVEP total-time regressions do not establish slower assembly:
most of the difference was annotator-invocation drift. This is precisely why
single-stage micro-optimizations do not meet the end-to-end promotion gate.

Early fallback avoids allocating/filling a resident buffer that demonstrably
cannot fit, then repeating lookup adaptively. It does **not** reduce annotation
cost per miss or imply a 21–27% gain when adaptive routing is already selected.
The predictor samples <=128 records, abstains on weak evidence and reports a
heuristic byte estimate rather than an RSS guarantee. It remains experimental.

All engineering jobs and collectors finished. Job **1368** reversed the A06
memory-screen order (candidate first, control second) and confirmed the largest
supported effect. No further speculative engine cells are scheduled. This
supports selecting the fallback guard for the release-candidate review; it does
not establish a gain when adaptive routing is already selected.

Final engineering evidence: `evidence-final-engineering-20261027/`. The older
`evidence-20261027T1122/` snapshot intentionally retains its then-pending states.
Launch receipts bind the immutable remote campaigns. Historical-output identity
audit **1367** finished separately; corrected reference resolution is being
audited read-only as **job 1370**, without repeating annotation.

## Decision before manuscript replay

Do not promote indexed spooling, smaller windows or trusted assembly based on
these results. Do not make the bounded router universal: a sparse-input gain
does not offset a dense-input regression. Review the final fallback confirmation
and its applicability when selecting the release candidate. The previous
adaptive/native engine results remain distinct from these incremental screens.

The cached-only manuscript runner is implemented, not yet submitted. Complete
per-row baseline approvals, preserve digest-only validation limitations, freeze
the candidate/helper/recipes, materialize the complete historical matrix plus
the pinned vep-rs extension, then launch the manuscript campaign. The Monday
deadline remains a target, not a claim that these future stages already ran.
