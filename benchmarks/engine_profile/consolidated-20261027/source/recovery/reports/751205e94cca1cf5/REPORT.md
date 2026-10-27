# Engine overhead profile

Engineering diagnostics, not new publication speedup claims.

Cells: 8. States: {"complete": 5, "incomplete_timeout": 3}.

Read `cells.tsv` for requested/effective backends, strict audits and resources; `phases.tsv` for inclusive stage timings.
Intervals are unioned within each phase. Different phases can overlap; do not stack or sum them as disjoint wall time.
Unattributed time excludes container/wait envelopes and includes process startup and uninstrumented boundaries. Compression and decoding CPU attribution requires a usable stack profile; coarse overlapping wall spans cannot establish it.

All attempts, commands, partial traces and raw audit failures are retained. Missing or running cells are not successes.
Cold storage is NOT claimed: staging and hashing warm files; first-use-after-copy is an annotation-access condition only.

## Integrity gates

Repeated cgroup paths (condition peak-memory values invalid if nonempty): []
Trace parse failures: 0.
Direct self-checks are reference sanity checks, not independent equivalence proof. `complete_with_direct_difference` remains a failed direct comparison; no fields are ignored.

## Optimization screening

Completed instrumentation pairs: 0/3. See `instrumentation.tsv`; absolute paired shifts above 5% require timing-noise/instrumentation review before interpreting stage fractions.
Rank measured cache-side phases only after checking instrumentation-on/off pairs and CPU/worker controls. A removable phase is an optimistic bound, not a promised speedup.
| Cell | Phase | Union wall s | Fraction of end-to-end |
| --- | --- | ---: | ---: |

Candidates require output equivalence, reproducible end-to-end benefit and a separate implementation review. Native fallback waste and memory pathology remain relevant even below a 5% phase threshold. No optimization is implemented by this campaign.

## Review days later

Run `python -m benchmarks.engine_profile.run collect CAMPAIGN` again after transferring the campaign directory if desired. No worker scratch or Slurm accounting is required for recorded cells.
Plot with `Rscript benchmarks/engine_profile/plot.R REPORT_DIRECTORY`; only R/ggplot2 is used.
