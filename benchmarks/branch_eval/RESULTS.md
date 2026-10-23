# Branch evaluation: where each engine generation stands

2026-10-23. Four VCFcache engine generations against four annotator recipes at
five input scales, each compared only with its **own** uncached run. 119 timed
cells, one run per cell, 16 CPUs, exclusive nodes, complete-record audits
outside timing. Engineering evaluation; no manuscript figure or published
measurement was changed.

Arms: `main` (136bafa9) · `legacy@head` (legacy path on 63757a07) ·
`fused` (c878d86c) · `adaptive` (8830bf7b3 lineage) · `native` (63757a07).
`uncached` runs first in every job and is the semantic reference for the rest.

Annotator wall times are **never** compared across annotators: they do different
work, and vep-rs in particular behaves differently.

## Backend-label correction (integration review)

The tables below retain the original **requested arms**, not necessarily the
executed backend. All four native WGS requests, native VEP WES, and the native
500-µs heavy request **fell back to adaptive** because of the resident memory
budget. Their differences from the adaptive arm are repeated measurements of
the same backend, not native gains. Native actually ran for fastVEP and both
vep-rs WES recipes and for the small inputs.

See [the backend audit](integration-review-20261023/BACKEND_AUDIT.md) and its
machine-readable table: all 142 archived rows remain intact (119 primary
matrix cells, three heavy cells, 20 superseded cells). Existing figures and
raw tables have not been overwritten. The historical benchmark arm `auto`
means adaptive lookup; it is not the new availability-based CLI engine `auto`.

## Wall seconds

| recipe | scale | N | uncached | main | legacy@head | original fused | adaptive | native request |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| fastVEP clinical | 1 var | 1 | 10.89 | 15.41 | 15.72 | 61.25 | 1.88 | 1.31 |
| | panel | 281 | 13.86 | 25.51 | 25.97 | 25.87 | 12.91 | 11.19 |
| | WES | 82,714 | 48.92 | 52.82 | 52.48 | 47.29 | 44.97 | 38.55 |
| | WGS | 4,858,092 | 458.34 | 184.39 | 187.79 | 120.17 | 117.86 | 111.15† |
| VEP 115.2 everything | 1 var | 1 | 5.23 | 25.96 | 22.26 | 109.35 | 1.94 | 1.35 |
| | panel | 281 | 33.27 | 33.77 | 32.23 | 36.95 | 10.31 | 11.56 |
| | WES | 82,714 | 282.33 | 149.58 | 154.32 | 148.72 | 137.38 | 140.68† |
| | WGS | 4,858,092 | 6004.04 | 718.10 | 793.89 | n/r | 701.54 | 695.43† |
| vep-rs core | 1 var | 1 | 2.89 | 12.02 | 12.02 | 35.51 | 1.97 | 1.34 |
| | panel | 281 | 7.53 | 15.69 | 16.20 | 14.18 | 6.56 | 5.96 |
| | WES | 82,714 | 9.14 | 24.70 | 23.79 | 17.15 | 16.45 | 12.73 |
| | WGS | 4,858,092 | 68.16 | 103.96 | 105.29 | 49.86 | 44.98 | 45.79† |
| vep-rs clinical | 1 var | 1 | 1.93 | 14.12 | 13.07 | 44.57 | 1.89 | 1.22 |
| | panel | 281 | 7.91 | 18.69 | 17.32 | 16.82 | 6.56 | 6.20 |
| | WES | 82,714 | 69.64 | 29.14 | 28.28 | 25.79 | 24.52 | 20.36 |
| | WGS | 4,858,092 | 976.77 | 303.00 | 308.60 | 255.83 | 260.52 | 268.74† |

† Adaptive fallback, not resident native. A requested-native column is not a
single-backend performance series.

## Speedup over the same recipe's uncached run

| recipe | panel | WES | WGS |
| --- | --- | --- | --- |
| fastVEP clinical | 0.54 / 0.54 / 1.07 / **1.24** | 0.93 / 1.03 / 1.09 / **1.27** | 2.49 / 3.81 / 3.89 / 4.12† |
| VEP everything | 0.99 / 0.90 / **3.23** / 2.88 | 1.89 / 1.90 / 2.06 / 2.01† | 8.36 / n/r / 8.56 / 8.63† |
| vep-rs core | 0.48 / 0.53 / 1.15 / **1.26** | 0.37 / 0.53 / 0.56 / **0.72** | 0.66 / 1.37 / 1.52 / 1.49† |
| vep-rs clinical | 0.42 / 0.47 / 1.21 / **1.28** | 2.39 / 2.70 / 2.84 / **3.42** | 3.22 / 3.82 / 3.75 / 3.63† |

Order within each cell: main / original fused / adaptive / native request.
† Adaptive fallback. Bold highlights measured resident improvements, not
statistical significance. All clinical ratios remain provisional (see audits).

## Fixed cost and per-variant cost

The one-variant cell is a measured floor, not a fitted constant: a complete
invocation whose per-variant work is negligible. Slopes are differences between
adjacent measured scales.

| recipe | arm | floor @1 var | µs/variant panel→WES | µs/variant WES→WGS |
| --- | --- | ---: | ---: | ---: |
| fastVEP | uncached | 10.89 s | 425 | 86 |
| | native | **1.31 s** | 332 | 15 |
| VEP | uncached | 5.23 s | 3021 | 1198 |
| | native | **1.35 s** | 1566 | 116 |
| vep-rs core | uncached | 2.89 s | 20 | 12.4 |
| | native | **1.34 s** | 82 | 6.9 |
| vep-rs clinical | uncached | 1.93 s | 749 | 190 |
| | native | **1.22 s** | 172 | 52 |

VEP is markedly non-linear — 3021 µs/variant panel→WES against 1198 µs/variant
WES→WGS — so no single per-variant constant describes it. These are exploratory
end-to-end scale differences, not isolated annotation costs. In particular,
the native-request slopes cross resident/adaptive backends, and an all-hit
one-variant run skips annotation: its floor is not annotator startup.

## Heavy VEP pipelines: archived uncached arms reused

VEP `--everything` plus a `SyntheticDelay` plugin, KPGP-00319 WGS, 4,795,706
records, 90.26% hits, **8 CPUs / 8 vep_forks** to match the original campaign
(`wgs-pipeline-spectrum-af001-cd87186-kpgp00319`). The uncached arms cost 4.2 h
and 22.7 h and were reused rather than rerun; `main` was remeasured in the same
job as the anchor.

| rung | uncached (archived) | main archived | main fresh | adaptive | native |
| --- | ---: | ---: | ---: | ---: | ---: |
| +500 µs/consequence | 15,140.38 s | 1,805.37 s | 1,796.21 s | 1,687.21 s | **1,673.57 s** |
| speedup | — | 8.39× | 8.43× | 8.97× | **9.05×** |
| +10 ms/consequence | 81,740.53 s | 8,103.33 s | cancelled | cancelled | cancelled |

The fresh anchor agrees with the archive to **0.5%**, which is what makes the
splice defensible.

For the 500 µs adaptive arm, the measured miss-annotation stage is 1,537.50 s
of 1,687.21 s (91.1%); the native-request arm actually used adaptive and spent
1,526.19 s of 1,673.57 s (91.2%) there. This is stronger evidence than
subtracting misses multiplied by the uncached average cost: consequences per
variant and miss composition need not match that average. The ~1.018× bound
previously assigned to the cancelled 10-ms rung was **model-dependent**, not a
measured or rigorous ceiling. No new native-specific heavy gain was measured.

## Findings

1. **Adaptive/native substantially reduce fixed overhead and cache traversal.**
   The small all-hit floor falls to ~1–2 s, but the core-vep-rs WGS reduction
   (103.96 to 44.98 s) is larger than that floor alone can explain. This is not
   evidence for a purely fixed-cost change or a universal per-variant slope.
2. **`main` and `legacy@head` use the same legacy algorithm, but are not within
   1–3% everywhere.** VEP single differs by ~14%, clinical single/panel by
   ~7%, VEP panel by ~5%, and VEP WGS by 10.6% in these single runs. A narrow,
   reversed-order WGS confirmation is appropriate; do not merge their numbers.
3. **The archived original fused implementation must not become a default.**
   Its 109.35-s single-variant result shows the old full-cache scan defect.
   Current `--engine fused` is a compatibility alias for **adaptive**, which
   fixes that defect. Retiring the original algorithm does not mean deleting
   this alias or adaptive's implementation.
4. **Caching pays where per-variant cost is high.** VEP pays at every scale
   including panel; vep-rs core, at 12–20 µs/variant, pays only at WGS.
5. **Beating uncached vep-rs at WES depends entirely on the recipe.** Core loses
   (0.72×). Clinical — genome-wide gnomAD + ClinVar, 749 µs/variant panel→WES —
   wins at **3.42×**, and even main wins at 2.39×. The earlier negative result
   was specific to core and should not have been generalised.
6. **The measured heavy rung is dominated by miss annotation (~91%).**
   Further routing work has limited leverage there. Cohort batching merits a
   separate, fairly controlled experiment, not another engine rewrite.

## Output discrepancies: annotator-side, not routing

| recipe | WES | WGS | cause |
| --- | --- | --- | --- |
| VEP everything | 2/82,714 `HGNC_ID` | 28/4,858,092 | [ensembl-vep#1959](https://github.com/Ensembl/ensembl-vep/issues/1959): `merge_features` back-fills HGNC IDs by gene symbol from other features in the same buffer, so any subset annotation can differ from a full-input run |
| vep-rs clinical | 22/82,714 ClinVar CSQ | 40/4,858,092 | dbNSFP input-context omission reproduced without VCFcache; precise internal mechanism unproven |

Every cached arm reported **byte-identical** mismatch counts and loci, so
routing is not implicated in either case. Both are rarer genome-wide than at
WES — the discrepancies concentrate in coding regions. The vep-rs clinical
timings are therefore sound measurements but **not** quotable as validated
speedups until that cause is understood.

Focused probes recovered the same 22 CSQs from the original cache payload.
Fresh annotation of those exact sample or blueprint records matches direct
annotation, including a repeat, one thread, and removal of GT/INFO/ID. All
22 still match direct with 1-kb blueprint neighborhoods but match the stored
cache with 100-kb neighborhoods (3,110 input records). The tool and ClinVar
data/index hashes agree with build provenance. This establishes recipe
input-context sensitivity, not a routing defect; it does not yet establish
the source-code mechanism. No comparison fields were ignored or weakened.

The [integration review](integration-review-20261023/REVIEW.md) includes a
three-variant reproducer: both preceding records suppress the target's three
ClinVar fields, while target alone or either preceding record does not. This
holds at one and 16 threads, including dbNSFP without gnomADc. Clinical speedups
remain provisional until a context-invariant recipe passes full-record audits.

## Limitations

- One run per cell. The two VEP WES direct measurements differ by 3%, but
  this is not a campaign-wide variance estimate. Differences of a few percent between engine arms
  are not an ordering. At vep-rs clinical WGS the three new engines span 5%
  (255.8–268.7 s) and should be read as tied.
- No cross-annotator comparison is valid.
- The heavy rungs use a different sample and 8 CPUs, not 16.
- `main` has no `--engine` option, so its arm omits that flag.

## Contents

- `RUN_MANIFEST.json` — every job including the four that failed and why
- `cells.tsv` — 142 timed cells, machine-readable, superseded jobs retained
- `SUMMARY.txt` — generated tables
- `inputs/INVENTORY.tsv` — record and contig counts per scale
- `<condition>-<job>/` — per-cell manifests, resources, comparisons, logs

Reproduce: `benchmarks/branch_eval/{prepare_inputs,eval,heavy}.sh`, summarised
with `benchmarks/branch_eval/summarize.py`.
