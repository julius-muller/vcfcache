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

## Wall seconds

| recipe | scale | N | uncached | main | legacy@head | fused | adaptive | native |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| fastVEP clinical | 1 var | 1 | 10.89 | 15.41 | 15.72 | 61.25 | 1.88 | 1.31 |
| | panel | 281 | 13.86 | 25.51 | 25.97 | 25.87 | 12.91 | 11.19 |
| | WES | 82,714 | 48.92 | 52.82 | 52.48 | 47.29 | 44.97 | 38.55 |
| | WGS | 4,858,092 | 458.34 | 184.39 | 187.79 | 120.17 | 117.86 | 111.15 |
| VEP 115.2 everything | 1 var | 1 | 5.23 | 25.96 | 22.26 | 109.35 | 1.94 | 1.35 |
| | panel | 281 | 33.27 | 33.77 | 32.23 | 36.95 | 10.31 | 11.56 |
| | WES | 82,714 | 282.33 | 149.58 | 154.32 | 148.72 | 137.38 | 140.68 |
| | WGS | 4,858,092 | 6004.04 | 718.10 | 793.89 | n/r | 701.54 | 695.43 |
| vep-rs core | 1 var | 1 | 2.89 | 12.02 | 12.02 | 35.51 | 1.97 | 1.34 |
| | panel | 281 | 7.53 | 15.69 | 16.20 | 14.18 | 6.56 | 5.96 |
| | WES | 82,714 | 9.14 | 24.70 | 23.79 | 17.15 | 16.45 | 12.73 |
| | WGS | 4,858,092 | 68.16 | 103.96 | 105.29 | 49.86 | 44.98 | 45.79 |
| vep-rs clinical | 1 var | 1 | 1.93 | 14.12 | 13.07 | 44.57 | 1.89 | 1.22 |
| | panel | 281 | 7.91 | 18.69 | 17.32 | 16.82 | 6.56 | 6.20 |
| | WES | 82,714 | 69.64 | 29.14 | 28.28 | 25.79 | 24.52 | 20.36 |
| | WGS | 4,858,092 | 976.77 | 303.00 | 308.60 | 255.83 | 260.52 | 268.74 |

## Speedup over the same recipe's uncached run

| recipe | panel | WES | WGS |
| --- | --- | --- | --- |
| fastVEP clinical | 0.54 / 0.54 / 1.07 / **1.24** | 0.93 / 1.03 / 1.09 / **1.27** | 2.49 / 3.81 / 3.89 / **4.12** |
| VEP everything | 0.99 / 0.90 / **3.23** / 2.88 | 1.89 / 1.90 / **2.06** / 2.01 | 8.36 / n/r / 8.56 / **8.63** |
| vep-rs core | 0.48 / 0.53 / 1.15 / **1.26** | 0.37 / 0.53 / 0.56 / **0.72** | 0.66 / 1.37 / **1.52** / 1.49 |
| vep-rs clinical | 0.42 / 0.47 / 1.21 / **1.28** | 2.39 / 2.70 / 2.84 / **3.42** | 3.22 / **3.82** / 3.75 / 3.63 |

Order within each cell: main / fused / adaptive / native.

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
WES→WGS — so no single per-variant constant describes it.

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

On the 500 µs rung the cache-side budget is 1,796.21 − (467,102 misses ×
3.157 ms) ≈ **321 s, 17.9%** of the cached run; native recovered **123 s of it,
38%**, against a 10.27× ceiling. The same arithmetic puts the 10 ms rung's
cache-side at ~144 s, **1.8%**, capping any engine gain at ~1.018×; that rung was
cancelled after 2 h 19 min rather than spend ~4.5 h more confirming a bound.

## Findings

1. **The engines removed a fixed cost, not a per-variant cost.** Every generation
   collapses the one-variant floor to ~1.3 s where main needs 12–26 s and fused
   needs 36–109 s. That single change is what turns VEP panel from 0.99× to
   3.23× and vep-rs core WGS from 0.66× to 1.52×.
2. **`main` ≈ `legacy@head`.** Within 1–3% at every cell except VEP WGS
   (718.10 vs 793.89 s, 10.6%), which is one run and needs a repeat. Archived
   "legacy" numbers can otherwise be read as main.
3. **`fused` must not become a default.** 109.35 s for one VEP variant against a
   5.23 s uncached run — the `workers == 1` full-cache scan. It is fine at WGS
   (3.81× on fastVEP), so the defect is confined to small inputs.
4. **Caching pays where per-variant cost is high.** VEP pays at every scale
   including panel; vep-rs core, at 12–20 µs/variant, pays only at WGS.
5. **Beating uncached vep-rs at WES depends entirely on the recipe.** Core loses
   (0.72×). Clinical — genome-wide gnomAD + ClinVar, 749 µs/variant panel→WES —
   wins at **3.42×**, and even main wins at 2.39×. The earlier negative result
   was specific to core and should not have been generalised.
6. **At the heavy end the engines are bounded by arithmetic, not engineering.**
   1.07× at the 500 µs rung, ≤1.02× at 10 ms, because cached runtime there is
   82–98% miss annotation.

## Output discrepancies: annotator-side, not routing

| recipe | WES | WGS | cause |
| --- | --- | --- | --- |
| VEP everything | 2/82,714 `HGNC_ID` | 28/4,858,092 | [ensembl-vep#1959](https://github.com/Ensembl/ensembl-vep/issues/1959): `merge_features` back-fills HGNC IDs by gene symbol from other features in the same buffer, so any subset annotation can differ from a full-input run |
| vep-rs clinical | 22/82,714 ClinVar CSQ | 40/4,858,092 | not established |

Every cached arm reported **byte-identical** mismatch counts and loci, so
routing is not implicated in either case. Both are rarer genome-wide than at
WES — the discrepancies concentrate in coding regions. The vep-rs clinical
timings are therefore sound measurements but **not** quotable as validated
speedups until that cause is understood.

## Limitations

- One run per cell; observed spread 1–3% (VEP WES uncached 273.98 s in job 1301
  against 282.33 s in 1315). Differences of a few percent between engine arms
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
