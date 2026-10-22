# VCFcache fixed overhead, profiled on real panel/WES scale

Slurm jobs **1253** (size sweep, real assays) and **1254** (span sweep),
2026-10-22, `sl-w1`/`sl-w2`, exclusive, 16 vCPU. One untimed warm-up then three
timed reps per cell; medians below.

Code: `feature/fused-routing` at `c878d86c3`, deployed to the cluster and
shadowed over the stale editable install via `PYTHONPATH`. `vcfcache --version`
reports 0.5.3 there — that is the old dist-info, not the code; both jobs assert
the branch by requiring `--engine` in the CLI.

Arm A and the span sweep draw their inputs **from the cache**, so every variant
is a hit, the annotation recipe is never invoked, and what is timed is
VCFcache's own cost and nothing else.

## Fixed overhead is essentially the whole cost at panel scale

Zero misses, inputs spanning all 22 autosomes (`N=1` is single-contig):

| N | legacy | fused |
|---|---|---|
| 1 | 24.0 s | 107.6 s |
| 250 | 23.7 s | 45.6 s |
| 2,000 | 22.0 s | 45.9 s |
| 75,000 | 27.0 s | 47.6 s |

Fitting the two ends: **legacy = 23.6 s + 45.1 µs/variant**, fused = 45.6 s +
26.8 µs/variant. The per-variant figures sit in the same range as the production
27.5 µs and reproduce the fused engine's ~2× per-variant win (45.1 → 26.8 µs).

But the split is the point:

- At **panel scale (250 variants)** the fixed term is ~**100%** of VCFcache's cost.
- At **WES scale (75,000)** it is still **87.5%** (legacy).

Independent confirmation of the model: legacy at `N=1` measured **24.0 s**
against the 23.6 s fitted intercept and the 21.7 s that
`benchmarks/analyze_breakeven.py` fit from production data. The fixed term is
real and correctly sized — now observed directly rather than extrapolated to.

## The decisive test: it is a cache scan, and region restriction already works

250 variants held constant, only the contig span varied (job 1254):

| contigs | legacy | fused |
|---|---|---|
| 1 | 19.6 s | **111.1 s** |
| 2 | 19.5 s | **10.2 s** |
| 4 | 19.7 s | **10.4 s** |
| 8 | 20.9 s | 18.6 s |
| 22 | 21.3 s | 45.7 s |

**Legacy is flat at ~20 s for 1 to 22 contigs.** It streams the whole 18.87 M-record
cache whatever the input touches; span cannot move it because nothing is
restricted. That ~20 s is a full-cache pass, and it is paid to annotate 250
variants.

**Fused's route stage is region-restricted and contig-parallel, and at 2–4
contigs it halves legacy: 10.2 s against 19.5 s.** That is the proof that the
lever is real — restricting cache access to the input's regions cuts the fixed
term by 2× already, with no further work.

It is wrecked at both ends:

- **Single-contig is pathological: 111.1 s, 5.7× worse than legacy.** The stage
  breakdown puts 109.9 s of it in `route` alone. The engine's own README notes a
  single-contig input cannot exploit contig parallelism, which explains the lost
  speedup but not a 5.7× regression against legacy. This needs reading the route
  path, not more measurement.
- **At 22 contigs the four workers are badly unbalanced** — per-worker route
  times 9.9 / 21.3 / 32.6 / 45.6 s, and wall time is set by the slowest. Total
  route work is ~109 s over 4 workers; balanced across 8 it would be ~14 s.
  Adjacent-contig grouping by record count is not achieving balance here.

So fused's 45.6 s fixed term at genome-wide span is **worse than legacy's 23.6 s**,
despite route being the right idea, purely because of scheduling.

## Real assays: where the cached second actually goes

Actual ACMG SF v3.3 panel and Twist core exome samples, which have misses and so
do invoke VEP. `VCFcache` column is the matched zero-miss cell; `VEP` is the
remainder:

| assay | variants | misses | legacy total | fused total | VCFcache | VEP |
|---|---|---|---|---|---|---|
| panel | 203 | 3 | 27.0 s | 41.7 s | 23.7 s | 3.3 s |
| WES | 75,135 | 1,339 | 156.3 s | 177.0 s | 27.0 s | 129.3 s |

**At panel scale, 88% of the cached run is VCFcache overhead and 12% is actual
annotation.** The campaign's panel loss (0.41–0.46×, cached 22.4 s against direct
10.5 s) is therefore almost entirely this full-cache pass — it is not annotator
start-up, and it is not per-variant work. At WES the balance inverts: VEP is 83%
and overhead 17%, which is why WES sits near parity.

## What this means for the earlier recommendation

The proposal to attack the fixed term rather than the per-variant slope is
supported, and sharpened:

1. **Region-restricted routing is worth ~2×, measured** (10.2 s vs 19.5 s at 2–4
   contigs). It is already implemented in the fused engine's route stage.
2. **Two defects stop it paying off** at the spans real assays actually have: the
   single-contig regression and worker imbalance at high span. Both are
   scheduling, not algorithmic.
3. **Fused routing must not be promoted to default on this evidence.** For a
   genome-wide panel it is 1.5–2× *slower* than legacy (41.7 s vs 27.0 s), and
   for single-contig input 5.7× slower. It wins on per-variant cost, which is
   the term that does not matter at these scales.

A fixed term near 10 s would move panel from 0.46× to roughly break-even against
a 10.5 s direct run; eliminating the full-cache pass entirely is what would make
panel caching clearly profitable.

## Caveats

- **Hit rates here are not cohort results.** These are 1000 Genomes samples
  against a gnomAD AF>=1% cache, giving 98.2–98.5% hits — they are not held out
  from gnomAD. Overhead is hit-rate independent, so the profiling stands, but
  these percentages must not be quoted as coverage findings. The held-out
  numbers remain KPGP/SGDP/PGP at 90–96%.
- **Arm C of job 1253 is excluded.** Its three inputs were all single-contig
  whatever span they claimed, because the genome-wide pool is contig-ordered and
  the selection took the first 250 matching rows. Caught from the inventory the
  job prints before timing. Job 1254 supersedes it and asserts realised span
  before measuring. The dead cells are left in the raw TSV rather than deleted.
- Two-point fits over N=250..75000; the `N=1` cell is reported separately because
  it is single-contig and so not on the same span curve.
- Job 1252 died silently on the `pipefail` + early-`exit` trap (`awk {exit}`
  SIGPIPEs `bcftools`); proven in isolation at rc=141 against rc=0. Fixed by
  draining, and the cache scan replaced with 22 indexed windows.
