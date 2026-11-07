# Rust / HTSlib 1.25 follow-up

Isolated pilot based on `feature/rust-wes-pilot-20261103` at `3e9f6a6ea`.
Main, previous source deployments, cache formats, immutable recipes, and prior
results are unchanged. No public default changes or full-cohort jobs are made.

## Questions and controls

Both helpers use identical Rust sources and Cargo dependency versions, including
the same libdeflate and zlib-ng dependencies. One uses bundled HTSlib 1.19.1;
the other uses upstream 1.25 from the SHA-bound bcftools release archive.
An isolated hts-sys copy supplies 1.25 sources; bindgen regenerates ABI bindings.
The Cargo registry is never edited. `build.py` records binary hashes, linked
versions/features, source hashes, dependency lock and commands. There is no
blind runtime library swap or experimental Cargo feature enabled.

Each WES block measures, once, these nine conditions:

| Condition | External bcftools | Linked helper HTSlib | Input path |
|---|---|---|---|
| adaptive-control121 | 1.21 | 1.19.1 | Original |
| adaptive-external125 | 1.25 | 1.19.1 | Original |
| adaptive-linked125 | 1.25 | 1.25 | Original |
| adaptive-prepared125 | 1.25 | 1.25 | No rewriting |
| resident-control121 | 1.21 | 1.19.1 | Original FIFO |
| resident-external125 | 1.25 | 1.19.1 | Original FIFO |
| resident-linked125 | 1.25 | 1.25 | Original FIFO |
| resident-prepared125 | 1.25 | 1.25 | Original BCF directly |
| direct-prepared125 | 1.25 | Not used | No rewriting |

These pairs isolate the external upgrade, linked-library upgrade, and finally
input preservation. No codec-pool tuning is added. Resident mode retains four
lookup workers and its original CPU-budgeted decoder pool. Rust-adaptive retains
its existing region-density planner and dynamic scheduling: no assay-specific
algorithm or selection of the fastest setting per sample.

`VCFCACHE_PILOT_PREPARED=1` is an experimental trusted-domain assertion. The
runner checks biallelic, non-star, unique exact variant keys and absence of the
required annotation **outside timing**. Candidate direct and cached runs use
the same immutable indexed BCF. Adaptive preparation still writes the indexed
spool needed by its current planner, but does not split/filter input and uses
exact cache/miss overlays. Resident preparation omits its normalizer process
and FIFO entirely. Unsupported inputs fail preflight rather than being silently
transformed. This is not generic production support for arbitrary multiallelic,
duplicate, or pre-annotated input. Original control preprocessing is retained;
all fold ratios are explicitly against **prepared direct**, not a normalized
direct arm that was never timed.

## Bounded design and durable review

- GRCh38, two existing public inputs: KPGP-00139 and lower-hit
  LP6005592-DNA_C05 WES; KPGP-00139 WGS. fastVEP clinical and vep-rs core.
- Four WES blocks × nine conditions, two WGS blocks × five conditions:
  **46 CLI timings, 40 cached-versus-direct comparisons**.
- Resident WGS is deliberately excluded: its memory fallback would measure
  adaptive rather than native and must not be labelled as a resident gain.
- Five Intel Broadwell workers (sl-w1/2/3/4/6), never AMD sl-w5. The sixth block
  follows the first on sl-w1; other idle workers after their assigned work finish
  do not imply an unscheduled full cohort. Every timing gets a fresh eight-CPU
  Slurm step and cgroup counters; nodes are exclusive. Order rotates by block.
- Recipes/caches/assets are frozen from the existing manuscript manifest.
  Source/binary/input/cache hashes, hardware/storage, commands, counts,
  resources, effective backend and comparisons are retained. Unexpected native
  fallback and strict comparison failure reject the block.
- Full record/relevant-header comparison uses no ignored fields and runs only
  after all timings in each block. Preparation, hashes, output audits and count
  verification are untimed; normal CLI operations remain timed. File caches are
  not flushed; shared NFS noise and single measurements limit inference.
- After-any collection writes `RESULTS.md`, `conditions.tsv`, `REVIEW.json`,
  `SNAPSHOT.json` and a terminal evidence archive. A bounded five-hour controller
  watcher renders R/ggplot2 PDF/SVG/600-dpi PNG panels and a combined figure into
  a new folder. No postprocessor submits annotation or full-cohort work.
- Two-hour job safety limit; expected completion including audits is roughly
  30–60 minutes after submission, not a guarantee. Scratch survives failures.
  The after-any review explicitly labels missing/nonterminal runner checkpoints
  as interrupted; it never silently treats a timed-out job as completed.

Build in a new owned temporary directory with `python -m
benchmarks.rust_hts125_pilot.build --help`. Submission is explicit via
`python -m benchmarks.rust_hts125_pilot.run prepare --help`; review can be repeated
readably days later using `collect ROOT`. Never overwrite accepted timings or
an existing figure folder. Only a convincing, strictly equivalent winner would
justify a separate 12-WES validation or full manuscript rerun decision.
