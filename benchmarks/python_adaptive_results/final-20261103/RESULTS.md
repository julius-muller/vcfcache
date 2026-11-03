# Final Python-adaptive manuscript results

## Decision-ready summary

The Python/bcftools implementation delivers strong acceleration for annotation-rich
WGS workflows. It does **not** accelerate every annotator/assay combination. In
particular, this cohort does not support a claim of improved vep-rs core WES or WGS
runtime, or a consistent fastVEP clinical WES benefit.

Median paired speedups for the GRCh38 assay set (12 independent samples per cell;
six KPGP and six SGDP, each relative to its own recipe's direct baseline):

| Recipe | Panel | WES | WGS | Median WGS time saved |
|---|---:|---:|---:|---:|
| VEP everything | 2.34× | 2.35× | 9.50× | 2.29 h |
| fastVEP clinical | 1.00× | 0.88× | 2.44× | 4.65 min |
| SnpEff clinical | 1.18× | 2.72× | 9.95× | 2.36 h |
| vep-rs core | 1.18× | 0.46× | 0.89× | −9.92 s |

Ratios below 1 mean caching was slower. These are medians of per-sample ratios,
not ratios of median runtimes. Median speedup and median time saving need not
describe the same sample. The full distributions and cohort-stratified bootstrap
intervals are in Figure 2 and the TSV summaries; do not interpret the table as
a cross-tool comparison of biological annotation content.

## Completion and equivalence

All 600 canonical conditions reached terminal checkpoints at 07:47 CET on
2026-11-03. There were no missing conditions or operational failures:

- **572 accepted conditions**, retaining their original comparator scope.
- **24 clinical vep-rs WES/WGS timings** with known CSQ differences. Non-CSQ records
  and relevant headers passed comparison, but these are **not validated speedups**.
  The 12 clinical vep-rs panel conditions passed and remain in the full assay data.
- **Four excluded SnpEff conditions**: hu6664F3, WES/WGS, core/clinical. Each failed
  strict comparison at chr2:11273459 because of
  `WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS` in ANN. These reproduce the previously
  observed discrepancy; they are not execution/truncation failures and are not
  silently accepted.

Acceptance is not universal byte identity with direct annotation: inherited VEP
HGNC_ID handling and unordered DOMAINS/CSQ equivalence remain documented in the
archived comparator evidence. Of the 600 conditions, 336 use complete retained
cached-reference records with an annotator-specific direct comparator; 252 use
retained paired direct controls with strict complete-record/header comparison;
12 synthetic-delay conditions use composed VEP semantic references rather than
retained per-delay complete outputs. No new field exceptions were introduced.

## Dataset and execution

- Frozen product: `ee2162473` on `feature/python-adaptive`; Python orchestrates
  bcftools with no compiled routing helper. The installed distribution-version
  string is not the build identity; use the frozen commit and source hashes.
- Eight configured CPUs per timed condition, an exclusive worker allocation and
  a new Slurm step per condition; aggregate CPU and peak-memory metrics are saved.
- No technical repeats. All **252 paired direct controls were reused** after
  input/recipe/resource/CPU and artifact-hash binding; other approved historical
  direct baselines retain their recorded validation scope. **Zero direct reruns**.
- 600 executed measurements yield **612 plotting views**. Twelve WGS aliases
  reuse already measured conditions; they are not extra measurements. The WGS
  strategy figure covers all **52 genomes** (20 KPGP, 20 SGDP, 12 PGP), three cache
  strategies and two annotator recipes, giving 312 views. No 1000 Genomes cohort
  is included. Cohort-cache seed-ID separation is not proof of absence from gnomAD.
- The main assay figure uses GRCh38. PGP/GRCh37 is retained in the assembly-separated
  supplementary figure, because identically named recipes across assemblies need
  not add identical information or perform identical work.
- Runtime excludes cache preparation/build/download. The throughput display is
  `86400 / observed per-sample elapsed time`: serial-processing capacity at the
  measured allocation, not an experiment running a multi-sample service for 24 h.

## Pipeline-complexity supplement

The 12 synthetic-delay conditions use one real WGS genome and two caches, with
six configured delay levels each. The delay plugin models additional computation;
it adds no biological annotation and is not a heavy fastVEP or vep-rs experiment.
At the largest configured wait, the observed ratios are 5.05× (AF ≥ 10%) and
11.63× (AF ≥ 1%), saving 18.21 and 20.75 hours respectively. The direct timings are
reused, and per-delay outputs are validated against a composed semantic reference.
The curves are descriptive; no biological replication or constant-hit-rate
asymptote is inferred from this one sample.

## Files and reproducibility

- [Figures and captions](CAPTIONS.md): final output directory is **figures-v2/**.
  Figure 1 is the unchanged archived workflow; all empirical plots are regenerated
  from the Python replay, not relabelled Rust measurements.
- [Machine-readable measurements](tables/cells.tsv), [plotting views](tables/plot_views.tsv),
  [excluded conditions](tables/exclusions.tsv), [primary summary](tables/primary_grch38_summary.tsv).
- [RStudio data](tables/manuscript_data.rds): named list of formatted tibbles.
- [Snapshot](source/SNAPSHOT.json) binds the complete compact terminal evidence
  archive by SHA256. Raw metadata remain locally in source/campaign/; the tracked
  archive contains those originals, including all task, comparison and resource
  documents. Large BCFs and sampling streams were not downloaded.
- `previous_rust_cached_seconds` is a historical Rust-adaptive cached measurement,
  **not** a former-main measurement. Date/node/filesystem-cache differences prevent
  interpreting that column alone as a randomized engine comparison.

To render into another **new** directory:

```bash
Rscript --vanilla benchmarks/python_adaptive_results/render.R \
  benchmarks/python_adaptive_results/final-20261103 figures-v3
```

All plots use R/ggplot2, vector PDF/SVG and 600-dpi PNG. Multi-panel figures also
have separately exported panels. The excluded clinical timings are diagnostic
only. No old figure or manuscript prose was overwritten. Phase-profile and
build-amortization plots were not transplanted: this replay did not measure the
necessary phases or new matched build costs. The next editorial step is to choose
the lead panels and update manuscript numbers/citations to this frozen dataset.
