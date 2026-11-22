# Repository Guidelines (untracked helper)

## Structure
- Code: `vcfcache/` (CLI, database orchestrators, utils, manifest helpers).
- Assets: `resources/`, Docker files in `docker/`, docs in `README.md`, `WIKI.md`, manifest in `public_caches.yaml`.
- Tests/fixtures: `tests/`, configs in `tests/config/`.

## Build/Test
- Dev env: `uv venv .venv && source .venv/bin/activate && uv pip install -e ".[dev]"`.
- Tests: `python -m pytest` (network/Zenodo integration gated by env).
- Lint/format/type: `ruff check --fix vcfcache`, `black vcfcache`, `mypy vcfcache`.
- Docker: multi-stage `docker/Dockerfile.vcfcache`; `--target test` installs dev deps and can run pytest; `--target final` is the lean runtime (`ghcr.io/julius-muller/vcfcache:<version>`).

## Manifest & Aliases
- Aliases: blueprint `bp-<genome>-<source>-<release>-<filt>`; cache `cache-<genome>-<source>-<release>-<filt>-<tool>-<tool_version>-<preset>`.
- Manifest entries include `version` (vcfcache used), no image tags. Fields: alias, type (blueprint|cache), genome/source/release/filt, doi, updated_at, md5; caches add tool/tool_version/preset and annotation_yaml_md5.
- CLI resolves aliases via `public_caches.yaml`; caches/blueprints live on Zenodo.

## Workflow tips
- Keep `annotation.yaml` immutable per cache; store `params.yaml` snapshots with caches.
- Record new aliases in `public_caches.yaml` when publishing.

## Docs
- `WIKI.md` is canonical for alias/manifest schema, Zenodo pull/push, Docker usage, cache layout. README mirrors high-level usage.

## Mission (manuscript benchmark)
- Goal: speed up the average annotation pipeline on WGS while output stays identical to direct annotation. The reference pipeline is an annotator with 10-20 well-known sources (gnomAD, dbSNP, ClinVar, REVEL, AlphaMissense, SpliceAI, conservation scores, ...) plus a few custom BED/bigWig-style tracks of about 100 million rows with a couple of columns each (e.g. in-house AI predictions).
- The comparison is always against the annotator as its authors ship or recommend it. Never change an annotator setting to make it slower or to favour vcfcache: the direct (uncached) baseline always runs as advised (e.g. fastVEP sources at 65 kb, `--chunk-bits 16`). vcfcache's own side may be tuned, and so may the miss step of a cached run: there a second copy of the same sources in narrower blocks is allowed (decided with the user, 2026-11-19; the fastVEP 0.5 recipe uses 1 kb), because it reads the same data and its output is identical.
- Out of scope: per-variant computation that few people run (e.g. fastVEP `--acmg`). Do not benchmark it.
- fastVEP 0.5.0 on WGS is the hard case and decides the project; a substantial speedup means 3x or more.

## What has been measured (2026-11, fastVEP 0.5.0, KPGP-00128, 4.6 M variants, 8 cores)
- Cached run = vcfcache pulls the hits + the annotator runs on the misses. It only wins if the annotator's run on the misses is much cheaper than its run on everything.
- vcfcache's own cost: about 17 s fixed (reading the 18.9 M-record cache) + 8.7 us wall per variant; 57 s for this genome.
- fastVEP per variant: about 20 us wall for base, `--hgvs`, `--everything`, ClinVar, gnomAD, REVEL; about 40 us with dbSNP and AlphaMissense added. Flags and sources add almost nothing per variant.
- fastVEP reads a source in blocks of genome and answers every variant in a block almost for free, so a source costs per block opened, not per variant. Cache misses (8.2% with the AF >= 0.1% gnomAD blueprint, 3.0% with AF >= 0.01%) are spread over the whole genome and open nearly every block, so the annotator's run on the misses costs almost as much as on all variants for the source part.
- Results: clinical stack (HGVS, ClinVar, gnomAD) 141 s uncached, 91 s cached (1.55x). Dense stack (adds dbSNP, AlphaMissense), extrapolated from chr1: about 1.3x. SpliceAI2 scores (5 columns, chr21 region): 3.5 s per 5 Mb whether 4,000 or 32,000 variants are annotated.
- fastVEP's generic `custom_vcf` source builder holds all records in memory (out of memory at 28 GB for 12.7 M records with 79 columns); dedicated readers stream.
- Record order: cached output equals direct output record for record since 0bd2030b3 (hits and misses sharing a position used to swap).
- Setup flaw found 2026-11-19 (affects the two fastVEP mission bullets below): fastVEP makes every `custom_vcf` file a source of its own, also when several files share one `--name` (fields `FV_TRACKA`, `FV_TRACKA_2`, ... `_24`, each repeating the value). The mock tracks were built per chromosome, so fastVEP looked up and printed each track 24 times; the uncached baseline was therefore too slow and those speedups are too high. A custom source must be one genome-wide file (`build_tracks1.sh`: 94.8 M rows, 5-7 min, 20.8 GB peak). gnomAD genomes as 166 custom pieces made it far worse (chr22 alone: 80 s instead of 21 s; whole genome over 40 min and 30,600 CPU-s vs 624 s and 2,400 CPU-s without it); it is dropped, since a genome-wide custom build of ~760 M records does not fit in memory. `wgs_mission.sh` now stops if any `FV_*_<n>` field appears. Corrected runs: jobs 3410-3415 (`mission[2]f_wgs_<sample>`), results below.
- Corrected mission runs (2026-11-19, `wgs_mission.sh`, one genome-wide file per mock track; 3 KPGP samples, 3 repeats, outputs identical to uncached in every run; median speedup over uncached):
  - 11 sources + 2 tracks, uncached 665-685 s: 0.1% 2.0-2.1x, 0.1%+narrow 2.3-2.5x, cohort 1.6-1.8x, cohort+narrow 3.1-3.4x, rerun 11.1-11.5x.
  - 15 sources + 2 tracks (adds 1000 Genomes, GERP, dbNSFP, ENCODE cCREs), uncached 911-996 s: 0.1% 1.75-1.83x, 0.1%+narrow 2.3-2.4x, cohort 1.55-1.73x, cohort+narrow 3.3-3.55x, rerun 13.6-15.2x.
  - Block width for the miss step (`diag_nb.sh`, chr1 of KPGP-00128, fastVEP alone, s, mean of 2; output identical at every width). All 369 k variants at 65 kb: 39.6 (11 sources) / 54.2 (15). 0.1% misses (38 k) at 65 kb / 4 kb / 1 kb / 256 bp: 25.6/19.6/15.8/29.6 (11) and 41.6/27.0/20.5/37.4 (15). Cohort misses (12.5 k): 26.0/11.3/11.8/25.2 (11) and 38.9/15.3/13.4/30.3 (15). At 1 kb the misses pay 4.7 s (0.1%) and 1.6 s (cohort) of the 14.6 s the four added sources cost on all variants, so with 1 kb miss sources each added source raises the speedup. 256 bp is too fine (13 GB, index-bound). Whole-genome 1 kb arms: builds 3434-3479, benchmarks 3480-3485 (`N10=1`, `results/mission[2]fn10_wgs_<sample>.tsv`).
  - Whole-genome result, miss step at 65 kb / 4 kb / 1 kb (3 samples x 3 repeats, median of per-sample median speedups; one digest per sample over all arms). 11 sources (uncached 669-694 s): 0.1% 2.03x/2.52x/3.26x, cohort 1.73x/3.40x/4.23x, rerun 11.2x. 15 sources (uncached 913-999 s): 0.1% 1.80x/2.38x/3.37x, cohort 1.64x/3.48x/4.30x, rerun 14.0x. Only at 1 kb do the four added sources raise the speedup in every arm: they add ~310 s uncached but ~90 s (0.1%, 205->300 s) and ~60 s (cohort, 161->225 s) cached, so the time saved grows from 470 to 680 s (0.1%) and 512 to 770 s (cohort). The ratio gains less than the saving, because misses still pay per-block cost and the dense added sources (GERP, 1000G) have little per-variant cost.
  - The added sources are dense (GERP: a score per base; 1000 Genomes): fastVEP's cost for them is mostly per block opened, which the miss run still pays at 65 kb (0.1% cache: +194 s of the +233 s they add uncached). Speedup grows with sources only where the miss run skips blocks (cohort+narrow) or there are no misses (rerun 11x -> 14x).
  - Widths below 1 kb (2026-11-20, `diag_nb.sh` with `BS="16 10 9 8" TAG=9`, chr1 of KPGP-00128, s, mean of 2; `results/diag_nb9_chr1_KPGP-00128.tsv`): 0.1% misses at 65 kb/1 kb/512 bp/256 bp 26.2/17.4/20.6/32.0 (11 sources) and 40.6/21.5/27.0/40.7 (15); cohort misses 28.4/12.0/17.6/29.8 and 42.2/14.5/18.2/34.8. 1 kb is the optimum; 15-source chr1 files 2.5/5.3/7.8/13 GB.
  - One source set or two (2026-11-20, `wgs_width.sh`, `results/width[2]_wgs_<sample>.tsv`, digests identical): direct whole genome at 65 kb vs 1 kb 676 vs 922 s (11 sources, KPGP-00128) and 992 vs 1428 s (15, KPGP-00006); cache build of the 18.9 M-record blueprint 785 vs 1298 s and 1281 vs 2021 s. Dense input must read 65 kb sources, so the recipe keeps two folders (13+19 GB, 18+28 GB).
  - Fixed cost of 1 kb sources (2026-11-20, `fixed_cost.sh`, 11 sources, random genome-wide subsets of KPGP-00128, s, 3 repeats; `results/fixed_cost_KPGP-00128.tsv`): 0 variants 5.3 vs 14.8 (65 kb vs 1 kb), 1 k 13-36 vs 14-28, 10 k 88 vs 47, 100 k 348 vs 102. 1 kb adds ~10 s per run (index load) and pays off from a few thousand scattered misses. Recipe test on chr22 alone (7,384 misses in one chromosome, `recipe_test.sh`, job 3505, identical digests): direct 14.8 s, misses at 65 kb 12.4 s, at 1 kb 18.0 s, so 1 kb is for whole genomes, not for small or single-chromosome inputs.
- Mission stack, whole genome (2026-11-17, superseded by the corrected runs; `wgs_mission.sh`; 11 sources incl. dbSNP, gnomAD, SpliceAI, PhyloP, AlphaMissense + 2 mock tracks of ~100 M rows; 3 KPGP samples, 3 repeats, sl-w1/3/5, 16 cores; output identical to uncached record for record in every run). Median speedup over uncached fastVEP (~705-730 s):
  - AF >= 0.1% blueprint: 2.0-2.15x (8.2% misses); with narrow-block sources for the misses only: 2.2-2.4x.
  - Cohort cache (blueprint + the other 22 KPGP samples, 2.9-3.1% misses): 1.7x; with narrow sources for the misses: 3.2-3.65x.
  - AF >= 0.01% blueprint (one sample): 1.56x, narrow 3.16x. The user considers it the least important (hard to roll out on Zenodo, little gain over 0.1% or a cohort cache); it runs only with `AF0001=1`.
  - Rerun of a sample already in the cache (no misses): 10.6-11.7x.
  - With few misses the annotator still opens nearly every 65 kb block of every source, so narrow blocks for the miss step are what turns a lower miss rate into time saved. The uncached baseline always uses fastVEP's advised 65 kb blocks.
  - Why the cohort cache is no faster than the 0.1% cache at 65 kb (`diag_cohort.sh`, KPGP-00128): fastVEP alone on the misses takes ~380 s for 376 k misses and for 132 k misses alike (more than half of its 690 s on all 4.6 M variants): the cost is opening source blocks genome-wide, not the variants. The cached run is bounded by that step (0.1%: 332/387 s, cohort: 387/379 s, equal within noise). With narrow sources the misses take 250-260 s (0.1%) and 150-155 s (cohort). The hit lookup alone: 70 s (18.9 M cache records), 75 s (23.8 M).
- vep-rs 0.5.0, same stack and samples (2026-11-18, `wgs_mission_veprs.sh`; vep-rs has no `--custom`, so the sources are its built-in plugins: generic dbNSFP-plugin tables for ClinVar, gnomAD exomes with ten groups, dbSNP, PhyloP and the two mock tracks; SpliceAI, REVEL, AlphaMissense, pLI, LoFtool, GWAS plugins; no RepeatMasker). Uncached 2,844-2,951 s. Median speedup, 3 repeats, outputs identical: 0.1% blueprint 6.3-6.6x, cohort 11.5-12.3x, rerun 52-58x. vep-rs looks up each variant in tabix tables, so its cost scales with the variants and few misses cost little.
  - Larger vep-rs stack (2026-11-20, `wgs_mission_veprs.sh STACK2=1`: adds 1000 Genomes AF table, dbNSFP SIFT4G/PolyPhen2 and GERP/phyloP via the dbNSFP plugin; jobs 3498-3503; `results/veprs_mission2_wgs_<sample>.tsv`, one digest per sample): uncached 5930-6472 s (vs 2844-2951 s with 11 sources); 0.1% 6.47-6.68x (vs 6.26-6.61x), cohort 12.8-14.0x (vs 11.5-12.3x), rerun 106-115x (vs 52-58x). vep-rs pays per variant, so added sources keep or raise the speedup without any tool-specific setting.
- Larger fastVEP stack (`wgs_mission_stack2.sh`, run with `STACK2=1`): adds 1000 Genomes 30x, GERP 91 mammals, dbNSFP 5.4a (fastVEP's reader keeps only SIFT/PolyPhen: 0.2 GB) and ENCODE cCREs: 15 known sources + 2 mock tracks. gnomAD genomes is built too but left out (see the setup flaw above). dbNSFP is academic (CC BY-NC-ND): no published cache from it; its download links are private files on the cluster.
- Cost-model figure and table (2026-11-22, `cost_model/figure.py`; `results/cost_model_table.{tsv,md}`, `results/cost_model_figure.{png,pdf}`; matplotlib needed only for the figure): per annotator x stack x cache x miss-step width, observed speedup D/C vs ideal D/(R + m D) (R = rerun time, m = miss fraction) and cost per miss (C - R)/(m D). vep-rs reaches 50-64% of ideal (cost per miss 1.7-2.2x); fastVEP 16-35% at 65 kb (cost per miss 5-18x) and 43-56% at 1 kb (2.6-5.4x). Panel B: time split into vcfcache's own cost and fastVEP on the misses, 15 sources.
- Scripts: `benchmarks/gather_engine/cost_model/`; per-run tables in `benchmarks/gather_engine/cost_model/results/`; everything else on the cluster under `/results/fastvep05/cost_model/`.

## Known annotator issues (do not investigate)
- Ensembl VEP reports `HGNC_ID` inconsistently for a few genes: the field is filled or empty for the same variant depending on which other variants are in the run. This is a known upstream bug (https://github.com/Ensembl/ensembl-vep/issues/1959), still present in 116.3; a fix exists but is unmerged (https://github.com/Ensembl/ensembl-vep/pull/2056).
- Cached-versus-direct comparisons for VEP therefore show a handful of `HGNC_ID`-only differences (seen for FGF7P3, GUSBP5, PRAMENP, NMRAL2P, LINC03025). Ignore them: they are not a vcfcache defect and need no further analysis. Report them as this known VEP bug.

## Git workflow
- After completing and verifying requested repository changes, commit them and push the current branch unless the user explicitly asks not to.
- Author commits as `Julius Müller <julius.mueller@dkfz-heidelberg.de>` rather than with the VM's default `ubuntu` identity.
- Inspect the worktree before staging, preserve unrelated user changes, and describe the committed scope in the final response.
