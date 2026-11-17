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
- The comparison is always against the annotator as its authors ship or recommend it. Never change an annotator setting to make it slower or to favour vcfcache (e.g. fastVEP source block width); only vcfcache's own side may be tuned.
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
- Mission stack, whole genome (2026-11-17, `wgs_mission.sh`; 11 sources incl. dbSNP, gnomAD, SpliceAI, PhyloP, AlphaMissense + 2 mock tracks of ~100 M rows; 3 KPGP samples, 3 repeats, sl-w1/3/5, 16 cores; output identical to uncached record for record in every run). Median speedup over uncached fastVEP (~705-730 s):
  - AF >= 0.1% blueprint: 2.0-2.15x (8.2% misses); with narrow-block sources for the misses only: 2.2-2.4x.
  - Cohort cache (blueprint + the other 22 KPGP samples, 2.9-3.1% misses): 1.7x; with narrow sources for the misses: 3.2-3.65x.
  - AF >= 0.01% blueprint (one sample): 1.56x, narrow 3.16x. The user considers it the least important (hard to roll out on Zenodo, little gain over 0.1% or a cohort cache); it runs only with `AF0001=1`.
  - Rerun of a sample already in the cache (no misses): 10.6-11.7x.
  - With few misses the annotator still opens nearly every 65 kb block of every source, so narrow blocks for the miss step are what turns a lower miss rate into time saved. The uncached baseline always uses fastVEP's advised 65 kb blocks.
- Scripts: `benchmarks/gather_engine/cost_model/`; per-run tables in `benchmarks/gather_engine/cost_model/results/`; everything else on the cluster under `/results/fastvep05/cost_model/`.

## Known annotator issues (do not investigate)
- Ensembl VEP reports `HGNC_ID` inconsistently for a few genes: the field is filled or empty for the same variant depending on which other variants are in the run. This is a known upstream bug (https://github.com/Ensembl/ensembl-vep/issues/1959), still present in 116.3; a fix exists but is unmerged (https://github.com/Ensembl/ensembl-vep/pull/2056).
- Cached-versus-direct comparisons for VEP therefore show a handful of `HGNC_ID`-only differences (seen for FGF7P3, GUSBP5, PRAMENP, NMRAL2P, LINC03025). Ignore them: they are not a vcfcache defect and need no further analysis. Report them as this known VEP bug.

## Git workflow
- After completing and verifying requested repository changes, commit them and push the current branch unless the user explicitly asks not to.
- Author commits as `Julius Müller <julius.mueller@dkfz-heidelberg.de>` rather than with the VM's default `ubuntu` identity.
- Inspect the worktree before staging, preserve unrelated user changes, and describe the committed scope in the final response.
