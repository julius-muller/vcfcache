# Benchmark record

The evidence base for the VCFcache manuscript: what was measured, when, whether
it is usable, and where the untouched output lives.

Nothing here is code. Scripts stay in `benchmarks/`, cluster-side scripts stay on
the cluster, and this directory holds only results, logs and the three documents
that make them interpretable:

| file | what it is |
|---|---|
| `README.md` | this file: goals, index, and the history log |
| `GUIDELINES.md` | how a benchmark is run and what has to be stored |
| `PLAN.md` | the annotation-configuration matrix and what is still outstanding |
| `raw/` | tool output exactly as written, never edited |

## Goals

1. **Establish when caching helps and when it does not.** Caching returns an
   annotator's per-variant work and never its start-up cost, so the useful claim
   is a break-even rule, not a single speedup number.
2. **Test that rule across a real range of per-variant cost**, using annotator
   configurations a practitioner would recognise rather than minimal ones chosen
   to look fast.
3. **Show the result is not an artefact of one cohort, one assembly, or one
   annotator.** Three cohorts, two assemblies, three annotators.
4. **Keep every number traceable** to a file produced by the tool that made it.

## What is being measured

- **Annotators** — VEP 115.2 (`--everything --hgvsg`), fastVEP 0.3.0,
  SnpEff 5.4c.
- **Assays** — ACMG SF v3.3 panel, Twist core exome, whole genome.
- **Cohorts** — KPGP (GRCh38), SGDP (GRCh38), PGP (GRCh37). All held out from
  the gnomAD AF>=1% blueprint the caches are built over.
- **Levels** — core (consequences only) and clinical (core + ClinVar + gnomAD
  allele frequencies). VEP's `--everything` cache already sits at the clinical
  level, which is what fixes the target for the other two.

GIAB genomes and the 1000 Genomes primary cohort appear in the early entries
below. **GIAB was dropped** and no GIAB result supports any manuscript claim;
the 1000 Genomes campaigns were superseded by the held-out external cohorts,
which do not share samples with gnomAD.

## Index of usable results

| what | where |
|---|---|
| SnpEff core GRCh38, WGS, 12 genomes | `raw/archive/snpeff-core-grch38-wgs-serial/` |
| SnpEff core GRCh38, Panel + WES | `raw/archive/snpeff-core-grch38-panel-wes-split/` |
| fastVEP + VEP external cohorts, 52 genomes | `raw/campaigns/external-*/` |
| Panel/WES assay extension (warm control) | `raw/campaigns/external-assay-v3-independent-af001-warm/` |
| fastVEP annotation-richness ladder | `raw/richness/results/` |
| Controlled per-variant cost sweep | `raw/campaigns/controlled-runtime-light-*/` |
| Cluster and tool versions at collection | `raw/COLLECTION_ENV.txt` |

Figure inputs and the derived tables the manuscript cites remain where the
figure pipeline expects them, under `benchmarks/figures/` and
`benchmarks/external_assay_v3/source_data/`. Every claim in
`drafts/gigascience_technical_note/claims.tsv` names its own evidence file.

## History log

One line per run that produced usable data. `superseded` means the measurements
are sound but a better-designed run replaced them; `excluded` means the design
was flawed and the numbers do not estimate what they appear to.

| date | run | outcome | status |
|---|---|---|---|
| 2026-07-31 | `primary-wgs-9f02f75fdc54` | First 1000 Genomes WGS campaign, 50 genomes, 150 rows | superseded |
| 2026-07-31 | `primary-wgs-25eec5863cd4` | Repeat under a corrected cache | superseded |
| 2026-08-01 | `primary-wgs-271155276f68` | Repeat with HGNC-aware semantic review | superseded |
| 2026-08-02 | `primary-wgs-32a760323fc5-bundled-af1` | 1000 Genomes WGS over the bundled AF>=1% cache | superseded |
| 2026-08-03 | `assay-singlepass-32a760323fc5` | Panel and WES, 20 genomes, single-pass design | superseded |
| 2026-08-04 | `external-wgs-6c422804f208-hg19` | **52 held-out genomes** (KPGP, SGDP, PGP), WGS, VEP | primary |
| 2026-08-06 | `external-fastvep-1089a2f58cc4` | First fastVEP arm over the same 52 genomes | superseded |
| 2026-08-06 | `external-fastvep-25fe9b44f91a` | fastVEP repeat after input repair | superseded |
| 2026-08-06 | `external-fastvep-804963323e11-light-strict` | fastVEP under the strict light recipe | superseded |
| 2026-08-06 | `external-fastvep-2cd2ab56d8ee-light-strict-v2` | **fastVEP, 52 genomes, matched recipe** | primary |
| 2026-08-06 | `external-vep-light-calibration-bbef41d7db8e` | VEP light-recipe calibration, 6 genomes | primary |
| 2026-08-07 | `external-vep-light-full52-bbef41d7db8e` | **VEP light recipe, 52 genomes** | primary |
| 2026-08-10 | `controlled-runtime-light-14ad84d-50-90` | Controlled per-variant cost sweep, first pass | superseded |
| 2026-08-10 | `controlled-runtime-light-cc100c0-50-90` | **Controlled sweep via `SyntheticDelay.pm`** | primary |
| 2026-08-10 | `wgs-pipeline-spectrum-4f854e5-kpgp00319` | WGS pipeline-complexity spectrum, single genome | superseded |
| 2026-08-11 | `wgs-pipeline-spectrum-af001-cd87186-kpgp00319` | **Pipeline spectrum over the AF>=1% cache** | primary |
| 2026-08-18 | `external-assay-v3-independent-af001` | Panel/WES extension; cold-state order confounded the short cells | **excluded** |
| 2026-08-18 | `external-assay-v3-independent-af001-warm` | **Panel/WES extension with an independent warm control** | primary |
| 2026-08-19 | fastVEP richness ladder | Per-variant cost ladder over supplementary annotations; `raw/richness/results/` | primary |
| 2026-08-20 | `snpeff-core-grch38-panel-wes-split` | SnpEff core, Panel + WES, 12 genomes, split invocation | superseded |
| 2026-08-21 | `snpeff-core-grch38-wgs-split` | SnpEff core WGS, 3 genomes, split invocation; retained as the optimised-invocation bound | superseded |
| 2026-08-21 | `snpeff-core-grch38-wgs-serial` | **SnpEff core WGS, 12 genomes, serial**; median 8.17x, median hit rate 92.58% | primary |
| 2026-08-22 | SnpEff clinical gates | GRCh38 and GRCh37 both passed: cached output identical to direct on ANN + CLNSIG + AF | primary |
| 2026-08-22 | SnpEff core GRCh37 gate | Passed; the newly installed GRCh37.75 database validated against VCFcache | primary |
| 2026-08-22 | fastVEP clinical annotations | ClinVar and gnomAD `.osa2` built for GRCh38 (GRCh37 in progress) | primary |

### Measured constants worth keeping in one place

| quantity | value | source |
|---|---|---|
| SnpEff core cache build, GRCh38, 18.87 M alleles | 2:55:28 | `raw/snpeff/results/snpeff_cache_build_time.txt` |
| SnpEff core per-variant cost | 0.558 ms | derived from the build above |
| SnpEff clinical throughput, GRCh38 | 327 records/s | `raw/snpeff/logs/` clinical gate |
| SnpEff clinical throughput, GRCh37 | 371 and 422 records/s | two runs on the same slice, 14% apart |
| SnpEff clinical per-variant cost, GRCh38 | ~2.37 ms | two-point fit, start-up separated |
| SnpEff clinical start-up | ~215 s | three JVM loads plus database reads |
| fastVEP cache overhead | 21.7 s + 27.5 us per variant | `benchmarks/analyze_breakeven.py` |
