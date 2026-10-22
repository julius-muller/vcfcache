# Benchmark record

The evidence base for the VCFcache manuscript: what was measured, when, whether
it is usable, and where the untouched output lives.

Nothing here is code. Scripts stay in `benchmarks/`, cluster-side scripts stay on
the cluster, and this directory holds only results, logs and the three documents
that make them interpretable:

| file | what it is |
|---|---|
| `RUN_LOG.md` | every cluster job run, in order, with its outcome |
| `RESULTS.md` | **the findings** — final campaign table and what it shows |
| `campaign_metrics.tsv` | all 360 campaign rows, one per genome/assay/condition |
| `campaign_summary.tsv` | validated medians generated from accepted direct/cached pairs |
| `concordance_failures.tsv` | the four explicitly excluded cells and their digests |
| `README.md` | this file: goals, index, and the history log |
| `GUIDELINES.md` | how a benchmark is run and what has to be stored |
| `PLAN.md` | the annotation-configuration matrix and what is still outstanding |
| `raw/` | results and provenance: metrics, manifests, timings, cluster config |

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
  SnpEff 5.4c. vep-rs 0.4.0 has been piloted but not campaigned.
- **Assays** — ACMG SF v3.3 panel, Twist core exome, whole genome.
- **Cohorts** — KPGP (GRCh38), SGDP (GRCh38), PGP (GRCh37). All held out from
  the gnomAD AF>=1% blueprint the caches are built over.
- **Configurations** — VEP `--everything`, fastVEP core and supplementary
  databases, and SnpEff core/clinical. They are realistic native configurations,
  not identical cross-tool outputs. Every reported speed-up compares direct and
  cached execution of the same tool and recipe.

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
| vep-rs 0.4.0 capability and cost pilot | `raw/archive/veprs-pilot-grch38-chr22-diagnostic/` |
| vep-rs annotation-heavy ladder, genome-wide gnomAD | `raw/archive/veprs-heavy-grch38-chr22-genomewide/` |
| Fused-routing engine threshold check | `raw/archive/fused-routing-threshold-check/` |
| VCFcache fixed-overhead profile, panel/WES | `raw/archive/vcfcache-fixed-overhead-profile/` |
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
| 2026-08-22 | fastVEP clinical annotations | ClinVar and gnomAD `.osa2` built for **both** assemblies | primary |
| 2026-08-23 | Campaigns v1 (jobs 1003-1023) | Invalidated: stale outputs across assays, duplicate PGP inputs, contig-order digest | **excluded** |
| 2026-08-24 | **Final campaign** (jobs 1102-1165) | 180 cells, 12 genomes x 3 assays x 2 conditions x 5 tool/level/assembly arms; 176/180 concordant. See `RESULTS.md` | primary |
| 2026-08-22 | fastVEP clinical gates | GRCh38 and GRCh37 both passed on CSQ + FV_CLINVAR + FV_GNOMAD | primary |
| 2026-08-22 | `fastvep_clinical_grch37_af001` | Genome-wide cache, 18,670,412 alleles, 625 MB, built in 6:05 | primary |
| 2026-08-22 | `fastvep_clinical_grch38_af001` | Genome-wide cache, 18,869,857 alleles, 1.9 GB, built in 28:49 | primary |
| 2026-08-22 | `snpeff_core_grch37_af001` | Genome-wide cache, 18,670,412 alleles, 682 MB, built in 3:01:38 | primary |
| 2026-08-22 | `snpeff_clinical_grch37_af001` | Genome-wide cache, 18,670,412 alleles, 721 MB, built in 7:20:53 | primary |
| 2026-08-22 | `snpeff_clinical_grch38_af001` | Genome-wide cache, 18,869,857 alleles, 1.5 GB, built in 13:11:20 | primary |
| 2026-10-22 | vep-rs 0.4.0 pilot (job 1242) | Perl->JSON cache conversion works; clinical level unreachable; 5.68-14.00 us/variant, so the break-even rule predicts caching loses at every scale | diagnostic |
| 2026-10-22 | vep-rs heavy ladder (job 1251) | 66 cells; clinical level IS reachable via the plugin path; one genome-wide gnomAD source takes 7.42 -> 31.95 us/variant and crosses break-even | primary |
| 2026-10-22 | fused-routing threshold check | VCFcache per-variant overhead 11.03 -> 5.55 us (1.99x), start-up flat; threshold ~29.7 -> ~15.0 us/variant | diagnostic |
| 2026-10-22 | fixed-overhead profile (jobs 1253, 1254) | Panel loss is a full-cache pass: 23.6 s fixed, ~100% of cost at panel scale. Region-restricted routing worth 2x; fused regresses at 1 and 22 contigs | primary |

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
| SnpEff core cache build, GRCh37, 18.67 M alleles | 3:01:38 | `raw/snpeff/results/build_time_snpeff_core_grch37_af001.txt` |
| fastVEP clinical cache build, GRCh37 | 6:05 | `raw/richness/results/build_time_fastvep_clinical_grch37_af001.txt` |
| fastVEP clinical cache build, GRCh38 | 28:49 | `raw/richness/results/build_time_fastvep_clinical_grch38_af001.txt` |
| Annotations per record, chr22, GRCh37 / GRCh38 | fastVEP 3.19 / 12.14; SnpEff 4.27 / 12.29 | `bcftools query` over the built caches |
| SnpEff clinical cache build, GRCh37 / GRCh38 | 7:20:53 / 13:11:20 | `raw/snpeff/results/build_time_snpeff_clinical_*.txt` |
| Clinical content, chr22 (CLNSIG / gnomAD AF) | SnpEff 3061 / 27581 (GRCh37), 3069 / 28591 (GRCh38); fastVEP 3437 / 28904 (GRCh37) | `bcftools query` over the built caches |
| vep-rs per-variant cost, GRCh38, `--fork 8` | 5.68 us (core) to 14.00 us (`--everything`) | `raw/archive/veprs-pilot-grch38-chr22-diagnostic/pilot_1242.tsv` |
| vep-rs start-up | 0.16-0.25 s (transcript index load) | same, two-point fit |
| Per-variant cost needed to reach break-even at f=0.9258 | 29.7 us | 27.5 us VCFcache overhead / f |
| vep-rs + one genome-wide gnomAD source | 31.95 us/variant | `raw/archive/veprs-heavy-grch38-chr22-genomewide/fit_ladder.txt` |
| vep-rs clinical level (ClinVar + genome-wide gnomAD AF) | 40.76 us/variant | same |
| VCFcache per-variant overhead, fused-routing engine | 5.55 us (1.99x below legacy) | `raw/archive/fused-routing-threshold-check/fit.txt` |
| Break-even threshold, fused-routing engine | ~15.0 us/variant | same, scaled from the synthetic harness |
| VCFcache fixed overhead, legacy, measured at N=1 | 23.6-24.0 s | `raw/archive/vcfcache-fixed-overhead-profile/fit_overhead.txt` |
| Fixed-term share of VCFcache cost, panel / WES | ~100% / 87.5% | same |
| Region-restricted routing, 2-4 contigs vs legacy | 10.2 s vs 19.5 s | `raw/archive/vcfcache-fixed-overhead-profile/span_1254/` |
