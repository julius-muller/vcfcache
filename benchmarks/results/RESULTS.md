# Results

Final campaign, 2026-08-24. Every cell is one genome, one assay, direct and
cached, serial invocation. Medians use up to 12 genomes; `n` shows the accepted
paired cells after explicit concordance exclusions.

Source: `campaign_metrics.tsv` (360 rows). The validated, machine-generated
outputs are `campaign_summary.tsv` and `concordance_failures.tsv`; regenerate
them with `../summarize_final_campaign.py`. The canonical archive index is
`raw/final_campaign/RUN_MANIFEST.json`.

## Break-even by annotator and assay

| tool | level | assembly | assay | n | variants | direct s | cached s | speedup | hit % |
|---|---|---|---|---|---|---|---|---|---|
| fastVEP | clinical | GRCh38 | panel | 12 | 255 | 10.5 | 22.4 | **0.46** | 90.4 |
| fastVEP | clinical | GRCh38 | wes | 12 | 80,190 | 48.6 | 50.7 | **0.93** | 93.2 |
| fastVEP | clinical | GRCh38 | wgs | 12 | 4,694,021 | 463.6 | 183.4 | **2.54** | 92.6 |
| fastVEP | clinical | GRCh37 | panel | 12 | 228 | 5.1 | 12.6 | **0.41** | 90.7 |
| fastVEP | clinical | GRCh37 | wes | 12 | 74,990 | 37.4 | 35.5 | **1.02** | 94.6 |
| fastVEP | clinical | GRCh37 | wgs | 12 | 4,196,713 | 216.4 | 114.9 | **1.70** | 95.8 |
| SnpEff | core | GRCh37 | panel | 12 | 228 | 33.5 | 38.3 | **1.07** | 90.7 |
| SnpEff | core | GRCh37 | wes | 11 | 69,911 | 67.4 | 57.1 | **1.16** | 94.7 |
| SnpEff | core | GRCh37 | wgs | 11 | 3,749,108 | 868.6 | 161.8 | **4.46** | 95.9 |
| SnpEff | clinical | GRCh38 | panel | 12 | 255 | 103.1 | 95.7 | **1.06** | 90.4 |
| SnpEff | clinical | GRCh38 | wes | 12 | 80,190 | 493.1 | 185.2 | **2.66** | 93.2 |
| SnpEff | clinical | GRCh38 | wgs | 12 | 4,694,021 | 8765.1 | 820.0 | **10.78** | 92.6 |
| SnpEff | clinical | GRCh37 | panel | 12 | 228 | 35.8 | 38.6 | **1.03** | 90.7 |
| SnpEff | clinical | GRCh37 | wes | 11 | 69,911 | 296.7 | 72.9 | **3.62** | 94.7 |
| SnpEff | clinical | GRCh37 | wgs | 11 | 3,749,108 | 3827.1 | 419.2 | **12.31** | 95.9 |

## What it shows

**Caching returns per-variant work, never start-up.** Panels are start-up
dominated: fastVEP is slower with caching, while SnpEff medians are effectively
at break-even (1.03--1.07x) and should not be interpreted as meaningful panel
acceleration. Every annotator wins at whole genome. Where it crosses over is set
by how much per-variant work the recipe does, not by the assay label.

- **fastVEP** (Rust, low per-variant cost) is still at 0.93-1.02x on exomes and
  only reaches 1.70-2.54x at whole genome.
- **SnpEff core** crosses just above 1.0 on exomes and reaches 4.46x at WGS.
- **SnpEff clinical** — the same variants, the same cache, plus ClinVar and
  gnomAD — already wins at exome scale (2.66-3.62x) and reaches 10.78-12.31x.

The clinical-versus-core pair on GRCh37 is the cleanest evidence: identical
genomes, identical cache membership, hit rates identical to the decimal
(94.7% / 95.9%). Only the annotator's per-variant cost changed, and the exome
speedup moves 1.16x -> 3.62x, the genome 4.46x -> 12.31x.

**Assembly is a cost axis for fastVEP, not for SnpEff.** `GRCh37.75` (Ensembl 75)
carries roughly a third of `GRCh38.115`'s transcripts. fastVEP's runtime tracks
that; SnpEff's does not, because its per-variant cost is dominated by fixed
per-record overhead rather than per-transcript work. So a lower fastVEP speedup
on the GRCh37 PGP genomes means less work to return, not a cohort effect.

**Hit rates are a property of the genomes, not the tool.** 90.4-95.9% across
three cohorts, two assemblies and two providers. PGP (GRCh37, independent
provider) lands in the same band as KPGP and SGDP, which argues the rate is
biological rather than an artefact of sample selection.

## Concordance

176 of 180 cells pass: cached output identical to direct for the fields defining
the tested annotation level, compared with a contig-order-independent digest
that canonicalizes record order within each chromosome-position block.

The four excluded cells are one genome, `hu6664F3`, at WES and WGS scale in
both SnpEff levels. Its two panel cells and all three fastVEP cells pass, and all
record counts match. Retained-output diagnosis found exactly one differing WES
allele, `chr2:11273459 C>T`, at both levels: cache-build SnpEff adds
`WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS` to two `C2orf50` consequences while
whole-exome SnpEff omits it. Consequence, transcript, HGVS, ClinVar and gnomAD
values otherwise agree. A direct-SnpEff-only two-record reproducer proved that
processing `chr2:11273458 A>AC` immediately before the target adds the warning,
while reversing the records removes it. Thus this is an upstream SnpEff 5.4c
processing-order side effect, not a VCFcache merge error. It is still an output
difference: no general warning exception was introduced and the four cells
remain excluded from the medians above. Exact inputs, outputs and values are
stored under `raw/diagnostics/hu6664F3/`.

The final GRCh37 arm uses the 12 samples frozen in
`../reproduction/final_clinical_campaign/samples.tsv`. It differs by two genomes
from the older 52-genome external campaign, so its 95.8% median hit rate must not
be substituted for the older campaign's 92.7% PGP estimate.

Three earlier campaign rounds were invalidated by harness defects, not by
VCFcache: stale output files across assay iterations with an unchecked return
code, duplicate input records in five PGP genomes, and a digest that assumed a
shared contig order. See `GUIDELINES.md`; the user-facing differences between
cached and direct output are recorded in `/KNOWN_LIMITATIONS.md`.
