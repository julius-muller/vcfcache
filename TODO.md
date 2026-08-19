# TODO: like-for-like VEP and fastVEP complexity assessment

## Current result

fastVEP was faster for the tested GRCh37/hg19 PGP genomes than for the GRCh38/hg38 KPGP and SGDP genomes:

| Cohort | Assembly | Median variants | Median direct fastVEP runtime | Runtime per million variants |
|---|---:|---:|---:|---:|
| PGP | GRCh37 | 4.20 M | 207.3 s | 58.1 s |
| KPGP | GRCh38 | 4.74 M | 426.7 s | 88.2 s |
| SGDP | GRCh38 | 4.84 M | 401.8 s | 82.1 s |

This does **not** show that GRCh37 is intrinsically faster. The runs used different cohorts and substantially different transcript resources. In the frozen fastVEP resources, the autosomal GRCh38 annotation contains about 2.75 times as many transcript/mRNA records as GRCh37 (252,846 versus 91,845), 3.07 times as many exon records, and 3.16 times as many CDS records. The GRCh38 transcript cache is also about twice as large (207 MB versus 105 MB). Input size explains only part of the runtime difference; annotation complexity is a likely contributor, while cohort/provider pipelines and primary-assembly versus toplevel FASTA remain confounders.

Therefore:

- Do not interpret the current cross-assembly difference as an assembly effect.
- Keep PGP out of the matched main comparison figure.
- Existing within-assembly VEP and fastVEP benchmarks remain useful.
- Source measurements: `benchmarks/figures/source_data/final/2026-08-09T223718Z-light-matched-final/fastvep_external_wgs_metrics.tsv`.

## Why this is important for the paper

The complete 52-genome fastVEP campaign used a core consequence/HGVS recipe. It
answers the main question for a fast annotator under that configuration, but it
does not yet establish how fastVEP behaves in a production-style pipeline with
many supplementary databases, prediction scores, or other annotation sources.

One separate dense fastVEP pilot did add ten supplementary databases. On a
4.33-million-variant WGS, direct annotation took 707.4 seconds; VCFcache took
196.2 seconds at 90% hits (3.60x speed-up) and 82.6 seconds at 100% hits
(8.57x). This shows that caching remains strongly useful when fastVEP has more
lookup and output work. However, it is an engineering stress test, including
six dense deterministic custom sources, rather than a matched biological
VEP-versus-fastVEP comparison.

Published fastVEP genome-wide data likewise do not show VEP becoming faster
than fastVEP as annotations are added: from consequence-only to gnomAD plus
ClinVar, fastVEP increases from 197.84 to 218.54 seconds in the one-thread
measurement (10.5%), while VEP increases from 4,621 to 4,905 seconds (6.1%).
The reported advantage changes only from 23.4x to 22.4x. This may reflect
different fixed and per-database costs, but cannot answer whether the two tools
deliver comparable information content or whether additional fastVEP threads
eliminate the extra work. Source:
`https://github.com/Huang-lab/fastVEP/blob/master/manuscript/data/vep_genomewide.csv`.

This is therefore a paper-critical follow-up, not merely an optimization task:

- It prevents an unfair cross-tool claim based on nominal command complexity.
- It shows where VCFcache adds value for both current VEP pipelines and fast,
  Rust-based annotators likely to become more important.
- It lets the manuscript distinguish three separate benefits: choosing a fast
  annotator, parallelising that annotator, and avoiding repeat annotation of
  cached variants.
- It provides the evidence needed for a credible statement about future
  annotation pipelines rather than a result limited to a minimal fastVEP run.

## Priority audit

Assess the actual annotation workload and information added per variant for **VEP and fastVEP on both GRCh37/hg19 and GRCh38/hg38**, so comparisons are like for like.

1. Freeze and record the exact toolchain and resources for every condition:
   - annotator version/commit and complete command;
   - transcript annotation release and checksum;
   - FASTA and index;
   - supplementary databases and indexes;
   - configured threads, CPU affinity, storage, sorting, compression and indexing.

2. Characterize identical or closely matched input sets:
   - variants and alleles after normalization;
   - SNVs, indels and multiallelic sites;
   - variants overlapping genes, transcripts, exons and CDS regions;
   - where feasible, use the same biological variants across assemblies and document liftover losses.

3. Quantify information content added per output variant:
   - CSQ entries and transcript consequences;
   - unique transcripts and genes annotated;
   - populated INFO/FORMAT fields and supplementary annotations;
   - compressed and uncompressed output bytes;
   - non-empty annotation values, not merely declared header fields.

4. Quantify annotation work and lookups:
   - transcript interval lookups and candidate transcripts examined;
   - consequences emitted;
   - reference FASTA fetches and bases fetched;
   - ClinVar, gnomAD, 1000 Genomes, TOPMed and other database lookups, hits and misses;
   - VCFcache lookups, hits, misses and merge/output work.

   Add lightweight counters or use an instrumented build/profiler where an annotator does not expose these metrics. Counters must not be included in timed publication runs unless their overhead is shown to be negligible.

5. Run a controlled comparison matrix:
   - VEP and fastVEP;
   - GRCh37 and GRCh38;
   - a minimal/core recipe and a matched richer recipe;
   - the same supplementary annotation set on each assembly;
   - a common transcript universe where possible, alongside each tool's normal production resources.

6. Report normalized outcomes:
   - seconds per million input variants;
   - seconds per million transcript consequences emitted;
   - output bytes and populated annotation values per variant;
   - time attributable to parsing, lookups, consequence calculation, output construction, compression and indexing;
   - VCFcache overhead and speed-up at controlled cache-hit rates.

## Completion criterion

We can explain whether an observed runtime difference comes from the annotator, assembly, input cohort, transcript universe, supplementary lookups or output volume. Only then should the manuscript compare VEP and fastVEP across assemblies or claim that one assembly/configuration is faster.
