# GRCh38 + fastVEP 0.5 recipe

A recipe for building and using a vcfcache cache with fastVEP 0.5 and its
supplementary annotation sources (`--sa-dir`).

## Why this recipe uses two source folders

fastVEP reads each source (`.osa2`) in blocks. `fastvep sa-build --chunk-bits 16`
(65 kb blocks) is what fastVEP advises, and it suits dense input: a whole genome
touches nearly every block, so each block is read once and used for many variants.

vcfcache sends fastVEP only the variants that are **not** in the cache. These are
few (for a whole genome, about 8% against the gnomAD AF >= 0.1% blueprint and 3%
against a cohort cache), and they are spread over the whole genome. With 65 kb
blocks they still open nearly every block of every source, so a cached run pays
almost the full source cost for a small fraction of the variants. The same
sources built with `--chunk-bits 10` (1 kb blocks) let the misses read only the
blocks they need. The output is identical: block width changes how fastVEP reads
a source, not what it returns.

Measured on three whole genomes (KPGP, 4.6-4.8 M variants), fastVEP 0.5.0 with 11
or 15 sources and two custom tracks, 3 repeats, outputs identical to direct
annotation at 65 kb:

| Cache | Misses at 65 kb | Misses at 4 kb | Misses at 1 kb |
|---|---|---|---|
| gnomAD AF >= 0.1%, 11 sources | 2.03x | 2.52x | 3.26x |
| gnomAD AF >= 0.1%, 15 sources | 1.80x | 2.38x | 3.37x |
| cohort, 11 sources | 1.73x | 3.40x | 4.23x |
| cohort, 15 sources | 1.64x | 3.48x | 4.30x |

Only with 1 kb sources for the misses does every added source raise the
speed-up. Narrower blocks are slower again (on chr1, 512 bp took 20-35% longer
than 1 kb and 256 bp about twice as long, with files 1.5 and 2.5 times larger).

Why not use the 1 kb sources for everything? Dense input reads them more slowly:

| Whole genome, direct | 65 kb | 1 kb |
|---|---|---|
| 11 sources | 676 s | 922 s |
| 15 sources | 992 s | 1,428 s |
| Cache build (18.9 M records), 15 sources | 1,281 s | 2,021 s |

So build the cache from the 65 kb sources and annotate the misses from the 1 kb
sources. The two folders take about 2.5 times the disk of one (18 + 28 GB for
15 sources).

**Use the 1 kb folder for genomes and exomes, not for panels.** Loading the
larger 1 kb indexes adds about 10-15 s to every run. That pays off when the
misses are spread over the genome, and costs time when they are few or local:

| Input (3 samples, median) | Direct | 0.1% cache, misses at 65 kb / 1 kb | cohort cache, 65 kb / 1 kb |
|---|---|---|---|
| Exome (Twist core, 74-81 k variants), 11 sources | 224-252 s | 3.60x / 3.79x | 3.53x / 4.96x |
| Exome, 15 sources | 330-358 s | 3.76x / 3.68x | 3.43x / 4.08x |
| Panel (ACMG SF v3.3, 207-289 variants), 11 sources | 6.2-6.6 s | 1.05x / 0.43x | 1.11x / 0.43x |
| Panel, 15 sources | 6.7-7.8 s | 1.12x / 0.34x | 1.23x / 0.36x |

A single chromosome behaves like a panel: chr22 of one genome (7,400 misses in
one chromosome) took 18.0 s at 1 kb and 12.4 s at 65 kb. For panels, single
chromosomes and other small or local inputs, annotate with `params.yaml`.

## Build the sources

Build every source twice, with the same file names, into two folders:

```bash
for bits in 16 10; do
  d=/data/fastvep/sa/GRCh38_$([ $bits = 16 ] && echo 65kb || echo 1kb)
  fastvep sa-build --assembly GRCh38 --source clinvar -i clinvar.vcf.gz -o $d/clinvar --chunk-bits $bits
  fastvep sa-build --assembly GRCh38 --source gnomad -i gnomad.exomes.v4.1.sites.chr1.vcf.bgz -o $d/gnomad_chr1 --chunk-bits $bits
  # ... every other source, identically
done
```

- Built-in sources (gnomad, dbsnp, clinvar, onekg, gerp, phylop, ...) may be built
  per chromosome.
- A custom source (`custom_vcf`, `custom_bed`) must be **one genome-wide file**:
  fastVEP treats every custom file as a source of its own, even with the same
  `--name` (a second file becomes `FV_<NAME>_2`, and so on), so per-chromosome
  custom files would make every variant look up and print the source once per file.
- Both folders must hold the same files; only `--chunk-bits` differs.

## Use

Edit the paths in `params.yaml` (cache build, 65 kb sources) and
`params.annotate.yaml` (genomes and exomes, 1 kb sources), then:

```bash
vcfcache cache-build --db <blueprint_dir> --name fastvep0.5_grch38 \
  -a vcfcache/recipes/grch38_fastvep0.5/annotation.yaml \
  -y vcfcache/recipes/grch38_fastvep0.5/params.yaml

vcfcache annotate -a <blueprint_dir>/cache/fastvep0.5_grch38 \
  -i sample.vcf.gz -o sample.annotated.bcf \
  -y vcfcache/recipes/grch38_fastvep0.5/params.annotate.yaml
```

`annotation.yaml` names the source folder only as `${params.sa_dir}`, so the
cache's recipe stays the same whichever folder a run uses.
