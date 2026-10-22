# vep-rs 0.4.0 pilot — what it establishes

Slurm job 1242, 2026-10-22, sl-w1, 16 vCPU. Total wall time about 100 s,
including the container pull and the cache conversion.

## 1. The cache converts, and it carries transcripts only

`storable_to_json.pl`, run inside `ensemblorg/ensembl-vep:release_116.2`,
converted chr22 of our Perl VEP 115 GRCh38 cache in 11 s: 51/51 files, 12,225
transcripts, 60 MB. vep-rs loads it and annotates.

It does **not** carry the co-located variant database. The Perl cache has it as
`all_vars.gz` per chromosome; the converter has no option to read it, and vep-rs
says so explicitly on every `--everything` run:

```
INFO vep_cli::json_cache: No variations directory found at .../GRCh38_chr22/variations
```

`veprs_field_probe.txt` confirms the consequence of that on real output. Over
231 missense consequences in a coding-dense region, **zero** carry SIFT,
PolyPhen, DOMAINS, `AF`, any `gnomAD*_AF`, `CLIN_SIG` or `Existing_variation`.
The long `## Extra column keys:` list in the output header is static
documentation and is not evidence that a field was populated — that mistake is
easy to make here and the probe exists to prevent it.

What a converted cache *does* reach: consequence and IMPACT, VARIANT_CLASS,
SYMBOL / SYMBOL_SOURCE / HGNC_ID, BIOTYPE, CANONICAL, TSL, CCDS, ENSP,
SWISSPROT / TREMBL, EXON / INTRON, and HGVSc / HGVSp. That is **L1 core plus
HGVS**. L2 clinical, as this project defines it (ClinVar + gnomAD AF), is not
reachable.

`--everything` also warns that `--regulatory` has no effect, so no TFBS or
regulatory-region consequence is produced. `--custom` and `--refseq` are
likewise accepted and unimplemented. A like-for-like reproduction of our VEP
`--everything --hgvsg` arm is therefore not possible, and vep-rs has to be
reported as its own native configuration.

## 2. Cost, from two sizes rather than one

16,385 and 314,314 records from the chr22 blueprint, `--fork 8`, fitting
`t = s + c·N`:

| config | per variant | start-up | throughput |
|---|---|---|---|
| core | 5.68 µs | 0.253 s | 154,187 rec/s |
| + symbol/canonical | 6.37 µs | 0.240 s | 140,245 rec/s |
| + HGVS | 10.59 µs | 0.201 s | 89,083 rec/s |
| `--everything` | 14.00 µs | 0.164 s | 68,842 rec/s |

Start-up is essentially the transcript index load, logged at about 170 ms. The
cache is read lazily per chromosome, so a genome-wide cache should not change it
much — but that should be measured, not assumed, before any absolute runtime is
quoted.

## 3. The break-even rule predicts caching loses here, at every scale

VCFcache's own overhead is 21.7 s + 27.5 µs per variant
(`benchmarks/analyze_breakeven.py`). With hit rate *f*:

```
T_direct − T_cached = f·N·t_ann − 21.7 − 27.5e−6·N
```

The per-variant coefficient is `f·t_ann − 27.5 µs`. At the measured *f* = 0.9258,
break-even needs `t_ann ≥ 29.7 µs`. vep-rs's **richest reachable** configuration
costs 14.00 µs — less than half the threshold. The coefficient is negative for
every configuration in the table, so the deficit grows with *N*: there is no
cohort or genome size at which caching turns profitable.

Projected at a 4.8 M-variant genome, `--everything`: direct 67 s, cached 159 s,
**0.42×**. Caching would roughly double the runtime.

Two caveats that do not change the direction. The 27.5 µs/variant figure was fit
against fastVEP; VCFcache's per-variant work is `bcftools` splitting, annotating
and merging, which is annotator-independent, but it deserves its own measurement
for this recipe. And vep-rs runs on 8 threads while the `bcftools` stages are
largely serial — which makes caching look *worse*, not better.

## 4. What this is worth

This is the first point in the matrix where the rule predicts a loss before the
measurement is taken, and predicts it for all N rather than below a threshold.
Confirming it is a stronger test of the rule than another win would be.

## Not quoting the authors' own comparisons

The README claims 187–343× faster than Perl VEP; the preprint abstract says
78–284×. The headline F1 of 1.000000000 is *adjusted*, excluding differences the
authors judge to be VEP defects; raw pooled F1 is 0.999998408, and structural
variants are much weaker (0.91 on GRCh38). Both authors are Natera employees.
None of that affects a direct-vs-cached measurement of vep-rs against itself,
which is the only comparison this project needs from it.
