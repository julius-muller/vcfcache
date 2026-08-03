# External WGS cohort similarity

This note records a representative variant-overlap check for the external WGS
benchmark inputs. It is a design diagnostic, not a population-genetics result:
one pair was sampled from each cohort and one pair across the two GRCh38
cohorts. The calculation was completed on 2026-08-03, before cache construction.

## Method

The first two lexicographically sorted evaluation source IDs from each cohort
were used for the within-cohort checks. The first KPGP and SGDP IDs were used
for the across-cohort GRCh38 check. Each prepared VCF was already split,
left-normalized, sorted, and restricted to autosomes 1-22.

Variants were compared as exact normalized `CHROM:POS:REF:ALT` alleles with
`bcftools isec -c none`. For variant sets A and B:

- intersection = |A intersect B|;
- union = |A| + |B| - |A intersect B|;
- Jaccard index = |A intersect B| / |A union B|; and
- directional sharing = |A intersect B| / |A| or |A intersect B| / |B|.

The directional fractions are included because a cache hit rate is directional:
it asks what fraction of an incoming sample is present in a cache. The symmetric
Jaccard index is therefore deliberately more conservative than the hit rate
expected from a cache built from a three-genome union.

## Results

| Comparison | Assembly | Records A | Records B | Intersection | Union | Jaccard | Shared from A | Shared from B |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| KPGP-00086 vs KPGP-00091 | GRCh38 | 4,789,204 | 4,690,148 | 2,993,035 | 6,486,317 | 46.14% | 62.50% | 63.82% |
| LP6005442-DNA_B03 vs LP6005442-DNA_B04 (SGDP) | GRCh38 | 4,754,466 | 4,823,693 | 2,982,654 | 6,595,505 | 45.22% | 62.73% | 61.83% |
| hu1D5A29 vs hu365511 (PGP) | GRCh37 | 3,641,537 | 5,005,430 | 2,294,775 | 6,352,192 | 36.13% | 63.02% | 45.85% |
| KPGP-00086 vs LP6005442-DNA_B03 (SGDP) | GRCh38 | 4,789,204 | 4,754,466 | 2,742,583 | 6,801,087 | 40.33% | 57.27% | 57.68% |

The KPGP and SGDP within-cohort pairs have similar Jaccard indices of about
45-46%, while the cross-cohort GRCh38 pair is lower at about 40%. The PGP pair
is more asymmetric because the two provider-derived callsets differ markedly
in size; nevertheless, 63.02% of the smaller callset is shared with the larger
one. These observations support benchmarking realistic, non-identical WGS
samples and measuring actual cache hit rates in the full campaign.

PGP is retained in its native GRCh37 coordinates, whereas KPGP and SGDP are
GRCh38. Raw-coordinate comparisons across those assemblies would not be
meaningful and were not performed. Such a comparison would require liftover,
renormalization, and a separate validation of reference-allele concordance.

## Reproduction

For each assembly-matched pair, the exact intersection count was obtained with:

```bash
bcftools isec -c none -n=2 -w1 A.vcf.gz B.vcf.gz \
  | bcftools query -f '%CHROM\n' \
  | wc -l
```

Input record counts were obtained with `bcftools index --nrecords`. The union
and ratios in the table were then calculated from those integer counts without
rounding; percentages are displayed to two decimal places.
