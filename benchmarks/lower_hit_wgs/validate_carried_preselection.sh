#!/usr/bin/env bash
set -euo pipefail
screen_dir=/mnt/data/vcfcache_benchmarks/lower_hit_wgs_20261014
source_vcf=/mnt/data/vcfcache_benchmarks/external_wgs/sources/sgdp/LP6005592-DNA_C05/LP6005592-DNA_C05.hs38DH.BQSR.autosome.g.vcf.gz
reference=/mnt/data/resources/reference/ucsc/hg38.fa.gz
rename_map=/mnt/data/vcfcache_benchmarks/external_wgs/work/normalization/numeric_to_chr.tsv
predicate='(GT="alt" && ALT!="<NON_REF>" && ALT!="*" && TYPE="snp") || (GT="alt" && TYPE="indel" && strlen(REF)<=50 && strlen(ALT)<=50)'
for mode in original preselected; do
  options=()
  if [[ "$mode" == preselected ]]; then options=(-i 'GT="alt"'); fi
  bcftools view -r chr1:1000000-3000000 -f PASS,. "${options[@]}" -Ou "$source_vcf" |
    bcftools annotate --rename-chrs "$rename_map" -Ou |
    bcftools norm -f "$reference" -m -any -Ou |
    bcftools view -i "$predicate" -Ou |
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' |
    LC_ALL=C sort -u > "$screen_dir/results/preselection_${mode}.keys.tsv"
done
cmp "$screen_dir/results/preselection_original.keys.tsv" "$screen_dir/results/preselection_preselected.keys.tsv"
test -s "$screen_dir/results/preselection_original.keys.tsv"
wc -l "$screen_dir/results/preselection_original.keys.tsv"
sha256sum "$screen_dir/results/preselection_original.keys.tsv" "$screen_dir/results/preselection_preselected.keys.tsv"
