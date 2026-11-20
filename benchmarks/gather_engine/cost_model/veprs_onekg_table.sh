#!/bin/bash
# 1000 Genomes 30x sites as one genome-wide vep-rs dbNSFP-plugin table (allele frequency, overall
# and per super-population), for the larger vep-rs stack (wgs_mission_veprs.sh STACK2=1).
# Input: /results/sources/onekg/chr*.sites.vcf.gz (dl_onekg.sh). Writes /results/veprs05/mission/tables/onekg.tsv.gz
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
ls /results /results/veprs05 /results/sources/onekg > /dev/null 2>&1
T=/results/veprs05/mission/tables; K=/results/sources/onekg
W=$(mktemp -d /mnt/data/tmp/veprs-onekg-XXXX); trap 'rm -rf $W' EXIT; cd $W
CH="$(seq 1 22) X"   # the 30x panel has no chrY
for c in $CH; do test -s $K/chr$c.sites.vcf.gz || { echo "missing chr$c"; exit 1; }; done
for c in $CH; do
  ( bcftools norm -m- $K/chr$c.sites.vcf.gz -Ou 2>/dev/null | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AF_EAS\t%AF_AMR\t%AF_EUR\t%AF_AFR\t%AF_SAS\n' \
      | awk -F'\t' '$4!="*"' | bgzip -l1 > chr$c.gz || echo "FAILED chr$c" > chr$c.failed ) &
  while [ $(jobs -r | wc -l) -ge 14 ]; do wait -n; done
done
wait
ls *.failed 2>/dev/null && { cat *.failed; exit 1; }
{ printf "#chr\tpos(1-based)\tref\talt\tKG_AF\tKG_AF_EAS\tKG_AF_AMR\tKG_AF_EUR\tKG_AF_AFR\tKG_AF_SAS\n"; for c in $CH; do zcat chr$c.gz; done; } | bgzip -@14 > onekg.tsv.gz \
  && tabix -s1 -b2 -e2 -S1 onekg.tsv.gz || { echo "FAILED join"; exit 1; }
echo "onekg: $(zcat onekg.tsv.gz | tail -n +2 | wc -l) rows, $(du -h onekg.tsv.gz | cut -f1)"
cp onekg.tsv.gz onekg.tsv.gz.tbi $T/ && touch $T/../ONEKG && echo DONE
