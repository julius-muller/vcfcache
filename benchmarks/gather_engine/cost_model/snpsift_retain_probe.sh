#!/bin/bash
# Does SnpSift keep the index of a chromosome after the input has moved past it? Smallest heap of the
# dbSNP group 1 stage (chr1-3) for 200 records of chr1 alone, of chr3 alone and of chr1 + chr3: if the
# stage frees a chromosome's index when the next one starts, chr1 + chr3 needs no more than the larger
# of the two (snpsift_heap_probe.sh: chr1-3 together need 8 GB). Needs the node-local sources.
# Usage: snpsift_retain_probe.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
OUT=/results/snpeff/mission/heap_probe; mkdir -p $OUT
SRC=/mnt/data/tmp/snpeff-mission/src; test -e $SRC/../SHARED || { echo "no node-local sources"; exit 1; }
W=$(mktemp -d /mnt/data/tmp/ser-XXXX); trap 'rm -rf $W' EXIT; cd $W
ERR=$W/err; CU=$SRC; T=$SRC; . /results/fastvep05/cost_model/snpeff_stack.sh
I=$KP/KPGP-00128.GRCh38.small_variants.vcf.gz
bcftools view -h $I > c1.vcf; bcftools view -H -r chr1 $I | head -200 >> c1.vcf
bcftools view -h $I > c3.vcf; bcftools view -H -r chr3 $I | head -200 >> c3.vcf
cat c1.vcf > c13.vcf; grep -v '^#' c3.vcf >> c13.vcf
echo "node $(hostname)"
for in in c1 c3 c13; do
  for h in 2g 3g 4g 5g 6g 8g; do
    if /usr/bin/time -v -o tv $JAVA -Xmx$h -jar $SE/SnpSift.jar annotate -noLog -c $SE/snpEff/snpEff.config -noId $SRC/dbsnp_g1.vcf - < $in.vcf > o.vcf 2> e \
       && [ "$(grep -vc '^#' o.vcf)" = "$(grep -vc '^#' $in.vcf)" ]; then
      printf "retain_%s\t%s\t%s\t%s\n" $in $h "$(awk -F': ' '/Elapsed/ {print $2}' tv)" "$(awk -F': ' '/Maximum resident/ {printf "%.2f", $2/1048576}' tv)" | tee -a $OUT/retain.tsv; break
    fi
    echo "$in: $h failed"
  done
done
echo DONE
