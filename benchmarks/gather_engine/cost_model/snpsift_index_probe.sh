#!/bin/bash
# SnpSift annotate reads a source either through its tabix index (bgzipped VCF) or through its own
# index (.sidx, built on first use of a plain VCF). Which is faster? chr22 of KPGP-00128 after
# SnpEff, against chr22 of the three slowest sources (dbSNP, gnomAD, PhyloP), each mode run twice
# after the .sidx is built. Usage: snpsift_index_probe.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
OUT=/results/snpeff/mission/index_probe; mkdir -p $OUT
W=$(mktemp -d /mnt/data/tmp/sidx-XXXX); trap 'rm -rf $W' EXIT; cd $W
ERR=$W/err; CU=/results/vep116/mission/custom; T=.; . /results/fastvep05/cost_model/snpeff_stack.sh
bcftools view -r chr22 $KP/KPGP-00128.GRCh38.small_variants.vcf.gz -Ov | eval "$SNPEFF" > se.vcf
SIFT="$JAVA -Xmx2g -jar $SE/SnpSift.jar annotate -noLog -c $SE/snpEff/snpEff.config -noId"
t() { local s=$(date +%s.%N); bash -c "$1" 2>> $ERR; echo "$(date +%s.%N) - $s" | bc; }
printf "source\tmode\trun\tseconds\tmd5\n" > $OUT/probe.tsv
for s in dbsnp gnomad phylop; do
  bcftools view -r chr22 $CU/$s.vcf.gz -Oz --write-index=tbi -o $s.vcf.gz
  bcftools view $s.vcf.gz -Ov -o $s.vcf
  e=$(t "$SIFT $s.vcf se.vcf > /dev/null"); printf "%s\tsidx\tbuild\t%s\t-\n" $s $e | tee -a $OUT/probe.tsv
  for r in 1 2; do
    for m in tabix sidx; do
      f=$([ $m = tabix ] && echo $s.vcf.gz || echo $s.vcf)
      e=$(t "$SIFT $f se.vcf > out.$m.vcf")
      printf "%s\t%s\t%s\t%s\t%s\n" $s $m $r $e $(grep -v '^##SnpSift' out.$m.vcf | md5sum | cut -c1-12) | tee -a $OUT/probe.tsv
    done
  done
done
tail -n 5 $ERR
echo DONE
