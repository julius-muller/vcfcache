#!/bin/bash
# Smallest Java heap per stage of the SnpEff pipeline on whole-genome input. SnpSift loads the index of
# a source one chromosome at a time as the input reaches it, so a chr22 input (snpsift_memory_probe.sh)
# understates the heap: whole-genome parts ran out of memory at 3 GB (jobs 3664-3669). Input: the first
# 200 records of every contig of KPGP-00128. For each stage, tries heaps from small to large and keeps
# the first that completes; peak RSS and wall time of that run. Needs the node-local sources.
# Usage: snpsift_heap_probe.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
OUT=/results/snpeff/mission/heap_probe; mkdir -p $OUT
SRC=/mnt/data/tmp/snpeff-mission/src; test -e $SRC/../SHARED || { echo "no node-local sources"; exit 1; }
W=$(mktemp -d /mnt/data/tmp/seh-XXXX); trap 'rm -rf $W' EXIT; cd $W
ERR=$W/err; CU=$SRC; T=$SRC; . /results/fastvep05/cost_model/snpeff_stack.sh
I=$KP/KPGP-00128.GRCh38.small_variants.vcf.gz
bcftools view -h $I > in.vcf
for c in $(bcftools index -s $I | cut -f1); do bcftools view -H -r $c $I | head -200 >> in.vcf; done
echo "node $(hostname), input $(grep -vc '^#' in.vcf) records on $(bcftools index -s $I | wc -l) contigs"
printf "stage\theap\twall\tpeak_rss_gb\n" > $OUT/probe.tsv
try() {  # stage heaps... command-with-HEAP
  local s=$1 cmd=$2 h; shift 2
  for h in "$@"; do
    if /usr/bin/time -v -o tv bash -c "${cmd//HEAP/$h} < in.vcf > o.vcf 2> e" && [ "$(grep -vc '^#' o.vcf)" = "$(grep -vc '^#' in.vcf)" ]; then
      printf "%s\t%s\t%s\t%s\n" $s $h "$(awk -F': ' '/Elapsed/ {print $2}' tv)" \
        "$(awk -F': ' '/Maximum resident/ {printf "%.2f", $2/1048576}' tv)" | tee -a $OUT/probe.tsv; return
    fi
    echo "$s: $h failed ($(grep -m1 -o 'OutOfMemoryError[^ ]*\|Exception[^ ]*' e))"
  done
  printf "%s\tnone\t\t\n" $s | tee -a $OUT/probe.tsv
}
try snpeff "$JAVA -XmxHEAP -jar $SE/snpEff/snpEff.jar -noLog -noStats -nodownload -c $SE/snpEff/snpEff.config -dataDir $SE/snpEff/data GRCh38.115 -" 3g 4g 5g 6g 8g
for f in clinvar gnomad dbsnp_g1 dbsnp_g2 dbsnp_g3 dbsnp_g4 dbsnp_g5 spliceai_snv spliceai_indel phylop trackA trackB; do
  try $f "$JAVA -XmxHEAP -jar $SE/SnpSift.jar annotate -noLog -c $SE/snpEff/snpEff.config -noId $SRC/$f.vcf -" 1g 2g 3g 4g 5g 6g 8g 10g 12g 16g
done
echo DONE
