#!/bin/bash
# Does the serial collector of the SnpSift stages (snpeff_stack.sh default, needed so the whole-genome
# pipeline fits a node) change output or time against the JVM's default collector? Both on chr1-4 of
# KPGP-00128 (fits either way), same node, alternating A B A B. Exits 1 if the outputs differ or the
# serial runs take more than 3% longer, so jobs that depend on it (afterok) do not start. Job 3721:
# serial 2.6% slower (3878 vs 3780 s), outputs identical; job 3722, G1 returning memory
# (ALT=-XX:G1PeriodicGCInterval=10000 -XX:MinHeapFreeRatio=10 -XX:MaxHeapFreeRatio=30): 1935 and
# 1978 s against 1917 and 2158 s, outputs identical. ALT (JVM options) and ALTN (its name) test
# another collector against the default one.
# Usage: [ALT=<jvm options> ALTN=<name>] snpeff_gc_check.sh [regions]
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
REG=${1:-chr1,chr2,chr3,chr4}
OUT=/results/snpeff/mission/pipe_probe; mkdir -p $OUT
SRC=/mnt/data/tmp/snpeff-mission/src; test -e $SRC/../SHARED || { echo "no node-local sources"; exit 1; }
W=$(mktemp -d /mnt/data/tmp/segc-XXXX); trap 'rm -rf $W' EXIT; cd $W
ERR=$W/err; CU=$SRC; T=$SRC
ALT=${ALT:--XX:+UseSerialGC -XX:MinHeapFreeRatio=10 -XX:MaxHeapFreeRatio=30}; ALTN=${ALTN:-serial}
echo "node $(hostname), regions $REG, $ALTN: $ALT"
for arm in g1 $ALTN g1 $ALTN; do
  if [ $arm = g1 ]; then SIFT_JOPT=; else SIFT_JOPT=$ALT; fi
  . /results/fastvep05/cost_model/snpeff_stack.sh
  s=$(date +%s)
  bcftools view -r $REG $KP/KPGP-00128.GRCh38.small_variants.vcf.gz -Ov | eval "$SPIPE" | bcftools view -Ob1 -W --threads 8 -o out.bcf
  st="${PIPESTATUS[*]}"; e=$(( $(date +%s) - s ))
  n=$(bcftools index -n out.bcf); md=$(bcftools view -H out.bcf | md5sum | cut -c1-12)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $(hostname) $REG $arm "$st" $e $n $md | tee -a $OUT/gc_check.tsv
  echo "$arm $e $md $st" >> runs
done
awk 'NR==1 {m=$3; ok=1} $3!=m {d=1} ($4 " " $5 " " $6) != "0 0 0" {ok=0} {t[$1]+=$2} END {
  printf "g1 %d s, '$ALTN' %d s (%.3fx), outputs %s, status %s\n", t["g1"], t["'$ALTN'"], t["'$ALTN'"]/t["g1"], d ? "DIFFER" : "identical", ok ? "ok" : "FAILED"
  exit (d || !ok || t["'$ALTN'"] > 1.03*t["g1"]) }' runs
