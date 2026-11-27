#!/bin/bash
# The uncached arm's whole pipeline (snpeff_stack.sh) on all of KPGP-00128 with given heaps and
# SnpSift JVM options, the resident memory of every stage sampled every 10 s. SnpSift keeps the
# database records it has read: at the smallest heaps that complete a whole genome
# (snpsift_wg_heap_probe.sh) the 13 processes would hold about 31 GB, more than a node has.
# Env: SNPEFF_HEAP, GNOMAD_HEAP, DBSNP_HEAP, SIFT_JOPT (see snpeff_stack.sh); TAG names the run.
# Usage: TAG=<name> [...] snpeff_pipe_probe.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
OUT=/results/snpeff/mission/pipe_probe; mkdir -p $OUT
SRC=/mnt/data/tmp/snpeff-mission/src; test -e $SRC/../SHARED || { echo "no node-local sources"; exit 1; }
W=$(mktemp -d /mnt/data/tmp/sep-XXXX); trap 'kill $mon 2>/dev/null; rm -rf $W' EXIT; cd $W
ERR=$W/err; CU=$SRC; T=$SRC; . /results/fastvep05/cost_model/snpeff_stack.sh
echo "node $(hostname), TAG=$TAG SNPEFF_HEAP=${SNPEFF_HEAP:-} GNOMAD_HEAP=${GNOMAD_HEAP:-} DBSNP_HEAP=${DBSNP_HEAP:-} SIFT_JOPT=${SIFT_JOPT:-}"
while sleep 10; do
  ps -u $(id -u) -o rss=,args= | awk -v t=$(date +%s) '/jre.bin.java/ {s="snpeff"; for (i=1;i<=NF;i++) if ($i ~ /src\/.*\.vcf$/) {s=$i; sub(/.*\//,"",s); sub(/\.vcf$/,"",s)}; printf "%s\t%s\t%d\n", t, s, $1}'
  awk -v t=$(date +%s) '/MemAvailable/ {printf "%s\tavail\t%d\n", t, $2}' /proc/meminfo
done > rss.tsv & mon=$!
s=$(date +%s)
bcftools view $KP/KPGP-00128.GRCh38.small_variants.vcf.gz -Ov | eval "$SPIPE" | bcftools view -Ob1 -W --threads 8 -o out.bcf
st="${PIPESTATUS[*]}"; e=$(( $(date +%s) - s )); kill $mon
n=$(bcftools index -n out.bcf 2>/dev/null); md=$(bcftools view -H out.bcf 2>/dev/null | md5sum | cut -c1-12)
# peak total resident memory of the pipeline, minimum MemAvailable, peak per stage
tot=$(awk '$2!="avail" {a[$1]+=$3} END {for (t in a) if (a[t]>m) m=a[t]; printf "%.1f", m/1048576}' rss.tsv)
av=$(awk '$2=="avail" {if (m=="" || $3<m) m=$3} END {printf "%.1f", m/1048576}' rss.tsv)
per=$(awk '$2!="avail" {if ($3>m[$2]) m[$2]=$3} END {for (k in m) printf "%s=%.1f ", k, m[k]/1048576}' rss.tsv)
echo "exit: $st"; grep -m3 -i 'exception\|error' $ERR
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $TAG "$st" $e "$n" $md $tot $av "$per" | tee -a $OUT/probe.tsv
cp rss.tsv $OUT/rss_$TAG.tsv; echo DONE
