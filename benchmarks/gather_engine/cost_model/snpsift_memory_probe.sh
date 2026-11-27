#!/bin/bash
# How much Java heap does SnpSift need for the whole-genome plain sources? A 2 GB heap runs out of
# memory while it builds the .sidx of the large sources. Builds each index with a 16 GB heap, one at a
# time, then annotates chr22 of KPGP-00128 (after SnpEff) with each source; peak RSS and time of
# every step. Runs on a node whose /mnt/data/tmp/snpeff-mission/src holds the decompressed sources.
# HEAP (default 16g) sets the heap, ONLY limits the run to the named sources (e.g. ONLY=dbsnp).
# Usage: [HEAP=26g] [ONLY="dbsnp ..."] snpsift_memory_probe.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
OUT=/results/snpeff/mission/memory_probe; mkdir -p $OUT
R=/mnt/data/tmp/snpeff-mission; SRC=$R/src
W=$(mktemp -d /mnt/data/tmp/sem-XXXX); trap 'rm -rf $W' EXIT; cd $W
ERR=$W/err; CU=$SRC; T=$SRC; . /results/fastvep05/cost_model/snpeff_stack.sh
ls -la $SRC > $OUT/src.txt
SIFT="$JAVA -Xmx${HEAP:-16g} -jar $SE/SnpSift.jar annotate -noLog -c $SE/snpEff/snpEff.config -noId"
bcftools view -h $SRC/clinvar.vcf > one.vcf; bcftools view -H $SRC/clinvar.vcf | head -1 >> one.vcf
bcftools view -r chr22 $KP/KPGP-00128.GRCh38.small_variants.vcf.gz -Ov | eval "$SNPEFF" > se.vcf
m() {  # step command: wall seconds and peak RSS (GB)
  /usr/bin/time -v -o tv bash -c "$2" 2>> $ERR || { echo "FAILED $1"; tail -n 5 $ERR; }
  printf "%s\t%s\t%s\n" "$1" "$(awk -F': ' '/Elapsed/ {print $2}' tv)" \
    "$(awk -F': ' '/Maximum resident/ {printf "%.2f", $2/1048576}' tv)" | tee -a $OUT/probe.tsv
}
[ -e $OUT/probe.tsv ] || printf "step\twall\tpeak_rss_gb\n" > $OUT/probe.tsv
for f in $SRC/*.vcf; do
  s=$(basename $f .vcf); [ -z "${ONLY:-}" ] || [[ " $ONLY " == *" $s "* ]] || continue
  [ -e $f.sidx ] || m build_${s}_${HEAP:-16g} "$SIFT $f one.vcf > /dev/null"
  m annotate_${s}_${HEAP:-16g} "$SIFT $f se.vcf > /dev/null"
done
ls -la $SRC >> $OUT/src.txt
echo DONE
