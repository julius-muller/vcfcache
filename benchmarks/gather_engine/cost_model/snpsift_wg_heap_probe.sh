#!/bin/bash
# Smallest Java heap of SnpSift stages on a whole genome (all 4.6 M records of KPGP-00128), the input
# of the uncached arm. SnpSift keeps the database records it has read, so the heap grows with how
# much of a source's chromosomes the input touches: part new-0 (chr1-3 of the new sites) ran the
# dbSNP group 1 stage out of a 2 GB heap twice (jobs 3710, 3711), and 200 records per contig
# (snpsift_heap_probe.sh) touch almost nothing. Two stages at a time; for each, heaps from small to
# large, the first that completes. PERCHROM=1 also runs dbSNP group 1 on chr1, chr2 and chr3 alone:
# if the whole genome needs no more than the largest of them, the stage frees a chromosome's records
# when the next one starts.
# Usage: [PERCHROM=1] snpsift_wg_heap_probe.sh <stage> [<stage>]
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
OUT=/results/snpeff/mission/heap_probe; mkdir -p $OUT
SRC=/mnt/data/tmp/snpeff-mission/src; test -e $SRC/../SHARED || { echo "no node-local sources"; exit 1; }
W=$(mktemp -d /mnt/data/tmp/sew-XXXX); trap 'rm -rf $W' EXIT; cd $W
ERR=$W/err; CU=$SRC; T=$SRC; . /results/fastvep05/cost_model/snpeff_stack.sh
I=$KP/KPGP-00128.GRCh38.small_variants.vcf.gz
bcftools view $I -Ov -o wg.vcf; N=$(grep -vc '^#' wg.vcf)
echo "node $(hostname), $N records"
try() {  # stage input heaps...
  local s=$1 in=$2 h n; shift 2; n=$(grep -vc '^#' $in)
  for h in "$@"; do
    if /usr/bin/time -v -o $s.$in.tv $JAVA -Xmx$h -jar $SE/SnpSift.jar annotate -noLog -c $SE/snpEff/snpEff.config -noId \
         $SRC/$s.vcf - < $in > $s.$in.out 2> $s.$in.err && [ "$(grep -vc '^#' $s.$in.out)" = "$n" ]; then
      printf "wg\t%s\t%s\t%s\t%s\t%s\n" $s $in $h "$(awk -F': ' '/Elapsed/ {print $2}' $s.$in.tv)" \
        "$(awk -F': ' '/Maximum resident/ {printf "%.2f", $2/1048576}' $s.$in.tv)" | tee -a $OUT/wg.tsv
      rm -f $s.$in.out; return
    fi
    echo "$s $in: $h failed ($(grep -m1 -o 'OutOfMemoryError' $s.$in.err))"; rm -f $s.$in.out
  done
  printf "wg\t%s\t%s\tnone\n" $s $in | tee -a $OUT/wg.tsv
}
H="2g 3g 4g 6g 8g 12g"
for s in "$@"; do try $s wg.vcf $H & done; wait
if [ "${PERCHROM:-0}" = 1 ]; then
  for c in chr1 chr2 chr3; do bcftools view -r $c $I -Ov -o $c.vcf; done
  try dbsnp_g1 chr1.vcf $H & try dbsnp_g1 chr2.vcf $H & wait; try dbsnp_g1 chr3.vcf $H
fi
echo DONE
