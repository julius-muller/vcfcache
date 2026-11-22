#!/bin/bash
# Fixed per-run cost of fastVEP's sources at 65 kb vs 1 kb blocks: annotate random genome-wide
# subsets (0, 1k, 10k, 100k variants) of one sample with the node-local sources of
# wgs_mission.sh N10=1 (sa16, sa10), 3 repeats each.
# Usage: [STACK2=1] fixed_cost.sh [sample]
set -uo pipefail
SAMPLE=${1:-KPGP-00128}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
FVD=/mnt/data/apps/fastvep/0.5.0
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
ls /results/fastvep05 /results/fastvep05/cost_model > /dev/null 2>&1
R=/mnt/data/tmp/mission-wgs${STACK2:+2}f
OUT=/results/fastvep05/cost_model/fixed${STACK2:+2}_$SAMPLE; mkdir -p $OUT
for f in $R/sa16 $R/sa10; do test -d $f || { echo "missing $f"; exit 1; }; done
W=$(mktemp -d /mnt/data/tmp/fixed-XXXX); trap 'rm -rf $W' EXIT; cd $W
N=$(bcftools index -n $KP/$SAMPLE.GRCh38.small_variants.vcf.gz)
for n in 0 1000 10000 100000; do
  bcftools view $KP/$SAMPLE.GRCh38.small_variants.vcf.gz | awk -v p=$(echo "$n / $N" | bc -l) 'BEGIN{srand(7)} /^#/{print; next} rand()<p' > in$n.vcf
done
fv() { env RAYON_NUM_THREADS=8 $FVD/bin/fastvep annotate --input $1 --output - --output-format vcf \
  --transcript-cache $FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache --fasta $FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
  --hgvs --no-progress --sa-dir $2 2>/dev/null | grep -vc "^#"; }
echo "node $(hostname); sources $(du -sh $R/sa16 | cut -f1) at 65 kb, $(du -sh $R/sa10 | cut -f1) at 1 kb" | tee $OUT/info.txt
fv in1000.vcf $R/sa16 > /dev/null; fv in1000.vcf $R/sa10 > /dev/null   # warm
printf "variants\tblocks\tseconds\trecords\n" > $OUT/runs.tsv
for rep in 1 2 3; do for n in 0 1000 10000 100000; do for b in 16 10; do
  s=$(date +%s.%N); r=$(fv in$n.vcf $R/sa$b); printf "%s\t%s\t%.2f\t%s\n" $n $b $(echo "$(date +%s.%N) - $s" | bc) $r | tee -a $OUT/runs.tsv
done; done; done
echo DONE
