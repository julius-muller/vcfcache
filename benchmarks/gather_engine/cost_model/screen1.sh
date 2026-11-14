#!/bin/bash
# Stage 1 screen with sources already built: fastVEP 0.5.0 alone on chr1 of
# one genome, at five input sizes and on the cache misses, for eight stacks.
# The slope over input size is the per-variant cost b; the intercept is
# start-up plus the cost of opening source blocks.
# Usage: screen1.sh [sample]
set -uo pipefail
SAMPLE=${1:-KPGP-00128}
OUT=/results/fastvep05/cost_model/screen1; mkdir -p $OUT
W=/mnt/data/tmp/screen1; rm -rf $W; mkdir -p $W; cd $W
SRC=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp/$SAMPLE.GRCh38.small_variants.vcf.gz
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
FV=/mnt/data/apps/fastvep/0.5.0
SA=/results/fastvep05/acmg/sa
BASE="--transcript-cache $FV/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache --fasta $FV/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa --no-progress"
bcftools view $SRC -r chr1 -Ob1 -W --threads 8 -o in1.bcf
for k in 2 4 8 16; do
  bcftools view -H in1.bcf | awk -v k=$k 'NR % k == 0' > body.txt
  { bcftools view -h in1.bcf; cat body.txt; } | bcftools view -Ob1 -W -o in$k.bcf
done
bcftools isec -n=1 -w1 -c none in1.bcf $BP -Ob -o inmiss.bcf 2>/dev/null; bcftools index -f inmiss.bcf
declare -A STACK=(
  [base]=""
  [hgvs]="--hgvs"
  [everything]="--everything"
  [clinical]="--hgvs --sa-dir /results/richness/sa/clinical_GRCh38"
  [enriched]="--hgvs --sa-dir /results/richness/sa/enriched"
  [acmg]="--hgvs --acmg --sa-dir $SA/stack"
  [dense]="--hgvs --sa-dir $SA/dense"
  [dense_acmg]="--hgvs --acmg --sa-dir $SA/dense"
)
: > $OUT/runs.tsv
echo "node $(hostname) chr1 records $(bcftools index -n in1.bcf) misses $(bcftools index -n inmiss.bcf)"
for rep in 1 2; do
  for s in base hgvs everything clinical enriched acmg dense dense_acmg; do
    for i in 1 2 4 8 16 miss; do
      s0=$(date +%s.%N)
      /usr/bin/time -f "%U\t%S\t%M" -o t.txt bash -c "bcftools view in$i.bcf -Ov | env RAYON_NUM_THREADS=8 $FV/bin/fastvep annotate --input - --output - --output-format vcf $BASE ${STACK[$s]} 2> err.txt | bcftools view -Ob -o d.bcf -W --threads 8" || tail -2 err.txt
      e=$(echo "$(date +%s.%N) - $s0" | bc)
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $s $i $rep "$(bcftools index -n in$i.bcf)" "$e" "$(tail -1 t.txt)" "$(bcftools index -n d.bcf)" "$(stat -c %s d.bcf)" >> $OUT/runs.tsv
    done
  done
done
cd /; rm -rf $W; echo DONE
