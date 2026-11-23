#!/bin/bash
# Exome and panel inputs: which fastVEP source folder should a cached run use for the misses?
# The sample's whole genome restricted to the Twist core exome targets (pad 125) and to the
# ACMG SF v3.3 CDS regions (pad 20); arms as in wgs_mission.sh N10=1: uncached (65 kb sources),
# the 1% and cohort caches with misses at 65 kb and at 1 kb, and the rerun cache; 3 repeats.
# Needs the node-local state of wgs_mission.sh N10=1 (same STACK2 and sample on this node).
# Usage: [STACK2=1] assay_test.sh <deployment> [sample]
set -uo pipefail
DEP=$1; SAMPLE=${2:-KPGP-00128}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
RG=/results/fastvep05/cost_model/regions   # copies of the head node's /mnt/data/vcfcache_benchmarks/regions/GRCh38
ls /results/fastvep05 /results/fastvep05/cost_model $RG > /dev/null 2>&1
R=/mnt/data/tmp/mission-wgs${STACK2:+2}f; W=$R/$SAMPLE
OUT=/results/fastvep05/cost_model/assay${STACK2:+2}_$SAMPLE; mkdir -p $OUT
for f in $R/SOURCES10 $W/CACHES10 $W/in.bcf $RG/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed \
         $RG/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed; do test -e $f || { echo "missing $f"; exit 1; }; done
cd $W
bcftools view -R $RG/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed in.bcf -Ob -W -o wes.bcf
bcftools view -R $RG/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed in.bcf -Ob -W -o panel.bcf
echo "node $(hostname); exome $(bcftools index -n wes.bcf) variants, panel $(bcftools index -n panel.bcf)" | tee $OUT/info.txt
digest() { bcftools view -H $1 | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12; }
vc() {  # assay arm cache [--uncached]
  local C=$W/caches/$3/cache/$3
  rm -rf sta outa.bcf*
  local s=$(date +%s.%N)
  env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $C -i $1.bcf -o outa.bcf \
    -y $C/params.snapshot.yaml --statistics light --stats-dir sta -v ${4:-} > loga.txt 2>&1 || { echo "FAILED $1 $2"; tail -5 loga.txt; return; }
  local e=$(echo "$(date +%s.%N) - $s" | bc)
  local miss=$(grep -ho "missing_variants: [0-9]*" -r sta | head -1 | cut -d' ' -f2)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$(bcftools index -n $1.bcf)" "${miss:-NA}" "$e" "$(digest outa.bcf)" | tee -a $OUT/runs.tsv
}
printf "assay\tarm\tvariants\tmisses\tseconds\tdigest\n" > $OUT/runs.tsv
for a in wes panel; do
  vc $a warmup af001 --uncached > /dev/null
  for rep in 1 2 3; do
    vc $a uncached af001 --uncached
    for k in af001 af001_n10 cohort cohort_n10 rerun; do vc $a $k $k; done
  done
done
rm -f wes.bcf* panel.bcf* outa.bcf*
echo DONE
