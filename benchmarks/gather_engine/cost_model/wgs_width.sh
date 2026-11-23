#!/bin/bash
# Does fastVEP get slower on dense input when its sources are stored in 1 kb blocks?
# If not, a user needs one source set (1 kb) for building the cache, for the misses and
# for direct annotation; if so, two (65 kb to build, 1 kb for the misses).
# On one whole genome, alternating, 3 repeats: the uncached pipeline (vcfcache --uncached)
# at 65 kb and at 1 kb; then once each the cache build (the 1% blueprint, 18.9 M records,
# 16 threads) at 65 kb and at 1 kb, with digests.
# Needs the node-local state of wgs_mission.sh N10=1 (same STACK2 and sample on this node).
# Usage: [STACK2=1] wgs_width.sh <deployment> [sample]
set -uo pipefail
DEP=$1; SAMPLE=${2:-KPGP-00128}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
ls /results/fastvep05 /results/fastvep05/cost_model > /dev/null 2>&1
R=/mnt/data/tmp/mission-wgs${STACK2:+2}f; W=$R/$SAMPLE
OUT=/results/fastvep05/cost_model/width${STACK2:+2}_wgs_$SAMPLE; mkdir -p $OUT
for f in $R/SOURCES10 $W/CACHES10 $R/sa16 $R/sa10 $W/in.bcf $W/bp_af001.bcf; do test -e $f || { echo "missing $f"; exit 1; }; done
cd $W; log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname), fastvep $($FV --version); sources $(du -sh $R/sa16 | cut -f1) at 65 kb, $(du -sh $R/sa10 | cut -f1) at 1 kb"
N=$(bcftools index -n in.bcf)
digest() { bcftools view -H $1 | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12; }
vc() {  # arm cache [--uncached]
  local C=$W/caches/$2/cache/$2
  rm -rf stw outw.bcf*
  local s=$(date +%s.%N)
  env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $C -i in.bcf -o outw.bcf \
    -y $C/params.snapshot.yaml --statistics light --stats-dir stw -v ${3:-} > logw.txt 2>&1 || { echo "FAILED $1"; tail -5 logw.txt; return; }
  local e=$(echo "$(date +%s.%N) - $s" | bc)
  [ -e digw.$1 ] || digest outw.bcf > digw.$1
  printf "%s\t%s\t%s\t%s\n" "$1" "$N" "$e" "$(cat digw.$1)" | tee -a $OUT/runs.tsv
}
: > $OUT/runs.tsv
for rep in 1 2 3; do vc uncached_65kb af001 --uncached; vc uncached_1kb af001_n10 --uncached; done
# the cache build: annotate the 1% blueprint at each width (as wgs_mission.sh's blueprint())
for b in 16 10; do
  s=$(date +%s.%N)
  bcftools view bp_af001.bcf -Ov | env RAYON_NUM_THREADS=16 $FV annotate --input - --output - --output-format vcf \
    --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $R/sa$b 2> annow.err | bcftools view -Ob -o bpw_$b.bcf -W --threads 8 \
    || { echo "FAILED build $b"; tail -3 annow.err; exit 1; }
  e=$(echo "$(date +%s.%N) - $s" | bc)
  printf "%s\t%s\t%s\t%s\n" build_$([ $b = 16 ] && echo 65kb || echo 1kb) $(bcftools index -n bpw_$b.bcf) $e "$(digest bpw_$b.bcf)" | tee -a $OUT/runs.tsv
done
rm -f bpw_*.bcf*
echo DONE
