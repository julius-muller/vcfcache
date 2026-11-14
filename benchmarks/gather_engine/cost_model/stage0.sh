#!/bin/bash
# Stage 0 anchors on one node: fastVEP 0.5.0 on a whole genome, cached and
# uncached; fastVEP alone on the full input and on the misses; vcfcache's own
# cost on inputs that are all cache hits; miss rates against two blueprints.
# Usage: stage0.sh <deployment> [sample]; CACHE, FLAGS and TAG in the
# environment select another cache, its fastVEP flags and the result folder.
set -uo pipefail
DEP=$1; SAMPLE=${2:-KPGP-00128}
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
OUT=/results/fastvep05/cost_model/${TAG:-stage0}; mkdir -p $OUT
W=/mnt/data/tmp/${TAG:-stage0}; rm -rf $W; mkdir -p $W; cd $W
SRC=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp/$SAMPLE.GRCh38.small_variants.vcf.gz
CD=${CACHE:-/results/fastvep05/caches/fastvep050_clinical_grch38_af001/cache/fastvep050_clinical_grch38_af001}
FLAGS=${FLAGS:---sa-dir /results/richness/sa/clinical_GRCh38}
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints
FV=/mnt/data/apps/fastvep/0.5.0
bcftools view $SRC -Ob1 -W --threads 8 -o wgs.bcf
N=$(bcftools index -n wgs.bcf)
sed "s#^temp_dir:.*#temp_dir: $W#" $CD/params.snapshot.yaml > p.yaml
echo "node $(hostname) records $N"
now() { date +%s.%N; }
vc() {  # label input [--uncached]
  rm -rf st out.bcf*
  local s=$(now)
  env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $CD -i $2 -o out.bcf \
    -y p.yaml --statistics light --stats-dir st -v ${3:-} > log.txt 2>&1 || { echo "FAILED $1"; tail -5 log.txt; return; }
  local e=$(echo "$(now) - $s" | bc)
  [ "$1" = wgs ] && [ -n "${3:-}" ] && [ ! -e dig.unc ] && bcftools view -H out.bcf | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12 > dig.unc
  [ "$1" = wgs ] && [ -z "${3:-}" ] && [ ! -e dig.c ] && bcftools view -H out.bcf | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12 > dig.c
  local miss=$(grep -ho "missing_variants: [0-9]*" -r st | head -1 | cut -d' ' -f2)
  printf "vcfcache\t%s\t%s\t%s\t%s\t%s\n" "$1" "${3:-cached}" "$(bcftools index -n $2)" "${miss:-NA}" "$e" | tee -a $OUT/runs.tsv
}
fv() {  # label input : fastVEP alone, same pipeline as the cache recipe
  local s=$(now)
  /usr/bin/time -f "%U %S %M" -o t.txt bash -c "bcftools view $2 -Ov | env RAYON_NUM_THREADS=8 $FV/bin/fastvep annotate --input - --output - --output-format vcf --transcript-cache $FV/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache --fasta $FV/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa --hgvs --no-progress $FLAGS 2> err.txt | bcftools view -Ob -o d.bcf -W --threads 8"
  local e=$(echo "$(now) - $s" | bc)
  printf "fastvep\t%s\talone\t%s\t%s\t%s\t%s\n" "$1" "$(bcftools index -n $2)" "$(cat t.txt | tr ' ' ',')" "$e" "$(bcftools index -n d.bcf)" | tee -a $OUT/runs.tsv
}
: > $OUT/runs.tsv
# 1. anchor
for rep in 1 2 3; do vc wgs wgs.bcf; vc wgs wgs.bcf --uncached; done
# misses of this cache, from the last cached run's work files
bcftools isec -n=1 -w1 -c none wgs.bcf $CD/vcfcache_annotated.bcf -Ob -o misses.bcf 2> isec.err; bcftools index -f misses.bcf
echo "digest cached $(cat dig.c) uncached $(cat dig.unc)" | tee -a $OUT/missrate.txt
# 2. fastVEP alone: full input, misses only
for rep in 1 2; do fv full wgs.bcf; fv misses misses.bcf; done
# 3. vcfcache alone: inputs made of cache hits only, at four sizes
bcftools isec -n=2 -w1 -c none wgs.bcf $CD/vcfcache_annotated.bcf -Ob -o hits_all.bcf 2> isec.err; bcftools index -f hits_all.bcf
for k in 1 2 4 8; do
  bcftools view -H hits_all.bcf | awk -v k=$k 'NR % k == 0' > body.txt
  { bcftools view -h hits_all.bcf; cat body.txt; } | bcftools view -Ob1 -W -o hits$k.bcf
  for rep in 1 2; do vc hits_1/$k hits$k.bcf; done
done
# 4. miss rates against both blueprints (exact biallelic match)
for b in af001 af0001; do
  f=$(ls $BP/gnomad_v4.1_GRCh38_joint_$b/blueprint/*.bcf | head -1)
  h=$(bcftools isec -n=2 -w1 -c none wgs.bcf $f 2>/dev/null | grep -vc '^#')
  echo "blueprint $b file $f records $(bcftools index -n $f) sample_hits $h of $N" | tee -a $OUT/missrate.txt
done
cd /; rm -rf $W; echo DONE
