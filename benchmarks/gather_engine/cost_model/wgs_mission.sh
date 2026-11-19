#!/bin/bash
# The mission pipeline on one whole genome: fastVEP 0.5.0 with eleven known
# sources (fifteen with STACK2=1) and two mock custom tracks, uncached against vcfcache with
#   - the gnomAD AF >= 0.1% and AF >= 0.01% blueprints,
#   - the misses annotated from the same sources stored in narrow blocks
#     (the uncached baseline always uses the sources as fastVEP advises),
#   - a cohort cache (blueprint plus the other KPGP samples), and
#   - a rerun cache that also holds the sample itself (no misses).
# Usage: wgs_mission.sh <deployment> [sample]
set -uo pipefail
DEP=$1; SAMPLE=${2:-KPGP-00128}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
A=/results/fastvep05/acmg; S=$A/sources
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
BPD=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints
OLD=/results/fastvep05/caches/fastvep050_clinical_grch38_af001
MS=/results/fastvep05/mission
ls /results/fastvep05 $MS > /dev/null 2>&1   # clears a stale "no such directory" on some nodes
# STACK2=1 adds four sources of wgs_mission_stack2.sh (1000 Genomes, GERP, dbNSFP, ENCODE cCREs).
# fastVEP makes every custom_vcf file a source of its own (a second file with the same --name
# becomes FV_<NAME>_2, ...), so a custom source must be ONE genome-wide file: the mock tracks
# come from build_tracks1.sh, and gnomAD genomes (custom, built in 20 Mb pieces) is left out.
MS2=/results/fastvep05/mission2; T1=$MS/tracks1
OUT=/results/fastvep05/cost_model/mission${STACK2:+2}f_wgs_$SAMPLE; mkdir -p $OUT
rm -rf /mnt/data/tmp/mission-wgs /mnt/data/tmp/mission-wgs2   # node-local copies of the earlier, flawed setup
R=/mnt/data/tmp/mission-wgs${STACK2:+2}f; mkdir -p $R/tmp; cd $R
M16=$R/sa16; M12=$R/sa12; mkdir -p $M16 $M12
log() { echo "[$(date -Is)] $*"; }
b() { local s=$(date +%s); $FV sa-build --no-progress --assembly GRCh38 "$@" > build.$$.log 2>&1 || { cat build.$$.log; return 1; }; log "sa-build $1 $2 ${*: -2} $(( $(date +%s) - s ))s"; }
echo "node $(hostname), fastvep $($FV --version)"

# ---- sources: built per chromosome by wgs_mission_stack.sh; local copies ----
if [ ! -e $R/SOURCES ]; then
  for c in $(seq 1 22) X Y; do test -e $MS/done/chr$c || { echo "sources of chr$c missing"; exit 1; }; done
  for b in 16 12; do
    for f in $MS/sa$b/*; do case $(basename $f) in tracka_chr*|trackb_chr*) ;; *) cp $f $R/sa$b/ || exit 1;; esac; done
    cp $T1/sa$b/track[ab].* $R/sa$b/ || exit 1
  done
  if [ -n "${STACK2:-}" ]; then
    for c in $(seq 1 22) X Y; do test -e $MS2/done/chr$c || { echo "further sources of chr$c missing"; exit 1; }; done
    for b in 16 12; do
      for f in $MS2/sa$b/*; do case $(basename $f) in gnomad_genomes_*) ;; *) cp $f $R/sa$b/ || exit 1;; esac; done
    done
  fi
  [ $(ls $M16 | wc -l) = $(ls $M12 | wc -l) ] || { echo "source sets differ"; exit 1; }
  ls $M16 | sed 's/_chr.*//' | sort | uniq -c | tee $OUT/sources.txt; du -sh $M16 $M12 | tee -a $OUT/sources.txt
  touch $R/SOURCES
fi

# sources and the annotated blueprints are shared by all samples on a node; the rest is per sample
W=$R/$SAMPLE; mkdir -p $W/tmp
# ---- inputs and caches -----------------------------------------------------
anno() {  # in.bcf out.bcf sa-dir threads
  bcftools view $1 -Ov | env RAYON_NUM_THREADS=$4 $FV annotate --input - --output - --output-format vcf \
    --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $3 2> anno.err | bcftools view -Ob -o $2 -W --threads 8
}
blueprint() {  # af001 | af0001: annotate the gnomAD blueprint once per node, link it into the sample directory
  local t=$1 s=$(date +%s)
  ln -sf $BPD/gnomad_v4.1_GRCh38_joint_$t/blueprint/vcfcache.bcf bp_$t.bcf; ln -sf $BPD/gnomad_v4.1_GRCh38_joint_$t/blueprint/vcfcache.bcf.csi bp_$t.bcf.csi
  if [ ! -e $R/c_$t.bcf.csi ]; then
    anno bp_$t.bcf $R/c_$t.part.bcf $M16 16 && mv $R/c_$t.part.bcf $R/c_$t.bcf && mv $R/c_$t.part.bcf.csi $R/c_$t.bcf.csi
    log "cache $t: $(bcftools index -n $R/c_$t.bcf) records, $(( $(date +%s) - s ))s, $(du -h $R/c_$t.bcf | cut -f1)"
    bcftools view -h $R/c_$t.bcf | grep -o '^##INFO=<ID=FV_[A-Z0-9_]*' | cut -d= -f3 | tee $OUT/fields.txt | grep -q '_[0-9][0-9]*$' && { echo "a source appears more than once"; exit 1; }
  fi
  ln -f $R/c_$t.bcf c_$t.bcf; ln -f $R/c_$t.bcf.csi c_$t.bcf.csi
}
mkcache() {  # name annotated.bcf sa-dir-for-misses
  local ROOT=$W/caches/$1 C=$W/caches/$1/cache/$1
  rm -rf $ROOT; mkdir -p $C $ROOT/blueprint; cp -r $OLD/workflow $ROOT/
  sed -e "s#--sa-dir /results/richness/sa/clinical_GRCh38#--sa-dir $3#" $OLD/cache/$(basename $OLD)/annotation.yaml > $C/annotation.yaml
  grep -q -- "--sa-dir $3" $C/annotation.yaml || { echo "recipe not rewritten"; exit 1; }
  sed "s#^temp_dir:.*#temp_dir: $W/tmp#" $OLD/cache/$(basename $OLD)/params.snapshot.yaml > $C/params.snapshot.yaml
  ln -f $2 $C/vcfcache_annotated.bcf; ln -f $2.csi $C/vcfcache_annotated.bcf.csi
}
cd $W
if [ ! -e $W/CACHES ]; then
  bcftools view $KP/$SAMPLE.GRCh38.small_variants.vcf.gz -Ob1 -W --threads 8 -o in.bcf
  blueprint af001
  # cohort: biallelic sites of the other samples, and of this sample, not in the 0.1% blueprint
  sites() { bcftools view -G $1 -Ou | bcftools norm -m- -Ou 2>/dev/null | bcftools annotate -x INFO,ID,QUAL,FILTER -Ob -W -o $2; }
  n=0; for f in $KP/*.small_variants.vcf.gz; do
    case $f in *$SAMPLE*) sites $f self.bcf;; *) n=$((n+1)); sites $f o$n.bcf;; esac
  done
  bcftools concat -a -d exact o*.bcf -Ob -W -o others.bcf --threads 8 2>/dev/null
  bcftools isec -n=1 -w1 -c none others.bcf bp_af001.bcf -Ob -o new_cohort.bcf 2>/dev/null; bcftools index -f new_cohort.bcf
  bcftools concat -a -d exact others.bcf self.bcf -Ob -W -o withself.bcf --threads 8 2>/dev/null
  bcftools isec -n=1 -w1 -c none withself.bcf bp_af001.bcf -Ob -o new_rerun.bcf 2>/dev/null; bcftools index -f new_rerun.bcf
  log "cohort: $n other samples, $(bcftools index -n new_cohort.bcf) sites beyond the blueprint; with the sample itself $(bcftools index -n new_rerun.bcf)"
  for k in cohort rerun; do
    anno new_$k.bcf a_$k.bcf $M16 16
    bcftools concat -a c_af001.bcf a_$k.bcf -Ob -W -o c_$k.bcf --threads 8 2>/dev/null
  done
  for k in af001 cohort rerun; do mkcache $k $W/c_$k.bcf $M16; mkcache ${k}_narrow $W/c_$k.bcf $M12; done
  rm -f o*.bcf o*.bcf.csi; touch $W/CACHES
fi

# ---- runs ------------------------------------------------------------------
N=$(bcftools index -n in.bcf)
digest() { bcftools view -H $1 | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12; }
vc() {  # arm cache [--uncached]
  local C=$W/caches/$2/cache/$2
  rm -rf st out.bcf*
  local s=$(date +%s.%N)
  env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $C -i in.bcf -o out.bcf \
    -y $C/params.snapshot.yaml --statistics light --stats-dir st -v ${3:-} > log.txt 2>&1 || { echo "FAILED $1"; tail -5 log.txt; return; }
  local e=$(echo "$(date +%s.%N) - $s" | bc)
  local miss=$(grep -ho "missing_variants: [0-9]*" -r st | head -1 | cut -d' ' -f2)
  [ -e dig.$1 ] || digest out.bcf > dig.$1
  printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$N" "${miss:-NA}" "$e" "$(cat dig.$1)" | tee -a $OUT/runs.tsv
}
: > $OUT/runs.tsv
for rep in 1 2 3; do
  vc uncached af001 --uncached
  for k in af001 af001_narrow cohort cohort_narrow rerun; do vc $k $k; done
done
# fastVEP alone, for reference
for i in 1 2; do s=$(date +%s.%N); anno in.bcf d.bcf $M16 8; echo "fastvep_alone $(echo "$(date +%s.%N) - $s" | bc)" | tee -a $OUT/runs.tsv; done
echo MAIN DONE
# the 0.01% blueprint only on request (AF0001=1): least important, and its cache takes longest to build
if [ "${AF0001:-0}" = 1 ]; then
  blueprint af0001
  mkcache af0001 $W/c_af0001.bcf $M16; mkcache af0001_narrow $W/c_af0001.bcf $M12
  for rep in 1 2 3; do vc af0001 af0001; vc af0001_narrow af0001_narrow; done
fi
echo DONE
