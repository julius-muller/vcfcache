#!/bin/bash
# The mission pipeline on whole genomes with Ensembl VEP 116.3: the eleven sources of the vep-rs
# mission stack (vep116_stack.sh), uncached against vcfcache with the gnomAD AF >= 0.1% blueprint, a
# cohort cache (blueprint plus the other KPGP samples) and a rerun cache (all 23 samples), as
# wgs_mission_veprs.sh. Ensembl VEP needs hours per genome, so the work is split into jobs:
#   prep                 sites to annotate: the blueprint, and every site of the 23 KPGP samples beyond
#                        it (one set for all samples), each split into parts of equal record count
#   part <set> <k>       annotate one part (16 forks); <set> is bp or new
#   measure <sample>     join the parts into the sample's caches on this node (once per node and sample),
#                        then ARM=uncached: one direct run; otherwise 3 repeats of the cached arms
# VEP annotates each record on its own, so the cohort and rerun caches are cut from the annotated
# set beyond the blueprint; the cached-versus-direct digests check this.
# Usage: wgs_mission_vep116.sh <deployment> prep | part <set> <k> | [ARM=uncached] measure <sample>
set -uo pipefail
DEP=$1; MODE=$2
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
TPL=/results/veprs05/caches/veprs05_clinical
M=/results/vep116/mission; P=$M/parts
NPB=${NPB:-12}; NPN=${NPN:-4}   # parts of the blueprint (18.9 M records) and of the set beyond it (5.1 M)
ls /results /results/vep116 $M > /dev/null 2>&1   # clears a stale "no such directory" on some nodes
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname), mode $MODE ${*:3}"

if [ $MODE = prep ]; then
  mkdir -p $P; W=$(mktemp -d /mnt/data/tmp/v116prep-XXXX); trap 'rm -rf $W' EXIT; cd $W
  sites() { bcftools view -G $1 -Ou | bcftools norm -m- -Ou 2>/dev/null | bcftools annotate -x INFO,ID,QUAL,FILTER -Ob -W -o $2; }
  n=0; for f in $KP/*.small_variants.vcf.gz; do n=$((n+1)); sites $f s$n.bcf & done; wait
  bcftools concat -a -d exact s*.bcf -Ob -W -o all.bcf --threads 8 2>/dev/null
  bcftools isec -n=1 -w1 -c none all.bcf $BP -Ob -o $P/new.bcf 2>/dev/null && bcftools index -f $P/new.bcf
  log "$n samples, $(bcftools index -n $P/new.bcf) sites beyond the blueprint"
  # parts by POS (bcftools view -R ... --regions-overlap pos): every record in exactly one part
  split_regions() {  # bcf parts prefix
    local tot=$(bcftools index -n $1)
    bcftools query -f '%CHROM\t%POS\n' $1 | awk -v per=$(( (tot + $2 - 1) / $2 )) -v pre=$3 'BEGIN{OFS="\t"; k=0; n=0}
      { if (n >= per && !($1 == c && $2 == last)) { print c, s, last > (pre "-" k ".regions"); k++; n=0 }   # never cut within one POS
        if (n == 0 || $1 != c) { if (n > 0) print c, s, last > (pre "-" k ".regions"); c=$1; s=$2 }
        last=$2; n++ }
      END { print c, s, last > (pre "-" k ".regions") }'
  }
  split_regions $BP $NPB $P/bp; split_regions $P/new.bcf $NPN $P/new
  for s in bp new; do ls $P/$s-*.regions | wc -l | sed "s/^/$s parts: /"; done
  touch $M/PREP; echo DONE; exit 0
fi

# node-local copies of the custom VCFs, the SpliceAI-style VCFs and the GWAS plugin's --dir
R=/mnt/data/tmp/vep116-mission; mkdir -p $R/tmp
if [ ! -e $R/DATA ]; then
  test -e $M/CUSTOM || { echo "custom VCFs missing"; exit 1; }
  mkdir -p $R/custom $R/tables $R/vepdir && cp /results/vep116/mission/custom/*.vcf.gz* $R/custom/ \
    && cp /results/veprs05/mission/tables/spliceai_*.vcf.gz* $R/tables/ || { echo "copy failed"; exit 1; }
  touch $R/DATA
fi
CU=$R/custom; T=$R/tables; VDIR=$R/vepdir
anno() {  # in.bcf out.bcf forks [region file]
  FORKS=$3; . /results/fastvep05/cost_model/vep116_stack.sh
  bcftools view ${4:+-R $4 --regions-overlap pos} $1 -Ov | $VEP $VFLAGS -i STDIN -o STDOUT 2> $R/anno.$$.err | bcftools view -Ob1 -W --threads 8 -o $2
}
if [ ! -e $R/GWAS ]; then   # the GWAS plugin writes its processed file into --dir on first use
  bcftools view -r chr22:20000000-20020000 $BP -Ob -W -o $R/h.bcf
  anno $R/h.bcf $R/h.out.bcf 1 && [ $(bcftools index -n $R/h.out.bcf) = $(bcftools index -n $R/h.bcf) ] && touch $R/GWAS \
    || { echo "GWAS priming failed"; tail -5 $R/anno.*.err; exit 1; }
fi

if [ $MODE = part ]; then
  S=$3; K=$4; test -e $M/PREP || { echo "prep missing"; exit 1; }
  IN=$([ $S = bp ] && echo $BP || echo $P/new.bcf)
  [ -e $P/$S-$K.bcf.csi ] && { echo "part $S-$K exists"; exit 0; }
  s=$(date +%s); anno $IN $R/$S-$K.bcf 16 $P/$S-$K.regions || { echo "FAILED part $S-$K"; tail -5 $R/anno.*.err; exit 1; }
  n=$(bcftools index -n $R/$S-$K.bcf); e=$(( $(date +%s) - s ))
  want=$(bcftools view -H -R $P/$S-$K.regions --regions-overlap pos $IN | wc -l)
  [ "$n" = "$want" ] || { echo "FAILED part $S-$K: $n of $want records"; exit 1; }
  cp $R/$S-$K.bcf $R/$S-$K.bcf.csi $P/ && rm -f $R/$S-$K.bcf* && log "part $S-$K: $n records, ${e}s" && echo -e "$S-$K\t$n\t$e\t$(hostname)" >> $M/parts.tsv
  echo DONE; exit 0
fi

[ $MODE = measure ] || { echo "unknown mode $MODE"; exit 1; }
SAMPLE=$3; OUT=$M/wgs_$SAMPLE; mkdir -p $OUT
for s in bp new; do for f in $P/$s-*.regions; do test -e ${f%.regions}.bcf.csi || { echo "missing ${f%.regions}.bcf"; exit 1; }; done; done
if [ ! -e $R/c_af001.bcf.csi ]; then   # one copy per node
  bcftools concat $(ls $P/bp-*.bcf | sort -t- -k2,2n) -Ob -W -o $R/c_af001.bcf --threads 8 &&
  bcftools concat $(ls $P/new-*.bcf | sort -t- -k2,2n) -Ob -W -o $R/a_new.bcf --threads 8 || { echo "FAILED join"; exit 1; }
  [ $(bcftools index -n $R/c_af001.bcf) = $(bcftools index -n $BP) ] || { echo "FAILED: blueprint parts incomplete"; exit 1; }
  log "joined: $(bcftools index -n $R/c_af001.bcf) + $(bcftools index -n $R/a_new.bcf) records"
fi
W=$R/$SAMPLE; mkdir -p $W/tmp; cd $W
mkcache() {  # name annotated.bcf
  local ROOT=$W/caches/$1 C=$W/caches/$1/cache/$1
  rm -rf $ROOT; mkdir -p $C $ROOT/blueprint; cp -r $TPL/workflow $ROOT/
  FORKS='${params.threads}'; . /results/fastvep05/cost_model/vep116_stack.sh
  cat > $C/annotation.yaml <<Y
annotation_cmd: \${params.bcftools_cmd} view \${INPUT_BCF} -Ov | \${params.annotation_tool_cmd} $VFLAGS
  -i STDIN -o STDOUT 2> \${AUXILIARY_DIR}/vep_stderr.txt | \${params.bcftools_cmd} view -Ob1 -W --threads \${params.threads} -o \${OUTPUT_BCF}
must_contain_info_tag: CSQ
genome_build: GRCh38
required_tool_version: '116.3'
Y
  cat > $C/params.snapshot.yaml <<Y
genome_build: GRCh38
bcftools_cmd: bcftools
annotation_tool_cmd: $VEP
tool_version_command: $VEP --help | grep ensembl-vep
threads: 8
temp_dir: $W/tmp
optional_checks: {}
Y
  ln -f $2 $C/vcfcache_annotated.bcf; ln -f $2.csi $C/vcfcache_annotated.bcf.csi
}
if [ ! -e $W/CACHES ]; then
  bcftools view $KP/$SAMPLE.GRCh38.small_variants.vcf.gz -Ob1 -W --threads 8 -o in.bcf
  sites() { bcftools view -G $1 -Ou | bcftools norm -m- -Ou 2>/dev/null | bcftools annotate -x INFO,ID,QUAL,FILTER -Ob -W -o $2; }
  n=0; for f in $KP/*.small_variants.vcf.gz; do case $f in *$SAMPLE*) ;; *) n=$((n+1)); sites $f o$n.bcf & ;; esac; done; wait
  bcftools concat -a -d exact o*.bcf -Ob -W -o others.bcf --threads 8 2>/dev/null
  # cohort: the annotated sites beyond the blueprint that the other samples carry
  bcftools isec -n=2 -w1 -c none $R/a_new.bcf others.bcf -Ob -o a_cohort.bcf 2>/dev/null && bcftools index -f a_cohort.bcf
  bcftools concat -a $R/c_af001.bcf a_cohort.bcf -Ob -W -o c_cohort.bcf --threads 8 2>/dev/null
  bcftools concat -a $R/c_af001.bcf $R/a_new.bcf -Ob -W -o c_rerun.bcf --threads 8 2>/dev/null
  ln -f $R/c_af001.bcf c_af001.bcf; ln -f $R/c_af001.bcf.csi c_af001.bcf.csi
  log "cohort: $n other samples, $(bcftools index -n a_cohort.bcf) annotated sites beyond the blueprint; rerun $(bcftools index -n $R/a_new.bcf)"
  for k in af001 cohort rerun; do mkcache $k $W/c_$k.bcf; done
  rm -f o*.bcf o*.bcf.csi; touch $W/CACHES
fi
N=$(bcftools index -n in.bcf)
digest() { bcftools view -H $1 | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12; }
vc() {  # arm cache [--uncached]
  local C=$W/caches/$2/cache/$2 o=out.$1.$$.bcf
  rm -rf st.$$ $o*
  local s=$(date +%s.%N)
  env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $C -i in.bcf -o $o \
    -y $C/params.snapshot.yaml --statistics light --stats-dir st.$$ -v ${3:-} > log.$1.txt 2>&1 || { echo "FAILED $1"; tail -5 log.$1.txt; return; }
  local e=$(echo "$(date +%s.%N) - $s" | bc)
  local miss=$(grep -ho "missing_variants: [0-9]*" -r st.$$ | head -1 | cut -d' ' -f2)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$N" "${miss:-NA}" "$e" "$(digest $o)" "$(hostname)" | tee -a $OUT/runs.tsv
  rm -rf st.$$ $o*
}
if [ "${ARM:-}" = uncached ]; then vc uncached af001 --uncached
else for rep in 1 2 3; do for k in af001 cohort rerun; do vc $k $k; done; done; fi
echo DONE
