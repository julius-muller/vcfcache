#!/bin/bash
# The 40-genome campaign: the mission stacks of fastVEP 0.5.0 (11 sources), vep-rs 0.5.0, Ensembl VEP
# 116.3 and SnpEff 5.4c on whole genomes of the 20 KPGP and 20 SGDP evaluation genomes
# (cohort40_prep.sh), uncached against vcfcache with
#   af001     the gnomAD AF >= 1% blueprint (the public cache),
#   cohort    blueprint plus the 6 training genomes (3 KPGP, 3 SGDP): one cache shared by all genomes,
#   rerun     cohort plus the genome itself (no misses),
#   fastVEP also af001_n10 and cohort_n10: misses annotated from 1 kb copies of the sources.
# Every arm of a genome, the direct run included, runs in one job, so on one node; the cached arms run
# first and the direct run last (it then finds the sources in the page cache at least as warm as the
# cached arms did). One repeat: the 3-sample missions varied by a few percent between repeats, and the
# genomes differ in their misses, which this campaign is about. Then every cached output is compared
# with the direct one: digests (INFO sorted) for fastVEP and vep-rs, vep116_compare.py and
# snpeff_compare.py for the two tools whose output varies with the run (HGNC_ID and CLIN_SIG order in
# VEP, INFO tag order and ANN warnings in SnpEff).
# The caches are annotated in parts, as in the VEP 116.3 mission: bp (the blueprint), new (sites of
# the 23 KPGP genomes beyond it) and extra (sites of the 23 SGDP genomes beyond both). VEP 116.3 and
# SnpEff reuse their mission's bp and new parts; each annotator annotates each record on its own, and
# the digests of the cached runs check that.
# Usage: cohort40.sh <deployment> <fastvep|veprs|vep116|snpeff> part <bp|new|extra> <k> | measure <sample>
set -uo pipefail
DEP=$1; AN=$2; MODE=$3
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
CM=/results/fastvep05/cost_model
PR=/results/cohort40/prep; A=/results/cohort40/$AN; PA=$A/parts; mkdir -p $PA
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
ls /results /results/cohort40 $PR $A > /dev/null 2>&1   # clears a stale "no such directory" on some nodes
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname), $AN, mode $MODE ${*:4}"
test -e $PR/READY || { echo "prep missing"; exit 1; }
partdir() { case $AN/$1 in vep116/bp|vep116/new) echo /results/vep116/mission/parts;; snpeff/bp|snpeff/new) echo /results/snpeff/mission/parts;; *) echo $PA;; esac; }
input() { case $1 in bp) echo $BP;; new) echo /results/vep116/mission/parts/new.bcf;; extra) echo $PR/extra.bcf;; esac; }

# ---- the annotators: node-local data, anno <in> <out> <threads> [regions], mkcache <name> <annotated> ----
case $AN in
fastvep)
  FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
  TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
  MS=/results/fastvep05/mission; T1=$MS/tracks1; OLD=/results/fastvep05/caches/fastvep050_clinical_grch38_af001
  R=/mnt/data/tmp/mission-wgsf; M16=$R/sa16; M10=$R/sa10; mkdir -p $R $M16 $M10   # as wgs_mission.sh
  ls /results/fastvep05 $MS /results/fastvep05/mission10 > /dev/null 2>&1
  if [ ! -e $R/SOURCES ]; then
    for c in $(seq 1 22) X Y; do test -e $MS/done/chr$c || { echo "sources of chr$c missing"; exit 1; }; done
    for f in $MS/sa16/*; do case $(basename $f) in tracka_chr*|trackb_chr*) ;; *) cp $f $M16/ || exit 1;; esac; done
    cp $T1/sa16/track[ab].* $M16/ || exit 1; touch $R/SOURCES
  fi
  if [ ! -e $R/SOURCES10 ]; then
    for c in $(seq 1 22) X Y; do test -e /results/fastvep05/mission10/done/chr$c || { echo "1 kb sources of chr$c missing"; exit 1; }; done
    for f in /results/fastvep05/mission10/sa12/*; do case $(basename $f) in tracka_chr*|trackb_chr*|gnomad_genomes_*) ;; *) cp $f $M10/ || exit 1;; esac; done
    cp $T1/sa10/track[ab].* $M10/ || exit 1; touch $R/SOURCES10
  fi
  [ $(ls $M16 | wc -l) = $(ls $M10 | wc -l) ] || { echo "source sets differ"; exit 1; }
  anno() {
    bcftools view ${4:+-R $4 --regions-overlap pos} $1 -Ov | env RAYON_NUM_THREADS=$3 $FV annotate --input - --output - --output-format vcf \
      --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $M16 2> $R/anno.$$.err | bcftools view -Ob -o $2 -W --threads 8
  }
  mkcache() {  # name annotated [sa-dir for the misses]
    local ROOT=$W/caches/$1 C=$W/caches/$1/cache/$1 SA=${3:-$M16}
    rm -rf $ROOT; mkdir -p $C $ROOT/blueprint; cp -r $OLD/workflow $ROOT/
    sed -e "s#--sa-dir /results/richness/sa/clinical_GRCh38#--sa-dir $SA#" $OLD/cache/$(basename $OLD)/annotation.yaml > $C/annotation.yaml
    grep -q -- "--sa-dir $SA" $C/annotation.yaml || { echo "recipe not rewritten"; exit 1; }
    sed "s#^temp_dir:.*#temp_dir: $W/tmp#" $OLD/cache/$(basename $OLD)/params.snapshot.yaml > $C/params.snapshot.yaml
    ln -f $2 $C/vcfcache_annotated.bcf; ln -f $2.csi $C/vcfcache_annotated.bcf.csi
  }
  ARMS="af001 af001_n10 cohort cohort_n10 rerun";;
veprs)
  VR=/results/veprs05/dist/vep-0.5.0-x86_64-unknown-linux-gnu/vep-0.5.0-x86_64-unknown-linux-gnu/vep
  JC=/results/veprs/json_cache/GRCh38_native_20261023
  FA=/mnt/data/apps/ensembl-vep/116/cachedir/homo_sapiens/116_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
  D=/results/veprs05/data; MS=/results/veprs05/mission; TPL=/results/veprs05/caches/veprs05_clinical
  R=/mnt/data/tmp/veprs-mission; mkdir -p $R   # as wgs_mission_veprs.sh
  ls /results/veprs05 $MS > /dev/null 2>&1
  test -e $MS/TABLES || { echo "tables missing"; exit 1; }
  if [ ! -e $R/TABLES ]; then cp -r $MS/tables $R/ || exit 1; touch $R/TABLES; fi
  T=$R/tables
  G=""; for f in AF AC AN nhomalt AF_afr AF_amr AF_asj AF_eas AF_fin AF_mid AF_nfe AF_remaining AF_sas; do G="$G,gnomAD_$f"; done
  PL="--plugin dbNSFP,$T/clinvar.tsv.gz,CV_CLNSIG,CV_REVSTAT,CV_CLNDN --plugin dbNSFP,$T/gnomad.tsv.gz$G --plugin dbNSFP,$T/dbsnp.tsv.gz,dbSNP_RS"
  PL="$PL --plugin REVEL,file=$D/new_tabbed_revel_grch38.tsv.gz --plugin AlphaMissense,$D/AlphaMissense_hg38.tsv.gz"
  PL="$PL --plugin SpliceAI,snv=$T/spliceai_snv.vcf.gz,indel=$T/spliceai_indel.vcf.gz --plugin dbNSFP,$T/phylop.tsv.gz,PhyloP"
  PL="$PL --plugin dbNSFP,$T/trackA.tsv.gz,TrackA_P,TrackA_Q --plugin dbNSFP,$T/trackB.tsv.gz,TrackB_P,TrackB_Q"
  PL="$PL --plugin pLI,$D/pLI_values.txt --plugin LoFtool,$D/LoFtool_scores.txt --plugin GWAS,file=$D/gwas-catalog-download-associations-alt-full.tsv"
  FLAGS="--hgvs --fasta $FA $PL"
  anno() {
    bcftools view ${4:+-R $4 --regions-overlap pos} $1 -Ov | $VR --assembly GRCh38 --json_cache $JC --fork $3 --vcf --no_stats --format vcf \
      --input_file STDIN --output_file STDOUT $FLAGS 2> $R/anno.$$.err | bcftools view -Ob1 -W --threads 8 -o $2
  }
  mkcache() {
    local ROOT=$W/caches/$1 C=$W/caches/$1/cache/$1
    rm -rf $ROOT; mkdir -p $C $ROOT/blueprint; cp -r $TPL/workflow $ROOT/
    cat > $C/annotation.yaml <<Y
annotation_cmd: \${params.bcftools_cmd} view \${INPUT_BCF} -Ov | \${params.annotation_tool_cmd}
  --assembly GRCh38 --json_cache \${params.transcript_cache} --fork \${params.threads}
  --vcf --no_stats --format vcf --input_file STDIN --output_file STDOUT $FLAGS | \${params.bcftools_cmd}
  view -Ob1 -W --threads \${params.threads} -o \${OUTPUT_BCF}
must_contain_info_tag: CSQ
genome_build: GRCh38
required_tool_version: '116.2'
Y
    sed -e "s#^temp_dir:.*#temp_dir: $W/tmp#" -e "s#^threads:.*#threads: 8#" $TPL/cache/$(basename $TPL)/params.snapshot.yaml > $C/params.snapshot.yaml
    ln -f $2 $C/vcfcache_annotated.bcf; ln -f $2.csi $C/vcfcache_annotated.bcf.csi
  }
  ARMS="af001 cohort rerun";;
vep116)
  VM=/results/vep116/mission; TPL=/results/veprs05/caches/veprs05_clinical
  R=/mnt/data/tmp/vep116-mission; mkdir -p $R/tmp   # as wgs_mission_vep116.sh
  ls /results/vep116 $VM > /dev/null 2>&1
  if [ ! -e $R/DATA ]; then
    test -e $VM/CUSTOM || { echo "custom VCFs missing"; exit 1; }
    mkdir -p $R/custom $R/tables $R/vepdir && cp $VM/custom/*.vcf.gz* $R/custom/ \
      && cp /results/veprs05/mission/tables/spliceai_*.vcf.gz* $R/tables/ || { echo "copy failed"; exit 1; }
    touch $R/DATA
  fi
  CU=$R/custom; T=$R/tables; VDIR=$R/vepdir
  anno() {
    FORKS=$3; . $CM/vep116_stack.sh
    bcftools view ${4:+-R $4 --regions-overlap pos} $1 -Ov | $VEP $VFLAGS -i STDIN -o STDOUT 2> $R/anno.$$.err | bcftools view -Ob1 -W --threads 8 -o $2
  }
  if [ ! -e $R/GWAS ]; then   # the GWAS plugin writes its processed file into --dir on first use
    bcftools view -r chr22:20000000-20020000 $BP -Ob -W -o $R/h.bcf
    anno $R/h.bcf $R/h.out.bcf 1 && [ $(bcftools index -n $R/h.out.bcf) = $(bcftools index -n $R/h.bcf) ] && touch $R/GWAS \
      || { echo "GWAS priming failed"; tail -5 $R/anno.*.err; exit 1; }
  fi
  mkcache() {
    local ROOT=$W/caches/$1 C=$W/caches/$1/cache/$1
    rm -rf $ROOT; mkdir -p $C $ROOT/blueprint; cp -r $TPL/workflow $ROOT/
    FORKS='${params.threads}'; . $CM/vep116_stack.sh
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
  compare() { $PY -I $CM/vep116_compare.py $1 $2; }
  ARMS="af001 cohort rerun";;
snpeff)
  M=/results/snpeff/mission; TPL=/results/veprs05/caches/veprs05_clinical
  R=/mnt/data/tmp/snpeff-mission; mkdir -p $R/tmp   # as wgs_mission_snpeff.sh
  ls /results/snpeff $M > /dev/null 2>&1
  if [ ! -e $R/SHARED ]; then
    test -e $M/src/READY || { echo "shared sources missing"; exit 1; }
    rm -rf $R/src; mkdir -p $R/src
    s=$(date +%s); cp --preserve=timestamps $M/src/*.vcf* $R/src/ || { echo "FAILED copy"; exit 1; }
    log "sources copied: $(du -sh $R/src | cut -f1), $(( $(date +%s) - s )) s"; touch $R/SHARED
  fi
  for x in $R/src/*.sidx; do [ $x -nt ${x%.sidx} ] || touch $x; done   # SnpSift rebuilds an index not newer than its VCF
  CU=$R/src; T=$R/src
  anno() {  # the thread count does not apply: one process per stage
    ERR=$R/anno.$$.err; . $CM/snpeff_stack.sh
    bcftools view ${4:+-R $4 --regions-overlap pos} $1 -Ov \
      | eval "$SPIPE"'; st=("${PIPESTATUS[@]}"); echo "stages exit: ${st[*]}" >> $ERR; for x in "${st[@]}"; do [ $x = 0 ] || exit $x; done' \
      | bcftools view -Ob1 -W --threads 8 -o $2
    local ps=("${PIPESTATUS[@]}"); echo "pipe exit: ${ps[*]}" >> $ERR
    for x in "${ps[@]}"; do [ $x = 0 ] || return 1; done
  }
  mkcache() {
    local ROOT=$W/caches/$1 C=$W/caches/$1/cache/$1
    rm -rf $ROOT; mkdir -p $C $ROOT/blueprint; cp -r $TPL/workflow $ROOT/
    ERR='${AUXILIARY_DIR}/snpeff_stderr.txt'; . $CM/snpeff_stack.sh
    cat > $C/annotation.yaml <<Y
annotation_cmd: \${params.bcftools_cmd} view \${INPUT_BCF} -Ov | $SPIPE |
  \${params.bcftools_cmd} view -Ob1 -W --threads \${params.threads} -o \${OUTPUT_BCF}
must_contain_info_tag: ANN
genome_build: GRCh38
required_tool_version: '5.4c'
Y
    cat > $C/params.snapshot.yaml <<Y
genome_build: GRCh38
bcftools_cmd: bcftools
annotation_tool_cmd: $JAVA
tool_version_command: $SNPEFF_VERSION
threads: 8
temp_dir: $W/tmp
optional_checks: {}
Y
    ln -f $2 $C/vcfcache_annotated.bcf; ln -f $2.csi $C/vcfcache_annotated.bcf.csi
  }
  compare() { $PY -I $CM/snpeff_compare.py $1 $2; }
  ARMS="af001 cohort rerun";;
*) echo "unknown annotator $AN"; exit 1;;
esac

if [ $MODE = part ]; then
  S=$4; K=$5; IN=$(input $S); RG=$PR/$S-$K.regions; test -e $RG || { echo "no $RG"; exit 1; }
  [ -e $(partdir $S)/$S-$K.bcf.csi ] && { echo "part $S-$K exists"; exit 0; }
  s=$(date +%s); anno $IN $R/c40.$S-$K.bcf 16 $RG || { echo "FAILED part $S-$K"; tail -n 5 $R/anno.$$.err; exit 1; }
  e=$(( $(date +%s) - s )); n=$(bcftools index -n $R/c40.$S-$K.bcf)
  want=$(bcftools view -H -R $RG --regions-overlap pos $IN | wc -l)
  [ "$n" = "$want" ] || { echo "FAILED part $S-$K: $n of $want records"; exit 1; }
  cp $R/c40.$S-$K.bcf $PA/$S-$K.bcf && cp $R/c40.$S-$K.bcf.csi $PA/$S-$K.bcf.csi && rm -f $R/c40.$S-$K.bcf* \
    && log "part $S-$K: $n records, ${e}s" && echo -e "$S-$K\t$n\t$e\t$(hostname)" > $PA/$S-$K.tsv
  echo DONE; exit 0
fi

[ $MODE = measure ] || { echo "unknown mode $MODE"; exit 1; }
SAMPLE=$4; OUT=$A/wgs_$SAMPLE; mkdir -p $OUT
IN=$(awk -F'\t' -v s=$SAMPLE '$2 == s && $5 == "GRCh38" {print $6}' $PR/samples.tsv); test -n "$IN" || { echo "unknown sample $SAMPLE"; exit 1; }
for s in bp new extra; do for f in $PR/$s-*.regions; do b=$(basename ${f%.regions})
  test -e $(partdir $s)/$b.bcf.csi || { echo "missing part $b"; exit 1; }; done; done
# ---- caches: the blueprint, the annotated sites beyond it and the cohort cache, once per node ----
L=$R/c40; mkdir -p $L
if [ ! -e $L/c_cohort.bcf.csi ]; then
  j() { local s=$1 o=$2; bcftools concat $(ls $(partdir $s)/$s-*.bcf | sort -t- -k2,2n) -Ob -W -o $o --threads 8; }
  j bp $L/c_af001.bcf && j new $L/a_new.bcf && j extra $L/a_extra.bcf || { echo "FAILED join"; exit 1; }
  bcftools concat -a $L/a_new.bcf $L/a_extra.bcf -Ob -W -o $L/a_all.bcf --threads 8 2> /dev/null
  [ $(bcftools index -n $L/c_af001.bcf) = $(bcftools index -n $BP) ] || { echo "FAILED: blueprint parts incomplete"; exit 1; }
  [ $(bcftools index -n $L/a_all.bcf) = $(bcftools index -n $PR/new46.bcf) ] || { echo "FAILED: parts beyond the blueprint incomplete"; exit 1; }
  bcftools isec -n=2 -w1 -c none $L/a_all.bcf $PR/train.bcf -Ob -o $L/a_train.bcf 2> /dev/null && bcftools index -f $L/a_train.bcf
  bcftools concat -a $L/c_af001.bcf $L/a_train.bcf -Ob -W -o $L/c_cohort.part.bcf --threads 8 2> /dev/null \
    && mv $L/c_cohort.part.bcf.csi $L/c_cohort.bcf.csi.tmp && mv $L/c_cohort.part.bcf $L/c_cohort.bcf && mv $L/c_cohort.bcf.csi.tmp $L/c_cohort.bcf.csi
  rm -f $L/a_new.bcf* $L/a_extra.bcf*
  log "caches: blueprint $(bcftools index -n $L/c_af001.bcf), beyond it $(bcftools index -n $L/a_all.bcf), of the training genomes $(bcftools index -n $L/a_train.bcf)"
fi
W=$R/c40_$SAMPLE; rm -rf $W; mkdir -p $W/tmp; cd $W
bcftools view $IN -Ob1 -W --threads 8 -o in.bcf
# rerun: the cohort cache plus the genome's own sites beyond it
bcftools concat -a -d exact $PR/train.bcf $PR/sites/$SAMPLE.bcf -Ob -W -o ts.bcf 2> /dev/null
bcftools isec -n=2 -w1 -c none $L/a_all.bcf ts.bcf -Ob -o a_rerun.bcf 2> /dev/null && bcftools index -f a_rerun.bcf
bcftools concat -a $L/c_af001.bcf a_rerun.bcf -Ob -W -o c_rerun.bcf --threads 8 2> /dev/null
for k in af001 cohort; do ln -f $L/c_$k.bcf c_$k.bcf; ln -f $L/c_$k.bcf.csi c_$k.bcf.csi; done
for k in af001 cohort rerun; do mkcache $k $W/c_$k.bcf; done
[ $AN = fastvep ] && for k in af001 cohort; do mkcache ${k}_n10 $W/c_$k.bcf $M10; done
log "caches of $SAMPLE: cohort $(bcftools index -n c_cohort.bcf), rerun $(bcftools index -n c_rerun.bcf) records"
# ---- runs ------------------------------------------------------------------
N=$(bcftools index -n in.bcf)
digest() { bcftools view -H $1 | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12; }
MT=$(awk '/MemTotal/ {print $2}' /proc/meminfo)
vc() {  # arm cache [--uncached]: wall, misses, peak memory in use on the node; keeps the output
  local C=$W/caches/$2/cache/$2 o=keep.$1.bcf
  rm -rf st $o*
  (while :; do awk '/MemAvailable/ {print $2}' /proc/meminfo; sleep 2; done) > mem.$1 & local mp=$!
  local s=$(date +%s.%N)
  env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $C -i in.bcf -o $o \
    -y $C/params.snapshot.yaml --statistics light --stats-dir st -v ${3:-} > log.$1.txt 2>&1 || { kill $mp; echo "FAILED $1"; tail -n 5 log.$1.txt; return 1; }
  local e=$(echo "$(date +%s.%N) - $s" | bc); kill $mp
  local miss=$(grep -ho "missing_variants: [0-9]*" -r st | head -1 | cut -d' ' -f2)
  local peak=$(sort -n mem.$1 | head -1 | awk -v t=$MT '{printf "%.1f", (t - $1) / 1048576}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$1" "$N" "${miss:-NA}" "$e" "$(digest $o)" "$(hostname)" "$peak" | tee -a $OUT/runs.tsv
}
: > $OUT/runs.tsv; : > $OUT/compare.tsv
ok=1; for k in $ARMS; do vc $k $k || ok=0; done
vc uncached af001 --uncached || ok=0
# ---- outputs: every cached arm against the direct run --------------------------------
if [ -e keep.uncached.bcf ]; then
  d=$(awk -F'\t' '$2 == "uncached" {print $6}' $OUT/runs.tsv)
  for k in $ARMS; do [ -e keep.$k.bcf ] || continue
    if [ $AN = vep116 ] || [ $AN = snpeff ]; then compare keep.uncached.bcf keep.$k.bcf | sed "s/^/$SAMPLE\t$k\t/" >> $OUT/compare.tsv
    else echo -e "$SAMPLE\t$k\t$([ "$(awk -F'\t' -v k=$k '$2 == k {print $6}' $OUT/runs.tsv)" = "$d" ] && echo identical || echo DIGEST_DIFFERS)" >> $OUT/compare.tsv; fi
  done
  cat $OUT/compare.tsv
fi
cp log.uncached.txt $OUT/ 2> /dev/null
cd $R; rm -rf $W
[ $ok = 1 ] && echo DONE || { echo "FAILED arms"; exit 1; }
