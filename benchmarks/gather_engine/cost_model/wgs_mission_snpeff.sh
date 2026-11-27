#!/bin/bash
# The mission pipeline with SnpEff 5.4c and SnpSift (snpeff_stack.sh: consequences plus the eight VCF
# sources of the Ensembl VEP 116.3 mission), uncached against vcfcache with the gnomAD AF >= 1%
# blueprint, a cohort cache and a rerun cache, as wgs_mission_vep116.sh, whose prep (the blueprint
# and the sites of the 23 KPGP samples beyond it, split into parts) this script reuses.
#   probe <sample>       chr22 only: the sample and the blueprint through the whole pipeline, the
#                        sample also through SnpEff alone and each SnpSift stage alone; checks that
#                        every record comes out and every source matches somewhere
#   part <set> <k>       annotate one part of the VEP prep (<set> bp or new)
#   measure <sample>     join the parts into the sample's caches on this node, then ARM=uncached: one
#                        direct run; ARM=compare: one run of each arm, outputs compared record by record
#                        (snpeff_compare.py); ARM=assay: exome and panel subsets as in
#                        wgs_mission_vep116.sh; otherwise 3 repeats of the cached arms
# Usage: wgs_mission_snpeff.sh <deployment> probe <sample> | part <set> <k> | [ARM=...] measure <sample>
set -uo pipefail
DEP=$1; MODE=$2
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
CM=/results/fastvep05/cost_model
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
TPL=/results/veprs05/caches/veprs05_clinical
VM=/results/vep116/mission; P=$VM/parts
M=/results/snpeff/mission; MP=$M/parts; mkdir -p $MP
ls /results /results/snpeff $VM $M > /dev/null 2>&1   # clears a stale "no such directory" on some nodes
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname), mode $MODE ${*:3}"
test -e $VM/PREP || { echo "VEP prep missing"; exit 1; }

# node-local copies of the shared sources and their SnpSift indexes (snpeff_sources.sh)
R=/mnt/data/tmp/snpeff-mission; mkdir -p $R/tmp
if [ ! -e $R/SHARED ]; then
  test -e $M/src/READY || { echo "shared sources missing"; exit 1; }
  rm -rf $R/custom $R/tables $R/DATA $R/PLAIN $R/src; mkdir -p $R/src
  s=$(date +%s); cp --preserve=timestamps $M/src/*.vcf* $R/src/ || { echo "FAILED copy"; exit 1; }
  log "sources copied: $(du -sh $R/src | cut -f1), $(( $(date +%s) - s )) s, dbSNP $(cat $M/src/DBSNP_MODE)"; touch $R/SHARED
fi
# SnpSift rebuilds an index that is not newer than its VCF (13 min and more than 8 GB of heap for a
# dbSNP group, which ran the 3 GB stages out of memory: jobs 3664-3669). The copy can give a VCF and
# its index the same time stamp; the indexes were built from these very files (snpeff_sources.sh,
# snpeff_dbsnp_split.sh; a rebuild gave byte-identical files), so they are marked newer.
for x in $R/src/*.sidx; do [ $x -nt ${x%.sidx} ] || touch $x; done
CU=$R/src; T=$R/src
anno() {  # in.bcf out.bcf [region file]
  ERR=$R/anno.$$.err; . $CM/snpeff_stack.sh
  # exit status of every process, to stderr: part new-0 (job 3710) failed with nothing on stderr
  bcftools view ${3:+-R $3 --regions-overlap pos} $1 -Ov \
    | eval "$SPIPE"'; st=("${PIPESTATUS[@]}"); echo "stages exit: ${st[*]}" >> $ERR; for x in "${st[@]}"; do [ $x = 0 ] || exit $x; done' \
    | bcftools view -Ob1 -W --threads 8 -o $2
  local ps=("${PIPESTATUS[@]}"); echo "pipe exit: ${ps[*]}" >> $ERR
  for x in "${ps[@]}"; do [ $x = 0 ] || return 1; done
}
check() {  # out.bcf expected-records: every record, and the header of every stage
  local h=$(bcftools view -h $1) t
  [ "$(bcftools index -n $1)" = "$2" ] || { echo "FAILED: $(bcftools index -n $1) of $2 records"; return 1; }
  for t in ANN CV_CLNSIG gnomAD_AF dbSNP_RS SNV_SpliceAI INDEL_SpliceAI PhyloP TrackA_P TrackB_P; do
    grep -q "^##INFO=<ID=$t," <<< "$h" || { echo "FAILED: no $t in header"; return 1; }
  done
}

if [ $MODE = probe ]; then
  SAMPLE=$3; OUT=$M/probe_$SAMPLE; mkdir -p $OUT; W=$R/probe; mkdir -p $W; cd $W
  bcftools view -r chr22 $KP/$SAMPLE.GRCh38.small_variants.vcf.gz -Ob -W -o in.bcf
  bcftools view -r chr22 $BP -Ob -W -o bp.bcf
  t() { local s=$(date +%s.%N); "$@"; echo "$(date +%s.%N) - $s" | bc; }
  printf "step\trecords\tseconds\n" > $OUT/probe.tsv
  anno in.bcf warm.bcf > /dev/null   # warm the page cache
  for f in in bp; do
    e=$(t anno $f.bcf $f.out.bcf); check $f.out.bcf $(bcftools index -n $f.bcf) || { tail -n 20 $R/anno.*.err; exit 1; }
    printf "pipeline_%s\t%s\t%s\n" $f $(bcftools index -n $f.bcf) $e | tee -a $OUT/probe.tsv
  done
  ERR=$R/probe.err; . $CM/snpeff_stack.sh
  e=$(t bash -c "bcftools view in.bcf -Ov | $SNPEFF > se.vcf"); printf "snpeff\t%s\t%s\n" $(bcftools index -n in.bcf) $e | tee -a $OUT/probe.tsv
  for s in clinvar gnomad dbsnp spliceai_snv spliceai_indel phylop trackA trackB; do
    f=$CU/$s.vcf; [ -e $f.sidx ] || f=$f.gz
    e=$(t bash -c "$JAVA -Xmx6g -jar $SE/SnpSift.jar annotate -noLog -c $SE/snpEff/snpEff.config -noId $f se.vcf > /dev/null 2>> $ERR")
    printf "snpsift_%s\t%s\t%s\n" $s $(bcftools index -n in.bcf) $e | tee -a $OUT/probe.tsv
  done
  for tag in ANN CV_CLNSIG gnomAD_AF dbSNP_RS SNV_SpliceAI INDEL_SpliceAI PhyloP TrackA_P TrackB_P; do
    echo -e "$tag\t$(bcftools view -H -i "INFO/$tag!=\".\"" bp.out.bcf 2>/dev/null | wc -l)"; done | tee $OUT/matches_bp_chr22.tsv
  bcftools view -h bp.out.bcf | grep -E "^##(INFO=<ID=ANN|SnpEff)" | head -3 > $OUT/header.txt
  rm -rf $W; echo DONE; exit 0
fi

if [ $MODE = part ]; then
  S=$3; K=$4
  IN=$([ $S = bp ] && echo $BP || echo $P/new.bcf)
  [ -e $MP/$S-$K.bcf.csi ] && { echo "part $S-$K exists"; exit 0; }
  s=$(date +%s); anno $IN $R/$S-$K.bcf $P/$S-$K.regions || { echo "FAILED part $S-$K"; tail -n 5 $R/anno.*.err; exit 1; }
  e=$(( $(date +%s) - s ))
  want=$(bcftools view -H -R $P/$S-$K.regions --regions-overlap pos $IN | wc -l)
  check $R/$S-$K.bcf $want || exit 1
  cp $R/$S-$K.bcf $R/$S-$K.bcf.csi $MP/ && rm -f $R/$S-$K.bcf* && log "part $S-$K: $want records, ${e}s" \
    && echo -e "$S-$K\t$want\t$e\t$(hostname)" >> $M/parts.tsv
  echo DONE; exit 0
fi

[ $MODE = measure ] || { echo "unknown mode $MODE"; exit 1; }
SAMPLE=$3; OUT=$M/wgs_$SAMPLE; mkdir -p $OUT
for s in bp new; do for f in $P/$s-*.regions; do b=$(basename ${f%.regions}); test -e $MP/$b.bcf.csi || { echo "missing $b.bcf"; exit 1; }; done; done
if [ ! -e $R/c_af001.bcf.csi ]; then   # one copy per node
  bcftools concat $(ls $MP/bp-*.bcf | sort -t- -k2,2n) -Ob -W -o $R/c_af001.bcf --threads 8 &&
  bcftools concat $(ls $MP/new-*.bcf | sort -t- -k2,2n) -Ob -W -o $R/a_new.bcf --threads 8 || { echo "FAILED join"; exit 1; }
  [ $(bcftools index -n $R/c_af001.bcf) = $(bcftools index -n $BP) ] || { echo "FAILED: blueprint parts incomplete"; exit 1; }
  log "joined: $(bcftools index -n $R/c_af001.bcf) + $(bcftools index -n $R/a_new.bcf) records"
fi
W=$R/$SAMPLE; mkdir -p $W/tmp; cd $W
mkcache() {  # name annotated.bcf
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
if [ ! -e $W/CACHES ]; then
  bcftools view $KP/$SAMPLE.GRCh38.small_variants.vcf.gz -Ob1 -W --threads 8 -o in.bcf
  sites() { bcftools view -G $1 -Ou | bcftools norm -m- -Ou 2>/dev/null | bcftools annotate -x INFO,ID,QUAL,FILTER -Ob -W -o $2; }
  n=0; for f in $KP/*.small_variants.vcf.gz; do case $f in *$SAMPLE*) ;; *) n=$((n+1)); sites $f o$n.bcf & ;; esac; done; wait
  bcftools concat -a -d exact o*.bcf -Ob -W -o others.bcf --threads 8 2>/dev/null
  bcftools isec -n=2 -w1 -c none $R/a_new.bcf others.bcf -Ob -o a_cohort.bcf 2>/dev/null && bcftools index -f a_cohort.bcf
  bcftools concat -a $R/c_af001.bcf a_cohort.bcf -Ob -W -o c_cohort.bcf --threads 8 2>/dev/null
  bcftools concat -a $R/c_af001.bcf $R/a_new.bcf -Ob -W -o c_rerun.bcf --threads 8 2>/dev/null
  ln -f $R/c_af001.bcf c_af001.bcf; ln -f $R/c_af001.bcf.csi c_af001.bcf.csi
  log "cohort: $n other samples, $(bcftools index -n a_cohort.bcf) annotated sites beyond the blueprint; rerun $(bcftools index -n $R/a_new.bcf)"
  for k in af001 cohort rerun; do mkcache $k $W/c_$k.bcf; done
  rm -f o*.bcf o*.bcf.csi; touch $W/CACHES
fi
digest() { bcftools view -H $1 | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12; }
vc() {  # arm cache [--uncached]; input $INP (default in.bcf)
  local C=$W/caches/$2/cache/$2 o=out.$1.$$.bcf I=${INP:-in.bcf}
  rm -rf st.$$ $o*
  local s=$(date +%s.%N)
  env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $C -i $I -o $o \
    -y $C/params.snapshot.yaml --statistics light --stats-dir st.$$ -v ${3:-} > log.$1.txt 2>&1 || { echo "FAILED $1"; tail -n 5 log.$1.txt; return; }
  local e=$(echo "$(date +%s.%N) - $s" | bc)
  local miss=$(grep -ho "missing_variants: [0-9]*" -r st.$$ | head -1 | cut -d' ' -f2)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$(bcftools index -n $I)" "${miss:-NA}" "$e" "$(digest $o)" "$(hostname)" | tee -a ${RUNS:-$OUT/runs.tsv}
  [ -n "${KEEP:-}" ] && mv $o keep.$1.bcf
  rm -rf st.$$ $o*
}
compare_out() { $PY -I $CM/snpeff_compare.py keep.$1.bcf keep.$2.bcf; }
if [ "${ARM:-}" = compare ]; then   # keep the outputs and compare each cache with the direct run
  KEEP=1; RUNS=$OUT/compare_runs.tsv; vc uncached af001 --uncached
  for k in af001 cohort rerun; do vc $k $k; compare_out uncached $k | sed "s/^/$k\t/" | tee -a $OUT/compare.tsv; done
  rm -f keep.*.bcf
elif [ "${ARM:-}" = uncached ]; then vc uncached af001 --uncached
elif [ "${ARM:-}" = assay ]; then
  RG=$CM/regions; RUNS=$OUT/assay_runs.tsv
  bcftools view -R $RG/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed in.bcf -Ob -W -o wes.bcf
  bcftools view -R $RG/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed in.bcf -Ob -W -o panel.bcf
  log "exome $(bcftools index -n wes.bcf) variants, panel $(bcftools index -n panel.bcf)"
  for a in wes panel; do
    INP=$a.bcf; vc $a.warmup af001 --uncached > /dev/null
    for rep in 1 2 3; do
      [ $rep = 1 ] && KEEP=1 || KEEP=
      vc $a.uncached af001 --uncached
      for k in af001 cohort rerun; do vc $a.$k $k
        [ $rep = 1 ] && compare_out $a.uncached $a.$k | sed "s/^/$a\t$k\t/" | tee -a $OUT/assay_compare.tsv; done
      rm -f keep.*.bcf
    done
  done
  rm -f wes.bcf* panel.bcf*
else for rep in 1 2 3; do for k in af001 cohort rerun; do vc $k $k; done; done; fi
echo DONE
