#!/bin/bash
# The mission pipeline on one whole genome with vep-rs 0.5.0: the fastVEP
# mission stack's sources as vep-rs plugins (veprs_mission_stack.sh / _join.sh),
# uncached against vcfcache with the gnomAD AF >= 1% blueprint, a cohort cache
# (blueprint plus the other KPGP samples) and a rerun cache (sample included).
# vep-rs reads tabix tables per variant, so there is no narrow-block variant.
# STACK2=1 adds 1000 Genomes (veprs_onekg_table.sh) and dbNSFP 5.4a read directly, twice: SIFT/PolyPhen,
# and GERP/phyloP (vep-rs has no bigwig or region sources, so GERP is only where dbNSFP scores it and
# ENCODE cCREs cannot be added). PREP=1 stops once the caches are built; a second run measures.
# ASSAY=1 measures instead the sample's genome restricted to the Twist core exome (pad 125) and the
# ACMG SF v3.3 CDS (pad 20), as assay_test.sh: a warm-up, then 3 repeats of every arm.
# Usage: [STACK2=1] [PREP=1|ASSAY=1] wgs_mission_veprs.sh <deployment> [sample]
set -uo pipefail
DEP=$1; SAMPLE=${2:-KPGP-00128}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
VR=/results/veprs05/dist/vep-0.5.0-x86_64-unknown-linux-gnu/vep-0.5.0-x86_64-unknown-linux-gnu/vep
JC=/results/veprs/json_cache/GRCh38_native_20261023
FA=/mnt/data/apps/ensembl-vep/116/cachedir/homo_sapiens/116_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
D=/results/veprs05/data
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
TPL=/results/veprs05/caches/veprs05_clinical
MS=/results/veprs05/mission
ls /results /results/veprs05 $MS > /dev/null 2>&1   # clears a stale "no such directory" on some nodes
OUT=/results/veprs05/mission/wgs${STACK2:+2}_$SAMPLE; mkdir -p $OUT
R=/mnt/data/tmp/veprs-mission${STACK2:+2}; mkdir -p $R; cd $R
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname), $($VR --version)"
# ---- tables: node-local copies ------------------------------------------------
test -e $MS/TABLES || { echo "tables missing"; exit 1; }
[ -z "${STACK2:-}" ] || test -e $MS/ONEKG || { echo "1000 Genomes table missing"; exit 1; }
DBN=/results/sources/dbnsfp/dbNSFP5.4a_grch38.gz
if [ ! -e $R/TABLES ]; then cp -r $MS/tables $R/ || exit 1
  if [ -n "${STACK2:-}" ]; then mkdir -m 700 -p $R/dbnsfp && cp $DBN $DBN.tbi $R/dbnsfp/ || exit 1; fi
  touch $R/TABLES
fi
T=$R/tables
G=""; for f in AF AC AN nhomalt AF_afr AF_amr AF_asj AF_eas AF_fin AF_mid AF_nfe AF_remaining AF_sas; do G="$G,gnomAD_$f"; done
PL="--plugin dbNSFP,$T/clinvar.tsv.gz,CV_CLNSIG,CV_REVSTAT,CV_CLNDN --plugin dbNSFP,$T/gnomad.tsv.gz$G --plugin dbNSFP,$T/dbsnp.tsv.gz,dbSNP_RS"
PL="$PL --plugin REVEL,file=$D/new_tabbed_revel_grch38.tsv.gz --plugin AlphaMissense,$D/AlphaMissense_hg38.tsv.gz"
PL="$PL --plugin SpliceAI,snv=$T/spliceai_snv.vcf.gz,indel=$T/spliceai_indel.vcf.gz --plugin dbNSFP,$T/phylop.tsv.gz,PhyloP"
PL="$PL --plugin dbNSFP,$T/trackA.tsv.gz,TrackA_P,TrackA_Q --plugin dbNSFP,$T/trackB.tsv.gz,TrackB_P,TrackB_Q"
PL="$PL --plugin pLI,$D/pLI_values.txt --plugin LoFtool,$D/LoFtool_scores.txt --plugin GWAS,file=$D/gwas-catalog-download-associations-alt-full.tsv"
if [ -n "${STACK2:-}" ]; then D5=$R/dbnsfp/$(basename $DBN)
  PL="$PL --plugin dbNSFP,$T/onekg.tsv.gz,KG_AF,KG_AF_EAS,KG_AF_AMR,KG_AF_EUR,KG_AF_AFR,KG_AF_SAS"
  PL="$PL --plugin dbNSFP,$D5,SIFT4G_score,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred"
  PL="$PL --plugin dbNSFP,$D5,GERP_92_mammals,phyloP100way_vertebrate"
fi
FLAGS="--hgvs --fasta $FA $PL"
anno() {  # in.bcf out.bcf threads
  bcftools view $1 -Ov | $VR --assembly GRCh38 --json_cache $JC --fork $3 --vcf --no_stats --format vcf \
    --input_file STDIN --output_file STDOUT $FLAGS 2> $R/anno.err | bcftools view -Ob1 -W --threads 8 -o $2
}
if [ ! -e $R/c_af001.bcf.csi ]; then
  s=$(date +%s); anno $BP $R/c_af001.part.bcf 16 && mv $R/c_af001.part.bcf $R/c_af001.bcf && mv $R/c_af001.part.bcf.csi $R/c_af001.bcf.csi || exit 1
  log "cache af001: $(bcftools index -n $R/c_af001.bcf) records, $(( $(date +%s) - s ))s, $(du -h $R/c_af001.bcf | cut -f1)"
  for f in $(grep -o "dbNSFP_[A-Za-z_]*\|SpliceAI_[A-Za-z_]*\|REVEL\|am_pathogenicity\|pLI\|LoFtool\|GWAS[A-Za-z_]*" <(bcftools view -h $R/c_af001.bcf | grep "ID=CSQ") | sort -u); do echo -n "$f "; done; echo
  bcftools view -H $R/c_af001.bcf | head -200000 | grep -o "CSQ=[^\t]*" | awk -F'|' '{for(i=1;i<=NF;i++) if($i!="") n[i]++} END{for(i in n) printf "%d:%d ", i, n[i]; print ""}' | tee $OUT/csq_fill.txt
fi
W=$R/$SAMPLE; mkdir -p $W/tmp; cd $W
mkcache() {  # name annotated.bcf
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
if [ ! -e $W/CACHES ]; then
  bcftools view $KP/$SAMPLE.GRCh38.small_variants.vcf.gz -Ob1 -W --threads 8 -o in.bcf
  ln -f $R/c_af001.bcf c_af001.bcf; ln -f $R/c_af001.bcf.csi c_af001.bcf.csi
  sites() { bcftools view -G $1 -Ou | bcftools norm -m- -Ou 2>/dev/null | bcftools annotate -x INFO,ID,QUAL,FILTER -Ob -W -o $2; }
  n=0; for f in $KP/*.small_variants.vcf.gz; do
    case $f in *$SAMPLE*) sites $f self.bcf;; *) n=$((n+1)); sites $f o$n.bcf;; esac
  done
  bcftools concat -a -d exact o*.bcf -Ob -W -o others.bcf --threads 8 2>/dev/null
  bcftools isec -n=1 -w1 -c none others.bcf $BP -Ob -o new_cohort.bcf 2>/dev/null; bcftools index -f new_cohort.bcf
  bcftools concat -a -d exact others.bcf self.bcf -Ob -W -o withself.bcf --threads 8 2>/dev/null
  bcftools isec -n=1 -w1 -c none withself.bcf $BP -Ob -o new_rerun.bcf 2>/dev/null; bcftools index -f new_rerun.bcf
  log "cohort: $n other samples, $(bcftools index -n new_cohort.bcf) sites beyond the blueprint; with the sample itself $(bcftools index -n new_rerun.bcf)"
  for k in cohort rerun; do
    anno new_$k.bcf a_$k.bcf 16
    bcftools concat -a c_af001.bcf a_$k.bcf -Ob -W -o c_$k.bcf --threads 8 2>/dev/null
  done
  for k in af001 cohort rerun; do mkcache $k $W/c_$k.bcf; done
  rm -f o*.bcf o*.bcf.csi; touch $W/CACHES
fi
[ "${PREP:-0}" = 1 ] && { echo PREP DONE; exit 0; }
# ---- runs ------------------------------------------------------------------
N=$(bcftools index -n in.bcf)
digest() { bcftools view -H $1 | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12; }
vc() {  # arm cache [--uncached]; input $INP (default in.bcf), table $RUNS (default runs.tsv)
  local C=$W/caches/$2/cache/$2 I=${INP:-in.bcf}
  rm -rf st out.bcf*
  local s=$(date +%s.%N)
  env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $C -i $I -o out.bcf \
    -y $C/params.snapshot.yaml --statistics light --stats-dir st -v ${3:-} > log.txt 2>&1 || { echo "FAILED $1"; tail -5 log.txt; return; }
  local e=$(echo "$(date +%s.%N) - $s" | bc)
  local miss=$(grep -ho "missing_variants: [0-9]*" -r st | head -1 | cut -d' ' -f2)
  [ -e dig.$1 ] || digest out.bcf > dig.$1
  printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$(bcftools index -n $I)" "${miss:-NA}" "$e" "$(cat dig.$1)" | tee -a ${RUNS:-$OUT/runs.tsv}
}
if [ "${ASSAY:-0}" = 1 ]; then
  RG=/results/fastvep05/cost_model/regions; RUNS=$OUT/assay_runs.tsv; : > $RUNS
  bcftools view -R $RG/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed in.bcf -Ob -W -o wes.bcf
  bcftools view -R $RG/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed in.bcf -Ob -W -o panel.bcf
  log "node $(hostname): exome $(bcftools index -n wes.bcf) variants, panel $(bcftools index -n panel.bcf)" | tee $OUT/assay_info.txt
  rm -f dig.wes.* dig.panel.*
  for a in wes panel; do
    INP=$a.bcf; vc $a.warmup af001 --uncached > /dev/null
    for rep in 1 2 3; do vc $a.uncached af001 --uncached; for k in af001 cohort rerun; do vc $a.$k $k; done; done
  done
  rm -f wes.bcf* panel.bcf*; echo DONE; exit 0
fi
: > $OUT/runs.tsv
for rep in 1 2 3; do
  vc uncached af001 --uncached
  for k in af001 cohort rerun; do vc $k $k; done
done
for i in $([ -z "${STACK2:-}" ] && echo 1 2); do s=$(date +%s.%N); anno in.bcf d.bcf 8; echo "veprs_alone $(echo "$(date +%s.%N) - $s" | bc)" | tee -a $OUT/runs.tsv; done
echo DONE
