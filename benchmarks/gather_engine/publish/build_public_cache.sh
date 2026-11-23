#!/bin/bash
# Build a public cache on the gnomAD v4.1 GRCh38 AF >= 1% blueprint (the af001 archive from Zenodo,
# 10.5281/zenodo.18190436) with vcfcache cache-build, then check it on chr22 of one KPGP genome:
# cached against --uncached, record by record (vep116_compare.py).
#   vep      Ensembl VEP 116.3 --everything --hgvsg (as the published VEP 115.2 caches)
#   fastvep  fastVEP 0.5.0 core: consequences and HGVS, no supplementary sources
# Only sources that may be redistributed. Usage: build_public_cache.sh <deployment> vep|fastvep
set -uo pipefail
DEP=$1; TOOL=$2
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
BPROOT=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp/KPGP-00128.GRCh38.small_variants.vcf.gz
OUT=/results/publish/$TOOL; ROOT=$OUT/gnomad_v4.1_GRCh38_joint_af001; RC=$OUT/recipe
ls /results /results/publish > /dev/null 2>&1
vc() { env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" "$@"; }
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname), $TOOL, vcfcache $(vc --version 2>&1)"
mkdir -p $OUT $RC; W=$(mktemp -d /mnt/data/tmp/pub-$TOOL-XXXX); trap 'rm -rf $W' EXIT

if [ $TOOL = vep ]; then
  NAME=vep116.3_e
  VS=/mnt/data/apps/ensembl-vep/116
  cat > $RC/annotation.yaml <<'Y'
## GRCh38 Ensembl VEP 116.3 annotation config (--everything, --hgvsg)
##
## Installation-agnostic (native, Docker or Apptainer): params.yaml sets annotation_tool_cmd.

genome_build: "GRCh38"
annotation_cmd: |
  ${params.bcftools_cmd} view \${INPUT_BCF} | \
  ${params.annotation_tool_cmd} \
    --offline \
    --buffer_size ${params.vep_buffer} \
    --fork ${params.vep_forks} \
    --cache \
    --dir_cache ${params.vep_cache} \
    --format vcf \
    --vcf \
    --force_overwrite \
    -i stdin \
    -o stdout \
    -a ${params.genome_build} \
    --hgvsg \
    --everything \
    --stats_file \${AUXILIARY_DIR}/vep_stats.html \
    --warning_file \${AUXILIARY_DIR}/vep_warnings.txt \
    2> \${AUXILIARY_DIR}/vep_stderr.txt | \
  ${params.bcftools_cmd} view -o \${OUTPUT_BCF} -Ob -W --threads ${params.threads}

must_contain_info_tag: CSQ
required_tool_version: "116.3"
Y
  cat > $RC/params.yaml <<Y
genome_build: "GRCh38"
annotation_tool_cmd: "apptainer exec -B /mnt/data:/mnt/data -B /results:/results $VS/vep.sif vep"
tool_version_command: "apptainer exec -B /mnt/data:/mnt/data $VS/vep.sif vep | awk '/ensembl-vep/ {print \$NF}'"
threads: 8
bcftools_cmd: "bcftools"
temp_dir: "$W"
vep_cache: "$VS/cachedir"
vep_buffer: 5000
vep_forks: 15
Y
else
  NAME=fastvep0.5.0_core
  FV=/mnt/data/apps/fastvep/0.5.0
  cat > $RC/annotation.yaml <<'Y'
## GRCh38 fastVEP 0.5.0 core annotation config: consequences and HGVS, no supplementary sources

genome_build: "GRCh38"
annotation_cmd: |
  ${params.bcftools_cmd} view \${INPUT_BCF} -Ov | \
  env RAYON_NUM_THREADS=${params.rayon_threads} ${params.annotation_tool_cmd} annotate \
    --input - \
    --output - \
    --output-format vcf \
    --transcript-cache ${params.transcript_cache} \
    --fasta ${params.fasta} \
    --hgvs --no-progress \
    2> \${AUXILIARY_DIR}/fastvep_stderr.txt | \
  ${params.bcftools_cmd} view -Ob -o \${OUTPUT_BCF} -W --threads ${params.threads}

must_contain_info_tag: CSQ
required_tool_version: "0.5.0"
optional_checks: {}
Y
  cat > $RC/params.yaml <<Y
genome_build: "GRCh38"
annotation_tool_cmd: "$FV/bin/fastvep"
tool_version_command: "$FV/bin/fastvep --version"
bcftools_cmd: "bcftools"
temp_dir: "$W"
threads: 8
rayon_threads: 16
transcript_cache: "$FV/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache"
fasta: "$FV/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"
optional_checks: {}
Y
fi

if [ ! -e $ROOT/cache/$NAME/.vcfcache_complete ]; then
  rm -rf $ROOT; mkdir -p $ROOT && cp -a $BPROOT/. $ROOT/ && rm -f $ROOT/zenodo_blueprint_provenance.json || { echo FAILED copy; exit 1; }
  s=$(date +%s)
  vc cache-build -d $ROOT -n $NAME -a $RC/annotation.yaml -y $RC/params.yaml > $OUT/build.log 2>&1 || { echo "FAILED build"; tail -20 $OUT/build.log; exit 1; }
  log "built $NAME in $(( $(date +%s) - s )) s: $(bcftools index -n $ROOT/cache/$NAME/vcfcache_annotated.bcf) of $(bcftools index -n $ROOT/blueprint/vcfcache.bcf) records"
fi
C=$ROOT/cache/$NAME
[ "$(bcftools index -n $C/vcfcache_annotated.bcf)" = "$(bcftools index -n $ROOT/blueprint/vcfcache.bcf)" ] || { echo "FAILED: record count"; exit 1; }

# check: chr22 of one genome, cached against direct
cd $W && bcftools view -r chr22 $KP -Ob -W -o in.bcf
for arm in cached uncached; do
  s=$(date +%s.%N)
  vc annotate -a $C -i in.bcf -o $arm.bcf -y $RC/params.yaml --stats-dir $arm.stats $([ $arm = uncached ] && echo --uncached) > $arm.log 2>&1 || { echo "FAILED $arm"; tail -10 $arm.log; exit 1; }
  log "$arm $(echo "$(date +%s.%N) - $s" | bc) s"
done
$PY -I /results/fastvep05/cost_model/vep116_compare.py uncached.bcf cached.bcf | tee $OUT/check_chr22_KPGP-00128.tsv
du -sh $ROOT | tee -a $OUT/check_chr22_KPGP-00128.tsv
echo DONE
