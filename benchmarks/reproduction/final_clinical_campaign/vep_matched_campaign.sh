#!/bin/bash
# One direct/cached matched-clinical VEP observation per sample/assay.
#SBATCH --job-name=vep-clin-campaign
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --time=12:00:00
#SBATCH --exclusive
#SBATCH --array=0-11%6
#SBATCH --output=/results/vep_clinical/logs/campaign-%A_%a.log
set -euo pipefail

ASM=${ASM:?ASM must be GRCh37 or GRCh38}
RUN_ROOT=${RUN_ROOT:?controller-created RUN_ROOT required}
VC=/mnt/data/vcfcache/.venv/bin/vcfcache
PY=/mnt/data/vcfcache/.venv/bin/python
BASE=/mnt/data/vcfcache_benchmarks
REG=/results/richness/regions
COMPARE=/results/reproduction/final_clinical_campaign/vep_compare_cli.py
DIRECT=/results/reproduction/final_clinical_campaign/vep_clinical_run.sh
NAME=vep_clinical_${ASM,,}_af001
CACHE=/results/vep_clinical/databases/$ASM/cache/$NAME
PARAM_SOURCE=$CACHE/params.snapshot.yaml
WORK=/mnt/data/tmp/vep-clin-campaign-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}
TASK_ROOT=$RUN_ROOT/tasks
mkdir -p "$WORK/tmp" "$TASK_ROOT"
export TMPDIR=$WORK/tmp PYTHONPATH=/mnt/data/vcfcache
trap 'rm -rf "$WORK"' EXIT

if [[ "$ASM" == GRCh38 ]]; then
  SAMPLES=(KPGP-00128 KPGP-00139 KPGP-00231 KPGP-00320 KPGP-00325 KPGP-00327
           LP6005442-DNA_B04 LP6005443-DNA_B07 LP6005519-DNA_A06 LP6005519-DNA_B04
           LP6005592-DNA_C05 LP6005619-DNA_B01)
  COHORTS=(kpgp kpgp kpgp kpgp kpgp kpgp sgdp sgdp sgdp sgdp sgdp sgdp)
  SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
  COHORT=${COHORTS[$SLURM_ARRAY_TASK_ID]}
  SRC=$BASE/external_wgs/samples/GRCh38/$COHORT/$SAMPLE.GRCh38.small_variants.vcf.gz
  PANEL_BED=$REG/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed
  WES_BED=$REG/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed
else
  SAMPLES=(hu1D5A29 hu24385B hu365511 hu57D7DD hu6664F3 hu6ABACE
           hu832966 hu876E3D hu9AF7CC huB333AA huC8266A huCDD5EE)
  SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
  COHORT=pgp
  SRC=$BASE/external_wgs/samples/GRCh37/pgp/$SAMPLE.GRCh37.small_variants.vcf.gz
  PANEL_BED=$REG/acmg_sf_v3.3.ensembl115_lifted.mane_select_cds_pad20_GRCh37.bed
  WES_BED=$REG/twist_human_core_exome_hg19.chr1-22-X.pad125.merged_GRCh37.bed
fi
for path in "$SRC" "$PANEL_BED" "$WES_BED" "$CACHE/vcfcache_annotated.bcf" \
            "$CACHE/vcfcache_annotated.bcf.csi" "$PARAM_SOURCE" "$COMPARE" "$DIRECT"; do
  [[ -s "$path" ]] || { echo "missing $path" >&2; exit 1; }
done
sed -e "s#^temp_dir:.*#temp_dir: \"$WORK/tmp\"#" "$PARAM_SOURCE" > "$WORK/params.yaml"

TSV=$TASK_ROOT/$SAMPLE.tsv
printf 'sample\tcohort\ttool\tlevel\tassembly\tassay\tcondition\twall_seconds\tinput_records\toutput_records\tmissing_variants\tdigest\tconcordant\n' > "$TSV"
for ASSAY in panel wes wgs; do
  case "$ASSAY" in
    panel) bcftools view --regions-file "$PANEL_BED" "$SRC" -Ou ;;
    wes) bcftools view --regions-file "$WES_BED" "$SRC" -Ou ;;
    wgs) bcftools view "$SRC" -Ou ;;
  esac \
    | bcftools norm -m- -Ou \
    | bcftools view -e 'ALT="*"' -Ou \
    | bcftools norm --rm-dup exact -Ob -o "$WORK/in.bcf" -W --threads 8
  N=$(bcftools index -n "$WORK/in.bcf")
  [[ "$N" =~ ^[0-9]+$ && "$N" -gt 0 ]] || { echo "bad input count $N"; exit 1; }
  ART=$TASK_ROOT/${SAMPLE}.${ASSAY}; mkdir -p "$ART"

  if [[ "$ASSAY" != wgs ]]; then
    "$DIRECT" "$WORK/in.bcf" "$WORK/warm.bcf" "$WORK/warm_aux" "$ASM"
    rm -f "$WORK/warm.bcf" "$WORK/warm.bcf.csi"; rm -rf "$WORK/warm_aux"
  fi

  START=$(date +%s.%N)
  /usr/bin/time -v -o "$ART/direct.time" \
    "$DIRECT" "$WORK/in.bcf" "$WORK/direct.bcf" "$ART/direct_aux" "$ASM"
  END=$(date +%s.%N); DIRECT_SECONDS=$(echo "$END - $START" | bc)
  DN=$(bcftools index -n "$WORK/direct.bcf")
  [[ "$DN" == "$N" ]] || { echo "direct count $DN != $N"; exit 1; }

  rm -rf "$WORK/st"
  START=$(date +%s.%N)
  /usr/bin/time -v -o "$ART/cached.time" \
    "$VC" annotate -a "$CACHE" -i "$WORK/in.bcf" -o "$WORK/cached.bcf" \
      --stats-dir "$WORK/st" --statistics light -y "$WORK/params.yaml" -q \
      > "$ART/cached.stdout" 2> "$ART/cached.stderr"
  END=$(date +%s.%N); CACHED_SECONDS=$(echo "$END - $START" | bc)
  CN=$(bcftools index -n "$WORK/cached.bcf")
  [[ "$CN" == "$N" ]] || { echo "cached count $CN != $N"; exit 1; }
  MISS=$(grep -ohP 'Found \K[0-9]+(?= variants not in cache)' "$WORK/st"/*/workflow.log | head -1)
  [[ "$MISS" =~ ^[0-9]+$ ]] || { echo "missing count unavailable"; exit 1; }

  "$VC" --version > "$ART/vcfcache.version"
  "$PY" "$COMPARE" "$WORK/direct.bcf" "$WORK/cached.bcf" "$ART/comparison.json"
  read -r DIRECT_DIGEST CACHED_DIGEST <<< "$("$PY" - "$ART/comparison.json" <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))
print(r['uncached_semantic_sha256'], r['cached_semantic_sha256'])
PY
)"
  printf '%s\t%s\tvep\tclinical\t%s\t%s\tdirect\t%s\t%s\t%s\t\t%s\t\n' \
    "$SAMPLE" "$COHORT" "$ASM" "$ASSAY" "$DIRECT_SECONDS" "$N" "$DN" "$DIRECT_DIGEST" >> "$TSV"
  printf '%s\t%s\tvep\tclinical\t%s\t%s\tcached\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
    "$SAMPLE" "$COHORT" "$ASM" "$ASSAY" "$CACHED_SECONDS" "$N" "$CN" "$MISS" "$CACHED_DIGEST" >> "$TSV"
  rm -f "$WORK/direct.bcf" "$WORK/direct.bcf.csi" "$WORK/cached.bcf" "$WORK/cached.bcf.csi"
done
echo "completed $SAMPLE $ASM"
