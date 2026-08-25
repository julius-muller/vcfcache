#!/bin/bash
# Retained-output diagnostic for the four excluded hu6664F3/SnpEff cells.
# This is not a timed publication run.
#SBATCH --job-name=diag-hu6664
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --time=04:00:00
#SBATCH --exclusive
#SBATCH --output=/results/diagnostics/hu6664F3-%j.log
set -euo pipefail

VC=/mnt/data/vcfcache/.venv/bin/vcfcache
DIST=/results/snpeff/dist
CLIN=/results/clinical_sources
REG=/results/richness/regions
SRC=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh37/pgp/hu6664F3.GRCh37.small_variants.vcf.gz
DEST=/results/diagnostics/hu6664F3-${SLURM_JOB_ID}
COMPARE=/results/reproduction/final_clinical_campaign/compare_snpeff_records.py
GENOME=$(cat "$CLIN/SNPEFF_GRCH37")
mkdir -p "$DEST/tmp" "$DEST/aux"
export TMPDIR="$DEST/tmp" SNPEFF_HEAP=8g

bcftools view --regions-file "$REG/twist_human_core_exome_hg19.chr1-22-X.pad125.merged_GRCh37.bed" "$SRC" -Ou \
  | bcftools norm -m- -Ou \
  | bcftools view -e 'ALT="*"' -Ou \
  | bcftools norm --rm-dup exact -Ob -o "$DEST/in.bcf" -W --threads 8

for LEVEL in core clinical; do
  mkdir -p "$DEST/$LEVEL/aux" "$DEST/$LEVEL/st"
  if [[ "$LEVEL" == core ]]; then
    ROOT=/results/snpeff/caches/snpeff_core_grch37_af001/cache/snpeff_core_grch37_af001
    RUN=$DIST/snpeff_parallel.sh
  else
    ROOT=/results/snpeff/caches/snpeff_clinical_grch37_af001/cache/snpeff_clinical_grch37_af001
    RUN=$DIST/snpeff_clinical_run.sh
    export SNPSIFT_JAR="$DIST/SnpSift.jar"
    export CLINVAR_VCF="$CLIN/clinvar_GRCh37.vcf.gz"
    export GNOMAD_VCF="$CLIN/gnomad_exomes_GRCh37.vcf.gz"
  fi
  sed -e "s#^temp_dir:.*#temp_dir: $DEST/tmp#" \
      -e "s#^snpeff_genome:.*#snpeff_genome: $GENOME#" \
      -e 's#^genome_build:.*#genome_build: GRCh37#' \
      "$DIST/params_snpeff_GRCh38.yaml" > "$DEST/$LEVEL/params.yaml"
  grep -q snpeff_clinical_run "$DEST/$LEVEL/params.yaml" \
    || echo "snpeff_clinical_run: $DIST/snpeff_clinical_run.sh" >> "$DEST/$LEVEL/params.yaml"

  "$RUN" "$DEST/in.bcf" "$DEST/$LEVEL/direct.bcf" "$DEST/$LEVEL/aux" \
    "$DIST/jre/bin/java" "$DIST/snpEff/snpEff.jar" "$DIST/snpEff/snpEff.config" \
    "$DIST/snpEff/data" "$GENOME"
  "$VC" annotate -a "$ROOT" -i "$DEST/in.bcf" -o "$DEST/$LEVEL/cached.bcf" \
    --stats-dir "$DEST/$LEVEL/st" -y "$DEST/$LEVEL/params.yaml" -q
  python3 "$COMPARE" "$DEST/$LEVEL/direct.bcf" "$DEST/$LEVEL/cached.bcf" \
    > "$DEST/$LEVEL/comparison.tsv"
  sha256sum "$DEST/$LEVEL/direct.bcf" "$DEST/$LEVEL/cached.bcf" \
    > "$DEST/$LEVEL/sha256.txt"
done

echo "diagnostic outputs retained at $DEST"
