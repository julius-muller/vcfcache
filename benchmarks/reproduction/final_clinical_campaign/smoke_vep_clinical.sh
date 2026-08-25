#!/bin/bash
# Untimed recipe gate before the two genome-wide matched VEP cache builds.
#SBATCH --job-name=vep-clin-smoke
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=16000M
#SBATCH --time=00:20:00
#SBATCH --exclusive
#SBATCH --array=0-1
#SBATCH --output=/results/vep_clinical/logs/smoke-%A_%a.log
set -euo pipefail

ASMS=(GRCh37 GRCh38)
ASM=${ASMS[$SLURM_ARRAY_TASK_ID]}
BASE=/mnt/data/vcfcache_benchmarks
REG=/results/richness/regions
WORK=/mnt/data/tmp/vep-clinical-smoke-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}
mkdir -p "$WORK" /results/vep_clinical/smoke
trap 'rm -rf "$WORK"' EXIT

if [[ "$ASM" == GRCh37 ]]; then
  SRC=$BASE/external_wgs/samples/GRCh37/pgp/hu1D5A29.GRCh37.small_variants.vcf.gz
  BED=$REG/acmg_sf_v3.3.ensembl115_lifted.mane_select_cds_pad20_GRCh37.bed
else
  SRC=$BASE/external_wgs/samples/GRCh38/kpgp/KPGP-00128.GRCh38.small_variants.vcf.gz
  BED=$REG/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed
fi
PARAM_SOURCE=$BASE/bundled_zenodo_caches/gnomad_v4.1_${ASM}_joint_af001/cache/vep115.2_everything/params.snapshot.yaml

bcftools view --regions-file "$BED" "$SRC" -Ou \
  | bcftools norm -m- -Ou \
  | bcftools view -e 'ALT="*"' -Ou \
  | bcftools norm --rm-dup exact -Ob -o "$WORK/in.bcf" -W --threads 8

bcftools view "$WORK/in.bcf" \
  | apptainer exec -B /mnt/data:/mnt/data /mnt/data/apps/ensembl-vep/115/vep.sif vep \
  --offline --buffer_size 500000 --fork 8 --cache \
  --dir_cache /mnt/data/apps/ensembl-vep/115/cachedir \
  --format vcf --vcf --force_overwrite -i stdin \
  -o stdout -a "$ASM" --hgvs --check_existing --af_gnomade \
  --no_stats 2> "$WORK/vep.stderr" \
  | bcftools view -Ob -o "/results/vep_clinical/smoke/${ASM}.bcf" -W --threads 8

BCF=/results/vep_clinical/smoke/${ASM}.bcf
HEADER=$(bcftools view -h "$BCF")
for field in CSQ CLIN_SIG gnomADe_AF; do
  grep -q "$field" <<< "$HEADER" || { echo "smoke output has no $field" >&2; exit 1; }
done
bcftools index -n "$BCF"
cp "$WORK/vep.stderr" "/results/vep_clinical/smoke/${ASM}.vep.stderr"
cp "$PARAM_SOURCE" "/results/vep_clinical/smoke/${ASM}.params.source.yaml"
echo "matched-clinical VEP smoke PASS: $ASM"
