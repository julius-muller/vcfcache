#!/bin/bash
# Direct matched-clinical VEP invocation used by the publication campaign.
set -euo pipefail
INPUT=$1
OUTPUT=$2
AUXILIARY=$3
ASM=$4
mkdir -p "$AUXILIARY"
bcftools view "$INPUT" \
  | apptainer exec -B /mnt/data:/mnt/data /mnt/data/apps/ensembl-vep/115/vep.sif vep \
      --offline --buffer_size 500000 --fork 8 --cache \
      --dir_cache /mnt/data/apps/ensembl-vep/115/cachedir \
      --format vcf --vcf --force_overwrite -i stdin -o stdout -a "$ASM" \
      --hgvs --check_existing --af_gnomade \
      --no_stats \
      2> "$AUXILIARY/vep_stderr.txt" \
  | bcftools view -Ob -o "$OUTPUT" -W --threads 8
