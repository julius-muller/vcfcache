#!/usr/bin/env bash
set -euo pipefail

# --- tiny argument parser -------------------------------------------------
while [[ $# -gt 0 ]]; do
  case $1 in
    --gnomad-af)     AF="$2";           shift 2;;
    --cache-dir)     CACHE_DIR="$2";    shift 2;;
    --threads)       THREADS="$2";      shift 2;;
    --tool-version)  TOOL_VER="$2";     shift 2;;
    --cache-name)    CNAME="$2";        shift 2;;
    --genome)        GENOME="$2";       shift 2;;
    *) echo "Unknown arg $1"; exit 1;;
  esac
done
mkdir -p "${CACHE_DIR}"

# --- example: download & filter gnomAD -----------------------------------
G_SRC="/tmp/gnomad.genomes.${GENOME}.vcf.gz"
curl -L "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/${GENOME}/gnomad.genomes.vcf.bgz" \
     -o "${G_SRC}"
tabix -p vcf "${G_SRC}"

bcftools view -i "AF>=${AF}" "${G_SRC}" -Ob -o /tmp/gnomad_af.bcf --threads "${THREADS}"
bcftools index /tmp/gnomad_af.bcf

# --- create the blueprint + annotation -----------------------------------
vcfstash stash-init    --vcf /tmp/gnomad_af.bcf --output "${CACHE_DIR}" \
                       -y ${VCFSTASH_ROOT}/tests/config/test_params.yaml

vcfstash stash-annotate --db "${CACHE_DIR}" --name "${CNAME}" \
        -a ${VCFSTASH_ROOT}/tests/config/test_annotation.config