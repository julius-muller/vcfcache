#!/bin/bash
# Freeze versions and hashes for every final-campaign annotator dependency.
#SBATCH --job-name=freeze-tools
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=2
#SBATCH --mem=4000M
#SBATCH --time=02:00:00
#SBATCH --output=/results/reproduction/final_clinical_campaign/toolchain/capture-%j.log
set -euo pipefail

OUT=/results/reproduction/final_clinical_campaign/toolchain
FV=/mnt/data/apps/fastvep/0.3.0-publication
SN=/results/snpeff/dist
mkdir -p "$OUT"

{
  date -Is
  hostname
  uname -a
  bcftools --version | sed -n '1,2p'
  apptainer --version
  "$FV/bin/fastvep" --version
  "$SN/jre/bin/java" -version 2>&1 | sed -n '1,3p'
  "$SN/jre/bin/java" -jar "$SN/snpEff/snpEff.jar" -version 2>&1 | sed -n '1,3p'
  apptainer exec /mnt/data/apps/ensembl-vep/115/vep.sif vep 2>&1 | grep ensembl-vep
  git -C /mnt/data/vcfcache rev-parse HEAD
  git -C /mnt/data/vcfcache status --short
  /mnt/data/vcfcache/.venv/bin/vcfcache --version
} > "$OUT/versions.txt"

FILES=(
  /usr/local/bin/bcftools
  /mnt/data/apps/ensembl-vep/115/vep.sif
  "$FV/bin/fastvep"
  "$FV/data/GRCh37/Homo_sapiens.GRCh37.fastvep.cache"
  "$FV/data/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
  "$FV/data/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai"
  "$FV/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache"
  "$FV/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"
  "$FV/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.fai"
  /results/richness/sa/clinical_GRCh37/clinvar.osa2
  /results/richness/sa/clinical_GRCh37/gnomad.osa2
  /results/richness/sa/clinical_GRCh38/clinvar.osa2
  /results/richness/sa/clinical_GRCh38/gnomad.osa2
  "$SN/snpEff/snpEff.jar"
  "$SN/SnpSift.jar"
  "$SN/snpEff/snpEff.config"
  "$SN/snpEff/data/GRCh37.75/snpEffectPredictor.bin"
  "$SN/snpEff/data/GRCh38.115/snpEffectPredictor.bin"
  /results/clinical_sources/clinvar_GRCh37.vcf.gz
  /results/clinical_sources/clinvar_GRCh37.vcf.gz.tbi
  /results/clinical_sources/clinvar_GRCh38.vcf.gz
  /results/clinical_sources/clinvar_GRCh38.vcf.gz.tbi
  /results/clinical_sources/gnomad_exomes_GRCh37.vcf.gz
  /results/clinical_sources/gnomad_exomes_GRCh37.vcf.gz.tbi
  /results/clinical_sources/gnomad_exomes_GRCh38.vcf.gz
  /results/clinical_sources/gnomad_exomes_GRCh38.vcf.gz.tbi
)
for path in "${FILES[@]}"; do
  [[ -f "$path" ]] || { echo "missing immutable asset: $path" >&2; exit 1; }
done
sha256sum "${FILES[@]}" > "$OUT/assets.sha256"
stat --printf='%s\t%n\n' "${FILES[@]}" > "$OUT/assets.bytes.tsv"

# VEP's transcript cache is a directory. Freeze a sorted file-level manifest
# and then hash that manifest, rather than reducing the directory to a path.
find /mnt/data/apps/ensembl-vep/115/cachedir -type f -print0 \
  | sort -z \
  | xargs -0 sha256sum > "$OUT/vep_cache_files.sha256"
sha256sum "$OUT/vep_cache_files.sha256" > "$OUT/vep_cache_manifest.sha256"
echo "toolchain capture complete"
