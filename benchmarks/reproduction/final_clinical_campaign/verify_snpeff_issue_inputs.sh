#!/bin/bash
# Verify the compact public inputs prepared for the upstream SnpEff issue.
#SBATCH --job-name=verify-snpeff-issue
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=2
#SBATCH --mem=12000M
#SBATCH --time=00:15:00
#SBATCH --output=/results/diagnostics/snpeff-issue-inputs/verify-%j.log
set -euo pipefail
D=/results/diagnostics/snpeff-issue-inputs
S=/results/snpeff/dist
G=$(cat /results/clinical_sources/SNPEFF_GRCH37)
printf 'file\twarning\tann_sha256\n' > "$D/verification.tsv"
for input in "$D"/*.vcf; do
  output=${input%.vcf}.annotated.vcf
  "$S/jre/bin/java" -Xmx8g -jar "$S/snpEff/snpEff.jar" \
    -noLog -noStats -nodownload -c "$S/snpEff/snpEff.config" \
    -dataDir "$S/snpEff/data" "$G" "$input" > "$output"
  ann=$(bcftools query -i 'POS=11273459 && REF="C" && ALT="T"' \
    -f '%INFO/ANN\n' "$output")
  warning=0
  [[ "$ann" == *WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS* ]] && warning=1
  digest=$(printf '%s' "$ann" | sha256sum | cut -d' ' -f1)
  printf '%s\t%s\t%s\n' "$(basename "$input")" "$warning" "$digest" \
    >> "$D/verification.tsv"
done
cat "$D/verification.tsv"
