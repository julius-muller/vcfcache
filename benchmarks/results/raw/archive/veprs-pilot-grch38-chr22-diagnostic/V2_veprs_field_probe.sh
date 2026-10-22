#!/bin/bash
# Which annotation fields a Perl->JSON converted vep-rs cache can actually
# populate. Run on the head node; seconds of work, deliberately not a timing
# measurement. Two regions because absence of SIFT/PolyPhen in an intergenic
# region proves nothing -- the second region carries missense consequences.
set -euo pipefail
VEP=/mnt/data/slurm-results/veprs/dist/vep
JC=/mnt/data/slurm-results/veprs/json_cache/GRCh38_chr22
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
W=$(mktemp -d /mnt/data/tmp/veprs-probe-XXXXXX)
trap 'rm -rf "$W"' EXIT

for R in chr22:16000000-16200000 chr22:19700000-20500000; do
  echo "######## region $R"
  bcftools view -r "$R" "$BP" -Ov -o "$W/probe.vcf"
  echo "records: $(grep -vc '^#' "$W/probe.vcf")"
  "$VEP" --assembly GRCh38 --force_overwrite --fork 4 --json_cache "$JC" \
    --input_file "$W/probe.vcf" --output_file "$W/out.txt" --everything 2>&1
  echo "---- consequence classes"
  grep -v '^#' "$W/out.txt" | cut -f7 | sort | uniq -c | sort -rn | awk 'NR<=10'
  echo "---- Extra keys present in the body (all rows)"
  grep -v '^#' "$W/out.txt" | cut -f14 | tr ';' '\n' | sed 's/=.*//' | sort | uniq -c | sort -rn
  echo "---- Extra keys present in missense rows"
  grep -v '^#' "$W/out.txt" | awk -F'\t' '$7 ~ /missense_variant/' | cut -f14 |
    tr ';' '\n' | sed 's/=.*//' | sort | uniq -c | sort -rn
  echo "---- one missense row"
  # awk must not exit early here: under pipefail an early exit SIGPIPEs the
  # upstream grep and the pipeline reports failure. See GUIDELINES.md.
  grep -v '^#' "$W/out.txt" | awk -F'\t' '$7 ~ /missense_variant/ {print $14}' | awk 'NR<=1'
  echo "---- full declared Extra key list (documentation header, not evidence of content)"
  sed -n '/Extra column keys/,/VEP command-line/p' "$W/out.txt"
done
