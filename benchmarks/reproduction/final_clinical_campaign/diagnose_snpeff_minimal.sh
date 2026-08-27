#!/bin/bash
# Minimal direct-SnpEff reproducer for context-sensitive transcript warnings.
#SBATCH --job-name=diag-snpeff-minimal
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --time=00:45:00
#SBATCH --exclusive
#SBATCH --output=/results/diagnostics/snpeff-minimal-%j.log
set -euo pipefail

DIST=/results/snpeff/dist
BLUEPRINT=/results/snpeff/caches/snpeff_core_grch37_af001/blueprint/vcfcache.bcf
GENOME=$(cat /results/clinical_sources/SNPEFF_GRCH37)
DEST=/results/diagnostics/snpeff-minimal-${SLURM_JOB_ID}
TARGET_POS=11273459
TARGET_EXPR='CHROM="chr2" && POS=11273459 && REF="C" && ALT="T"'
mkdir -p "$DEST/inputs" "$DEST/outputs" "$DEST/aux" "$DEST/tmp"
export TMPDIR="$DEST/tmp" SNPEFF_HEAP=8g SNPEFF_BCFTOOLS_THREADS=8
printf 'label\tinput_order\tinput_records\twarning_present\tann_sha256\tann\n' \
  > "$DEST/minimal_matrix.tsv"

annotate_bcf () {
  local label=$1 input=$2 order=$3 ann warning digest records
  local output="$DEST/outputs/${label}.bcf"
  records=$(bcftools view -H "$input" | wc -l)
  "$DIST/snpeff_parallel.sh" "$input" "$output" "$DEST/aux/$label" \
    "$DIST/jre/bin/java" "$DIST/snpEff/snpEff.jar" \
    "$DIST/snpEff/snpEff.config" "$DIST/snpEff/data" "$GENOME"
  ann=$(bcftools query -i "$TARGET_EXPR" -f '%INFO/ANN\n' "$output")
  [[ -n "$ann" ]] || { echo "target absent from $label" >&2; exit 1; }
  warning=0; [[ "$ann" == *WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS* ]] && warning=1
  digest=$(printf '%s' "$ann" | sha256sum | cut -d' ' -f1)
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$label" "$order" "$records" "$warning" "$digest" "$ann" \
    >> "$DEST/minimal_matrix.tsv"
}

make_region () {
  local label=$1 region=$2
  bcftools view -r "$region" "$BLUEPRINT" -Ob \
    -o "$DEST/inputs/${label}.bcf" -W --threads 8
}

# Target alone and three sorted two-record controls.
make_region target_only chr2:11273459-11273459
make_region adjacent_insertion_then_target chr2:11273458-11273459
make_region upstream_repeat_then_target chr2:11273436-11273459
make_region target_then_sample_downstream chr2:11273459-11273813
annotate_bcf target_only "$DEST/inputs/target_only.bcf" 'target'
annotate_bcf adjacent_insertion_then_target \
  "$DEST/inputs/adjacent_insertion_then_target.bcf" '11273458:A>AC,target'
annotate_bcf upstream_repeat_then_target \
  "$DEST/inputs/upstream_repeat_then_target.bcf" '11273436 alleles,target'
annotate_bcf target_then_sample_downstream \
  "$DEST/inputs/target_then_sample_downstream.bcf" 'target,11273813:T>C'

# Reverse the minimal adjacent pair. This is deliberately unsorted and is used
# only to determine whether warning emission depends on processing order.
reversed_vcf="$DEST/inputs/target_then_adjacent_insertion.vcf"
bcftools view -h "$BLUEPRINT" > "$reversed_vcf"
bcftools view -H -r chr2:11273459-11273459 "$BLUEPRINT" >> "$reversed_vcf"
bcftools view -H -r chr2:11273458-11273458 "$BLUEPRINT" >> "$reversed_vcf"
output="$DEST/outputs/target_then_adjacent_insertion.bcf"
mkdir -p "$DEST/aux/target_then_adjacent_insertion"
bcftools view -Ov "$reversed_vcf" \
  | "$DIST/jre/bin/java" -Xmx8g -jar "$DIST/snpEff/snpEff.jar" \
      -noLog -noStats -nodownload -c "$DIST/snpEff/snpEff.config" \
      -dataDir "$DIST/snpEff/data" "$GENOME" - \
      2> "$DEST/aux/target_then_adjacent_insertion/snpeff_stderr.txt" \
  | bcftools view -Ob -o "$output"
ann=$(bcftools query -i "$TARGET_EXPR" -f '%INFO/ANN\n' "$output")
warning=0; [[ "$ann" == *WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS* ]] && warning=1
digest=$(printf '%s' "$ann" | sha256sum | cut -d' ' -f1)
printf '%s\t%s\t2\t%s\t%s\t%s\n' target_then_adjacent_insertion \
  'target,11273458:A>AC' "$warning" "$digest" "$ann" \
  >> "$DEST/minimal_matrix.tsv"

printf 'genome\t%s\njob\t%s\nnode\t%s\ncompleted\t%s\n' \
  "$GENOME" "$SLURM_JOB_ID" "$(hostname)" "$(date -Is)" > "$DEST/run.tsv"
echo "minimal diagnostic complete: $DEST"
