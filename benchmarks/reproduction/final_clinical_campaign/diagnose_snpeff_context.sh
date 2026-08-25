#!/bin/bash
# Untimed direct-SnpEff diagnostic for chr2:11273459 C>T.
#
# The accepted cache contains WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS on two
# C2orf50 ANN entries, whereas direct annotation of hu6664F3 WES does not. This
# matrix runs SnpEff without VCFcache on the allele in progressively larger
# contexts derived independently from the sample and the cache blueprint.
#SBATCH --job-name=diag-snpeff-context
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --output=/results/diagnostics/snpeff-context-%j.log
set -euo pipefail

DIST=/results/snpeff/dist
SAMPLE=/results/diagnostics/hu6664F3-1205/in.bcf
BLUEPRINT=/results/snpeff/caches/snpeff_core_grch37_af001/blueprint/vcfcache.bcf
CACHE=/results/snpeff/caches/snpeff_core_grch37_af001/cache/snpeff_core_grch37_af001/vcfcache_annotated.bcf
RETAINED=/results/diagnostics/hu6664F3-1205/core/direct.bcf
GENOME=$(cat /results/clinical_sources/SNPEFF_GRCH37)
DEST=/results/diagnostics/snpeff-context-${SLURM_JOB_ID}
TARGET=chr2:11273459-11273459
EXPR='REF="C" && ALT="T"'
mkdir -p "$DEST/inputs" "$DEST/outputs" "$DEST/aux" "$DEST/tmp"
export TMPDIR="$DEST/tmp" SNPEFF_HEAP=8g SNPEFF_BCFTOOLS_THREADS=8

for required in "$SAMPLE" "$SAMPLE.csi" "$BLUEPRINT" "$BLUEPRINT.csi" \
                "$CACHE" "$CACHE.csi" "$RETAINED" "$RETAINED.csi"; do
  [[ -s "$required" ]] || { echo "missing required file: $required" >&2; exit 1; }
done

printf 'label\tsource\tregion\tinput_records\twarning_present\tann_sha256\tann\n' \
  > "$DEST/context_matrix.tsv"

record_result () {
  local label=$1 source=$2 region=$3 input_records=$4 output=$5 ann warning digest
  ann=$(bcftools query -r "$TARGET" -i "$EXPR" -f '%INFO/ANN\n' "$output")
  [[ -n "$ann" ]] || { echo "target absent from output: $label" >&2; exit 1; }
  warning=0
  [[ "$ann" == *WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS* ]] && warning=1
  digest=$(printf '%s' "$ann" | sha256sum | cut -d' ' -f1)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$label" "$source" "$region" "$input_records" "$warning" "$digest" "$ann" \
    >> "$DEST/context_matrix.tsv"
}

record_existing () {
  local label=$1 source=$2 output=$3 records
  records=$(bcftools index -n "$output.csi")
  record_result "$label" "$source" existing "$records" "$output"
}

run_region () {
  local label=$1 source_name=$2 source_bcf=$3 region=$4 input output aux records target_count
  input="$DEST/inputs/${label}.bcf"
  output="$DEST/outputs/${label}.bcf"
  aux="$DEST/aux/${label}"
  mkdir -p "$aux"
  bcftools view -r "$region" "$source_bcf" -Ob -o "$input" -W --threads 8
  records=$(bcftools index -n "$input.csi")
  target_count=$(bcftools view -r "$TARGET" -i "$EXPR" -H "$input" | wc -l)
  [[ "$target_count" -eq 1 ]] || {
    echo "$label contains $target_count copies of target" >&2
    exit 1
  }
  "$DIST/snpeff_parallel.sh" "$input" "$output" "$aux" \
    "$DIST/jre/bin/java" "$DIST/snpEff/snpEff.jar" \
    "$DIST/snpEff/snpEff.config" "$DIST/snpEff/data" "$GENOME"
  record_result "$label" "$source_name" "$region" "$records" "$output"
  echo "completed $label records=$records" >&2
}

# Freeze the two observations that motivated the diagnostic.
record_existing retained_direct_wes sample "$RETAINED"
record_existing stored_cache_record blueprint "$CACHE"

# Progressively larger, sorted contexts around the allele. The smallest entry
# contains only records at the target position; the largest sample entry is the
# complete retained WES input. Blueprint contexts use the unannotated blueprint
# that originally produced the cached record.
for source_name in sample blueprint; do
  if [[ "$source_name" == sample ]]; then source_bcf=$SAMPLE; else source_bcf=$BLUEPRINT; fi
  run_region "${source_name}_target" "$source_name" "$source_bcf" "$TARGET"
  run_region "${source_name}_1kb" "$source_name" "$source_bcf" chr2:11272459-11274459
  run_region "${source_name}_10kb" "$source_name" "$source_bcf" chr2:11263459-11283459
  run_region "${source_name}_100kb" "$source_name" "$source_bcf" chr2:11173459-11373459
  run_region "${source_name}_1mb" "$source_name" "$source_bcf" chr2:10273459-12273459
  run_region "${source_name}_10mb" "$source_name" "$source_bcf" chr2:1273459-21273459
done

run_region sample_chr2 sample "$SAMPLE" chr2

sha256sum "$DIST/snpEff/snpEff.jar" "$DIST/snpEff/snpEff.config" \
  > "$DEST/tool.sha256"
"$DIST/jre/bin/java" -jar "$DIST/snpEff/snpEff.jar" -version \
  > "$DEST/snpeff.version.txt" 2>&1 || true
printf 'genome\t%s\njob\t%s\nnode\t%s\ncompleted\t%s\n' \
  "$GENOME" "$SLURM_JOB_ID" "$(hostname)" "$(date -Is)" > "$DEST/run.tsv"
echo "diagnostic complete: $DEST"
