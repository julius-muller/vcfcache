#!/bin/bash
# fastVEP annotation, one streaming process, mirroring the recipe the caches were
# built from exactly. The direct arm calls this bare -- without vcfcache in the
# path at all -- so it measures what you would run if you did not have a cache.
set -euo pipefail
IN=$1; OUT=$2; AUX=$3; FVBIN=$4; TCACHE=$5; FASTA=$6; SADIR=${7:-}
THREADS=${FASTVEP_THREADS:-8}
mkdir -p "$AUX"
RECORDS=$(bcftools index -s "$IN" 2>/dev/null | awk '{s+=$3} END {print s+0}')
[[ "$RECORDS" =~ ^[0-9]+$ ]] || { echo "fastvep_run: bad record count" >&2; exit 1; }
SA_ARG=(); [[ -n "$SADIR" ]] && SA_ARG=(--sa-dir "$SADIR")
echo "fastvep_run: serial, $RECORDS records, ${THREADS} rayon threads, sa=${SADIR:-none}" >&2
bcftools view "$IN" -Ov \
  | env RAYON_NUM_THREADS="$THREADS" "$FVBIN" annotate \
      --input - --output - --output-format vcf \
      --transcript-cache "$TCACHE" --fasta "$FASTA" \
      --hgvs --no-progress "${SA_ARG[@]}" 2>> "$AUX/fastvep_stderr.txt" \
  | bcftools view -Ob -o "$OUT" -W --threads "$THREADS"
set +o pipefail
for tag in CSQ ${SADIR:+FV_CLINVAR FV_GNOMAD}; do
  bcftools view -H "$OUT" 2>/dev/null | grep -qm1 "$tag=" \
    || { echo "fastvep_run: no $tag in output; a pipeline stage failed" >&2
         tail -5 "$AUX"/*stderr* 2>/dev/null >&2; exit 1; }
done
