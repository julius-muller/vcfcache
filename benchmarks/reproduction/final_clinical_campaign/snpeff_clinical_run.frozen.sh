#!/bin/bash
# Level-2 SnpEff pipeline: consequences, then ClinVar, then gnomAD frequencies.
#
# Each stage is a single streaming process, matching the documented SnpEff and
# SnpSift invocations. This is the minimum a clinical pipeline carries, and it
# is the level VEP's --everything already occupies via its offline cache.
set -euo pipefail
IN=$1; OUT=$2; AUX=$3; JAVA=$4; JAR=$5; CONFIG=$6; DATA=$7; GENOME=$8
SNPSIFT=${SNPSIFT_JAR:?SNPSIFT_JAR required}
CLINVAR=${CLINVAR_VCF:?CLINVAR_VCF required}
GNOMAD=${GNOMAD_VCF:?GNOMAD_VCF required}
HEAP=${SNPEFF_HEAP:-8g}
THREADS=${SNPEFF_BCFTOOLS_THREADS:-8}
for f in "$CLINVAR" "$GNOMAD"; do
  [[ -s "$f" ]] || { echo "snpeff_clinical_run: missing source $f" >&2; exit 1; }
done
mkdir -p "$AUX"
RECORDS=$(bcftools index -s "$IN" 2>/dev/null | awk '{s+=$3} END {print s+0}')
[[ "$RECORDS" =~ ^[0-9]+$ ]] || { echo "snpeff_clinical_run: bad record count" >&2; exit 1; }
echo "snpeff_clinical_run: serial, $RECORDS records, $HEAP heap, genome $GENOME" >&2
bcftools view "$IN" -Ov \
  | "$JAVA" "-Xmx$HEAP" -jar "$JAR" -noLog -noStats -nodownload \
      -c "$CONFIG" -dataDir "$DATA" "$GENOME" - 2>> "$AUX/snpeff_stderr.txt" \
  | "$JAVA" "-Xmx$HEAP" -jar "$SNPSIFT" annotate -noLog -c "$CONFIG" "$CLINVAR" - \
      2>> "$AUX/snpsift_clinvar_stderr.txt" \
  | "$JAVA" "-Xmx$HEAP" -jar "$SNPSIFT" annotate -noLog -c "$CONFIG" -info AF,AC,AN "$GNOMAD" - \
      2>> "$AUX/snpsift_gnomad_stderr.txt" \
  | bcftools view -Ob -o "$OUT" -W --threads "$THREADS"

# A stage that dies mid-stream still leaves a well-formed but under-annotated
# BCF, which would be timed and gated as if it were a clinical pass. Assert the
# three fields the level is defined by are actually present.
set +o pipefail
for tag in ANN CLNSIG AF; do
  bcftools view -H "$OUT" 2>/dev/null | grep -qm1 "$tag=" \
    || { echo "snpeff_clinical_run: no $tag in output; a pipeline stage failed" >&2
         tail -5 "$AUX"/*stderr* 2>/dev/null >&2; exit 1; }
done
set -o pipefail
