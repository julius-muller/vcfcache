#!/bin/bash
# Annotate a BCF with SnpEff exactly as the SnpEff documentation prescribes:
# one process, one database load, streaming VCF in and out.
set -euo pipefail
IN=$1; OUT=$2; AUX=$3; JAVA=$4; JAR=$5; CONFIG=$6; DATA=$7; GENOME=$8
HEAP=${SNPEFF_HEAP:-8g}
THREADS=${SNPEFF_BCFTOOLS_THREADS:-8}
mkdir -p "$AUX"
RECORDS=$(bcftools index -s "$IN" 2>/dev/null | awk '{s+=$3} END {print s+0}')
echo "snpeff_run: serial, $RECORDS records, $HEAP heap" >&2
bcftools view "$IN" -Ov \
  | "$JAVA" "-Xmx$HEAP" -jar "$JAR" -noLog -noStats -nodownload \
      -c "$CONFIG" -dataDir "$DATA" "$GENOME" - 2>> "$AUX/snpeff_stderr.txt" \
  | bcftools view -Ob -o "$OUT" -W --threads "$THREADS"
