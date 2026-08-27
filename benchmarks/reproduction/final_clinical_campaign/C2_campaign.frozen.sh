#!/bin/bash
#SBATCH --job-name=C2-campaign
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --time=12:00:00
#SBATCH --exclusive
#SBATCH --array=0-11%6
#SBATCH --output=/results/campaign_logs/C2-%A_%a.out
#
# One campaign cell block: one genome, every assay, direct and cached, for a
# given tool/level/assembly. Serial invocation throughout -- no contig splitting,
# because splitting gives the direct arm a parallel speedup the smaller cached
# arm cannot use and understates the result.
set -uo pipefail
TOOL=${TOOL:?TOOL required (snpeff|fastvep)}
LEVEL=${LEVEL:?LEVEL required (core|clinical)}
ASM=${ASM:?ASM required (GRCh38|GRCh37)}
VC=/mnt/data/vcfcache/.venv/bin/vcfcache
DIST=/results/snpeff/dist
FV=/mnt/data/apps/fastvep/0.3.0-publication
CLIN=/results/clinical_sources
REG=/results/richness/regions
BENCH=/mnt/data/vcfcache_benchmarks
BS=/results/snpeff/scripts/blocksort.awk
source /results/benchmark_archive/archive_lib.sh

RUN_ID=${TOOL}-${LEVEL}-${ASM,,}-${RUN_TAG:-campaign}-serial
OUT=$(archive_dir "$RUN_ID")/tasks; mkdir -p "$OUT"
W=/mnt/data/tmp/c2-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID
mkdir -p "$W/tmp" "$W/aux"; export TMPDIR="$W/tmp"
trap 'rm -rf "$W"' EXIT

if [[ "$ASM" == "GRCh38" ]]; then
  SAMPLES=(KPGP-00128 KPGP-00139 KPGP-00231 KPGP-00320 KPGP-00325 KPGP-00327
           LP6005442-DNA_B04 LP6005443-DNA_B07 LP6005519-DNA_A06 LP6005519-DNA_B04
           LP6005592-DNA_C05 LP6005619-DNA_B01)
  COHORTS=(kpgp kpgp kpgp kpgp kpgp kpgp sgdp sgdp sgdp sgdp sgdp sgdp)
  SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}; COHORT=${COHORTS[$SLURM_ARRAY_TASK_ID]}
  SRC=$BENCH/external_wgs/samples/GRCh38/$COHORT/$SAMPLE.GRCh38.small_variants.vcf.gz
  ASSAYS=(panel wes wgs)
  GENOME=GRCh38.115
  PANEL_BED=$REG/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed
  WES_BED=$REG/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed
else
  mapfile -t PGP < <(ls "$BENCH/external_wgs/samples/GRCh37/pgp"/*.vcf.gz | head -12)
  SRC=${PGP[$SLURM_ARRAY_TASK_ID]}; COHORT=pgp
  SAMPLE=$(basename "$SRC" | cut -d. -f1)
  ASSAYS=(panel wes wgs)
  GENOME=$(cat "$CLIN/SNPEFF_GRCH37")
  # Lifted from the GRCh38 originals with UCSC liftOver; see LIFTOVER_PROVENANCE.txt.
  PANEL_BED=$REG/acmg_sf_v3.3.ensembl115_lifted.mane_select_cds_pad20_GRCh37.bed
  WES_BED=$REG/twist_human_core_exome_hg19.chr1-22-X.pad125.merged_GRCh37.bed
fi
[[ -s "$SRC" ]] || { echo "MISSING INPUT $SRC"; exit 1; }
[[ -n "${ASSAYS_OVERRIDE:-}" ]] && { read -r -a ASSAYS <<< "$ASSAYS_OVERRIDE"; }
for b in "$PANEL_BED" "$WES_BED"; do
  printf '%s\n' "${ASSAYS[@]}" | grep -qE '^(panel|wes)$' && { [[ -s "$b" ]] || { echo "MISSING REGIONS $b"; exit 1; }; }
done

case "$TOOL/$LEVEL" in
  snpeff/core)     [[ "$ASM" == GRCh38 ]] \
                     && ROOT=/results/snpeff/cache_af001/cache/snpeff_af001 \
                     || ROOT=/results/snpeff/caches/snpeff_core_grch37_af001/cache/snpeff_core_grch37_af001
                   QF='%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' ;;
  snpeff/clinical) N=snpeff_clinical_${ASM,,}_af001
                   ROOT=/results/snpeff/caches/$N/cache/$N
                   QF='%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t%INFO/CLNSIG\t%INFO/AF\n'
                   export SNPSIFT_JAR="$DIST/SnpSift.jar"
                   export CLINVAR_VCF="$CLIN/clinvar_${ASM}.vcf.gz"
                   export GNOMAD_VCF="$CLIN/gnomad_exomes_${ASM}.vcf.gz" ;;
  fastvep/clinical) N=fastvep_clinical_${ASM,,}_af001
                   ROOT=/results/richness/caches/$N/cache/$N
                   QF='%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\t%INFO/FV_CLINVAR\t%INFO/FV_GNOMAD\n' ;;
  *) echo "UNSUPPORTED $TOOL/$LEVEL"; exit 1 ;;
esac
[[ -s "$ROOT/vcfcache_annotated.bcf" ]] || { echo "CACHE MISSING $ROOT"; exit 1; }
[[ -s "$BS" ]] || { echo "MISSING blocksort.awk"; exit 1; }

if [[ "$TOOL" == snpeff ]]; then
  export SNPEFF_HEAP=8g
  sed -e "s#^temp_dir:.*#temp_dir: $W/tmp#" -e "s#^snpeff_genome:.*#snpeff_genome: $GENOME#" \
      -e "s#^genome_build:.*#genome_build: $ASM#" \
      "$DIST/params_snpeff_GRCh38.yaml" > "$W/params.yaml"
  grep -q snpeff_clinical_run "$W/params.yaml" \
    || echo "snpeff_clinical_run: $DIST/snpeff_clinical_run.sh" >> "$W/params.yaml"
  if [[ "$LEVEL" == clinical ]]; then
    direct_run () { "$DIST/snpeff_clinical_run.sh" "$1" "$2" "$W/aux" "$DIST/jre/bin/java" \
      "$DIST/snpEff/snpEff.jar" "$DIST/snpEff/snpEff.config" "$DIST/snpEff/data" "$GENOME"; }
  else
    direct_run () { "$DIST/snpeff_parallel.sh" "$1" "$2" "$W/aux" "$DIST/jre/bin/java" \
      "$DIST/snpEff/snpEff.jar" "$DIST/snpEff/snpEff.config" "$DIST/snpEff/data" "$GENOME"; }
  fi
else
  sed -e "s#^temp_dir:.*#temp_dir: $W/tmp#" "$FV/config/params_$ASM.yaml" > "$W/params.yaml"
  TC=$(grep '^transcript_cache:' "$W/params.yaml" | awk '{print $2}')
  FA=$(grep '^fasta:' "$W/params.yaml" | awk '{print $2}')
  direct_run () { /results/richness/scripts/fastvep_run.sh "$1" "$2" "$W/aux" \
    "$FV/bin/fastvep" "$TC" "$FA" "/results/richness/sa/clinical_$ASM"; }
fi

digest () {
  local src="$1" td; td=$(mktemp -d "$W/dg.XXXXXX")
  bcftools query -f "$QF" "$src" 2>/dev/null | awk -f "$BS" \
    | awk -v d="$td" '{ f=d"/"$1; print >> f }'
  for c in $(ls "$td" 2>/dev/null | sort); do
    printf '%s %s\n' "$c" "$(md5sum < "$td/$c" | cut -d' ' -f1)"
  done | md5sum | cut -d' ' -f1
  rm -rf "$td"
}

write_manifest "$(dirname "$OUT")" "{
  \"run_id\": \"$RUN_ID\", \"tool\": \"$TOOL\", \"level\": \"$LEVEL\",
  \"assembly\": \"$ASM\", \"genome_db\": \"$GENOME\", \"cache\": \"$ROOT\",
  \"assays\": [\"${ASSAYS[*]}\"], \"genomes\": 12, \"cohorts\": [\"$COHORT\"],
  \"invocation\": \"serial (documented single-process form)\",
  \"slurm_job\": ${SLURM_ARRAY_JOB_ID:-0}, \"status\": \"primary\",
  \"note\": \"Direct arm calls the annotator bare, without vcfcache in the path. Concordance uses a CHROM+POS block-local digest; a plain order-sensitive digest fails on same-position record reordering and an external sort outruns the annotation at WGS scale.\"
}"

TSV=$OUT/${SAMPLE}.tsv
printf 'sample\tcohort\ttool\tlevel\tassembly\tassay\tcondition\twall_seconds\tinput_records\toutput_records\tmissing_variants\tdigest\tconcordant\n' > "$TSV"

for ASSAY in "${ASSAYS[@]}"; do
  case $ASSAY in
    panel) bcftools view --regions-file "$PANEL_BED" "$SRC" -Ou ;;
    wes)   bcftools view --regions-file "$WES_BED" "$SRC" -Ou ;;
    wgs)   bcftools view "$SRC" -Ou ;;
  esac | bcftools norm -m- -Ou | bcftools view -e 'ALT="*"' -Ou \
      | bcftools norm --rm-dup exact -Ob -o "$W/in.bcf" -W --threads 8
  N=$(bcftools index -n "$W/in.bcf.csi")
  NDUP=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$W/in.bcf" | sort | uniq -d | wc -l)
  [[ "$NDUP" -eq 0 ]] || { echo "INPUT STILL HAS $NDUP DUPLICATE KEYS for $SAMPLE/$ASSAY"; continue; }
  echo "===== [$(date -Is)] $SAMPLE/$ASSAY n=$N ====="

  # Cold state dominates a sub-minute cell, so Panel and WES get one untimed
  # pre-run. At WGS scale it costs an entire extra pass and changes nothing.
  [[ "$ASSAY" != wgs ]] && direct_run "$W/in.bcf" "$W/warm.bcf" >/dev/null 2>&1

  s=$(date +%s.%N); direct_run "$W/in.bcf" "$W/d.bcf"; rc=$?; e=$(date +%s.%N)
  [[ $rc -ne 0 ]] && { echo "DIRECT FAILED $SAMPLE/$ASSAY"; continue; }
  DW=$(echo "$e - $s" | bc); DM=$(digest "$W/d.bcf")
  DN=$(bcftools index -n "$W/d.bcf.csi")
  [[ "$DN" == "$N" ]] || { echo "DIRECT RECORD COUNT $DN != INPUT $N for $SAMPLE/$ASSAY"; continue; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\tdirect\t%s\t%s\t%s\t\t%s\t\n' \
    "$SAMPLE" "$COHORT" "$TOOL" "$LEVEL" "$ASM" "$ASSAY" "$DW" "$N" "$DN" "$DM" >> "$TSV"
  echo "[$(date -Is)] direct=${DW}s"

  rm -rf "$W/st" "$W/c.bcf" "$W/c.bcf.csi"; mkdir -p "$W/st"
  s=$(date +%s.%N)
  $VC annotate -a "$ROOT" -i "$W/in.bcf" -o "$W/c.bcf" --stats-dir "$W/st" -y "$W/params.yaml" -q 2>&1 | tail -2
  rc=${PIPESTATUS[0]}
  e=$(date +%s.%N)
  [[ $rc -ne 0 ]] && { echo "CACHED FAILED $SAMPLE/$ASSAY rc=$rc"; continue; }
  CW=$(echo "$e - $s" | bc); CM=$(digest "$W/c.bcf")
  CN=$(bcftools index -n "$W/c.bcf.csi" 2>/dev/null || echo 0)
  MISS=$(grep -ohP 'Found \K[0-9]+(?= variants not in cache)' "$W/st"/*/workflow.log 2>/dev/null | head -1)
  if [[ "$CN" != "$N" ]]; then CONC="FAIL_COUNT:$CN/$N"
  elif [[ "$DM" == "$CM" ]]; then CONC=PASS
  else CONC=FAIL_CONTENT; fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\tcached\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$SAMPLE" "$COHORT" "$TOOL" "$LEVEL" "$ASM" "$ASSAY" "$CW" "$N" "$CN" "${MISS:-}" "$CM" "$CONC" >> "$TSV"
  echo "[$(date -Is)] cached=${CW}s miss=${MISS:-NA} concordance=$CONC"
done
cat "$TSV"; echo "[$(date -Is)] done $SAMPLE"
