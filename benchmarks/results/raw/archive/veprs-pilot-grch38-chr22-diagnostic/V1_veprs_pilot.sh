#!/bin/bash
#SBATCH --job-name=V1-veprs-pilot
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=16
#SBATCH --mem=30000M
#SBATCH --time=08:00:00
#SBATCH --exclusive
#SBATCH --output=/results/veprs/logs/V1-pilot-%j.out
#
# vep-rs 0.4.0 pilot. Establishes, before any campaign is committed:
#   1. that a usable JSON cache can be produced from the Perl VEP 115 cache;
#   2. what annotation configuration is actually reachable with that cache;
#   3. start-up and per-variant cost, measured at two input sizes rather than
#      extrapolated from one slice.
#
# Nothing here is a timed direct-vs-cached comparison. The point is to place
# vep-rs on the per-variant cost axis so the break-even rule can be tested.
set -uo pipefail
D=/results/veprs
DIST=$D/dist
PERLCACHE=/mnt/data/apps/ensembl-vep/115/cachedir/homo_sapiens
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
SIF=$D/dist/ensembl-vep_116.2.sif
JSONCACHE=$D/json_cache/GRCh38_chr22
W=/mnt/data/tmp/veprs-$SLURM_JOB_ID; mkdir -p "$W"; export TMPDIR="$W"
trap 'rm -rf "$W"' EXIT

say () { echo "[$(date -Is)] $*"; }
say "vep-rs $($DIST/vep --version 2>&1 | head -1) on $(hostname), $(nproc) CPUs"

# --- 1. container -----------------------------------------------------------
if [[ ! -s "$SIF" ]]; then
  say "pulling ensembl-vep container"
  apptainer pull --force "$SIF" docker://ensemblorg/ensembl-vep:release_116.2 \
    || { echo "CONTAINER PULL FAILED"; exit 1; }
fi
ls -lh "$SIF"

# --- 2. convert chr22 of the Perl cache -------------------------------------
# chr22 only: enough records to time against, a fraction of the conversion cost.
if [[ ! -s "$JSONCACHE/info.json" ]]; then
  rm -rf "$JSONCACHE"; mkdir -p "$JSONCACHE"
  say "converting Perl cache chr22 -> JSON"
  s=$(date +%s)
  apptainer exec -B "$PERLCACHE/115_GRCh38:/input:ro" -B "$JSONCACHE:/output" \
    -B "$DIST:/scripts:ro" "$SIF" \
    perl /scripts/storable_to_json.pl /input /output 22 \
    || { echo "CONVERSION FAILED"; exit 1; }
  say "conversion took $(( $(date +%s) - s ))s"
fi
[[ -s "$JSONCACHE/info.json" ]] || { echo "NO info.json AFTER CONVERSION"; exit 1; }
say "json cache: $(du -sh "$JSONCACHE" | cut -f1), $(find "$JSONCACHE" -name '*.json' | wc -l) files"
cat "$JSONCACHE/info.json" | head -20

# --- 3. two input sizes -----------------------------------------------------
# The same chr22 blueprint slices the SnpEff clinical gate used, so the inputs
# are directly comparable across tools.
bcftools view -r chr22 "$BP" -Ob -o "$W/big.bcf" -W --threads 8
bcftools view -r chr22:16000000-17500000 "$BP" -Ob -o "$W/small.bcf" -W --threads 8
for f in small big; do
  n=$(bcftools index -s "$W/$f.bcf" | awk '{s+=$3} END {print s+0}')
  [[ "$n" =~ ^[0-9]+$ && "$n" -gt 0 ]] || { echo "BAD RECORD COUNT for $f"; exit 1; }
  echo "$n" > "$W/$f.n"
  bcftools view "$W/$f.bcf" -Ov -o "$W/$f.vcf"
  say "$f slice: $n records"
done

# --- 4. time each configuration at both sizes -------------------------------
# --fork 8 matches the VEP arm. Levels probe how much the converted cache
# actually supports; a configuration that errors is recorded, not hidden.
TSV=$D/results/pilot_${SLURM_JOB_ID}.tsv
printf 'config\tsize\trecords\twall_seconds\trc\toutput_lines\tcsq_fields\n' > "$TSV"

run_cfg () {  # run_cfg <name> <extra vep flags...>
  local name="$1"; shift
  for size in small big; do
    local n out; n=$(cat "$W/$size.n"); out="$W/${name}_${size}.txt"
    local s e rc
    s=$(date +%s.%N)
    "$DIST/vep" --json_cache "$JSONCACHE" --assembly GRCh38 --fork 8 \
      --force_overwrite -i "$W/$size.vcf" -o "$out" "$@" \
      > "$W/${name}_${size}.stdout" 2> "$W/${name}_${size}.stderr"
    rc=$?; e=$(date +%s.%N)
    local wall lines fields
    wall=$(echo "$e - $s" | bc)
    lines=$(grep -vc '^#' "$out" 2>/dev/null || echo 0)
    fields=$(grep -m1 '^#Uploaded_variation' "$out" 2>/dev/null | awk -F'\t' '{print NF}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$name" "$size" "$n" "$wall" "$rc" "$lines" "${fields:-0}" >> "$TSV"
    say "$name/$size rc=$rc wall=${wall}s lines=$lines"
    [[ $rc -ne 0 ]] && tail -3 "$W/${name}_${size}.stderr"
  done
}

run_cfg core
run_cfg symbol_canonical --symbol --canonical
run_cfg hgvs --symbol --canonical --hgvs
run_cfg everything --everything

# Keep the tools' own output, not just numbers derived from it.
mkdir -p "$D/results/raw_${SLURM_JOB_ID}"
cp "$W"/*.stderr "$D/results/raw_${SLURM_JOB_ID}/" 2>/dev/null
head -40 "$W/core_small.txt" > "$D/results/raw_${SLURM_JOB_ID}/core_small.head.txt" 2>/dev/null
head -40 "$W/everything_small.txt" > "$D/results/raw_${SLURM_JOB_ID}/everything_small.head.txt" 2>/dev/null

say "--- results ---"; cat "$TSV"
say "done"
