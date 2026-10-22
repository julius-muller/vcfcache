#!/bin/bash
#SBATCH --job-name=P3-span
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=6:00:00
#SBATCH --output=/results/veprs/logs/P3-span-%j.out
#
# The decisive test for where VCFcache's fixed overhead comes from:
# hold input size CONSTANT at 250 variants and vary only how many contigs they
# span (1, 2, 4, 8, 22). Inputs are drawn from the cache, so every variant is a
# hit, the annotation recipe is never invoked, and nothing but VCFcache's own
# cost is being timed.
#
# A head-node pre-flight timed 250 single-contig variants at 76.9 s, of which
# the `route` stage was 76.5 s. Route streams the cache overlay. If the fixed
# cost is a cache scan rather than a start-up constant, it must move with span
# while N stays fixed. That is what this measures, and it decides whether
# region-restricted cache access is the lever worth building.
#
# This arm exists separately because job 1253's version of it was wrong: the
# genome-wide pool is ordered chr1..chr22, so taking the first 250 matching rows
# always landed entirely in chr1 whatever the requested span. All three cells
# reported one contig. Here the rows are taken round-robin, PER CONTIG, and the
# realised span is asserted before anything is timed.

set -euo pipefail

DEP=/results/deployment/vcfcache-fused-c878d86c
VENV=/results/deployment/vcfcache-recovery-95a01bf3e305/.venv
ENGINE=$DEP/native/target/release/vcfcache-engine
CACHE=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything
CBCF=$CACHE/vcfcache_annotated.bcf
J=${SLURM_JOB_ID:-manual}
OUT=/results/veprs/results/span_$J
RAW=$OUT/raw
W=/mnt/data/tmp/span-$J
mkdir -p "$OUT" "$RAW" "$W/in" "$W/run"; chmod -R 777 "$W"
cd "$W"

log(){ echo "[$(date -Is)] $*"; }
for f in "$ENGINE" "$VENV/bin/vcfcache" "$CBCF"; do
  [ -e "$f" ] || { echo "FATAL: missing $f" >&2; exit 1; }
done
VC=( env HOME=/home/appuser PYTHONPATH="$DEP" VCFCACHE_ENGINE="$ENGINE" "$VENV/bin/vcfcache" )
"${VC[@]}" annotate --help 2>&1 | grep -q -- "--engine" \
  || { echo "FATAL: CLI lacks --engine; PYTHONPATH shadow not in effect" >&2; exit 1; }

{ bcftools view -h "$CBCF" | grep -E '^##(fileformat|contig|reference)'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'; } > "$W/hdr.txt"

N=250
printf 'input\trecords\tcontigs\n' > "$RAW/input_inventory.tsv"
# No awk `exit` anywhere below: an early exit SIGPIPEs bcftools and pipefail
# then fails the script. That killed job 1252. See GUIDELINES.md.
for K in 1 2 4 8 22; do
  per=$(( N / K + 1 ))
  : > "$W/body.txt"
  set +o pipefail
  for i in $(seq 1 "$K"); do
    bcftools view -H "$CBCF" -r "chr$i:20000000-25000000" 2>/dev/null \
      | awk -v n="$per" 'BEGIN{OFS="\t"} c<n {$6=".";$7=".";$8=".";NF=8; print; c++}' >> "$W/body.txt"
  done
  set -o pipefail
  # Trim to exactly N while keeping the spread: drop from the tail of each
  # contig's block rather than truncating whole contigs away.
  awk -v n="$N" -v k="$K" 'BEGIN{keep=int(n/k)} {c[$1]++; if(c[$1]<=keep) print}' "$W/body.txt" > "$W/trim.txt"
  cat "$W/hdr.txt" "$W/trim.txt" > "$W/in/span${K}c.vcf"
  bcftools view "$W/in/span${K}c.vcf" -Ob -o "$W/in/span${K}c.bcf" -W 2>/dev/null
  got=$(bcftools index -n "$W/in/span${K}c.bcf")
  gotk=$(cut -f1 "$W/trim.txt" | sort -u | wc -l)
  printf 'span%sc\t%s\t%s\n' "$K" "$got" "$gotk" >> "$RAW/input_inventory.tsv"
  [ "$gotk" -eq "$K" ] || { echo "FATAL: span${K}c spans $gotk contigs, wanted $K" >&2; exit 1; }
  [ "$got" -gt 0 ] || { echo "FATAL: span${K}c is empty" >&2; exit 1; }
  log "input span${K}c: $got records across $gotk contigs"
done
cat "$RAW/input_inventory.tsv"

TS=$OUT/span_$J.tsv
printf 'cell\tengine\trecords\tcontigs\trep\twall_seconds\trc\n' > "$TS"
printf 'cell\tengine\tstage\twall_seconds\n' > "$RAW/stage_timings.tsv"

for K in 1 2 4 8 22; do
  inp=$W/in/span${K}c.bcf
  nrec=$(bcftools index -n "$inp")
  for eng in legacy fused; do
    tag=span${K}c_$eng
    # One untimed warm-up: cold state would otherwise land entirely on rep 1.
    sd=$W/run/warm_$tag; rm -rf "$sd"; mkdir -p "$sd"; chmod 777 "$sd"
    "${VC[@]}" annotate -a "$CACHE" -i "$inp" -o "$sd/out.bcf" --engine "$eng" \
        --stats-dir "$sd" > "$sd/run.log" 2>&1 \
      || { echo "FATAL: $tag warm-up failed" >&2; tail -30 "$sd/run.log" >&2; exit 1; }
    for rep in 1 2 3; do
      sd=$W/run/${tag}_$rep; rm -rf "$sd"; mkdir -p "$sd"; chmod 777 "$sd"
      st=$(date +%s.%N)
      set +e
      "${VC[@]}" annotate -a "$CACHE" -i "$inp" -o "$sd/out.bcf" --engine "$eng" \
          --stats-dir "$sd" > "$sd/run.log" 2>&1
      rc=$?
      set -e
      en=$(date +%s.%N)
      wall=$(awk -v a="$st" -v b="$en" 'BEGIN{printf "%.6f", b-a}')
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "span${K}c" "$eng" "$nrec" "$K" "$rep" "$wall" "$rc" >> "$TS"
      log "$tag rep$rep rc=$rc wall=${wall}s"
      [ "$rc" -eq 0 ] || { echo "FATAL: $tag rc=$rc" >&2; tail -30 "$sd/run.log" >&2; exit 1; }
      if [ "$rep" -eq 1 ]; then
        grep -hoE "[A-Za-z-]+ (stage )?completed in [0-9.]+s" "$sd"/*_vcstats/workflow.log "$sd/run.log" 2>/dev/null \
          | sed -E "s/^/span${K}c\t${eng}\t/; s/ completed in /\t/; s/s$//" >> "$RAW/stage_timings.tsv" || true
        mkdir -p "$RAW/stats/$tag"; cp -r "$sd"/*_vcstats "$sd/run.log" "$RAW/stats/$tag/" 2>/dev/null || true
      fi
    done
  done
done

cp "$TS" "$RAW/" 2>/dev/null || true
log "--- stage timings ---"; cat "$RAW/stage_timings.tsv"
log "--- results ---";       cat "$TS"
log "done; outputs in $OUT"
