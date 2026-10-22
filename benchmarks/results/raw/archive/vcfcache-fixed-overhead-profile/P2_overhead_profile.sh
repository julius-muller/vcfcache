#!/bin/bash
#SBATCH --job-name=P2-overhead
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=8:00:00
#SBATCH --output=/results/veprs/logs/P2-overhead-%j.out
#
# Profile VCFcache's FIXED overhead at real panel/WES scale, legacy vs fused.
#
# Why: the campaign shows caching losing at panel (0.41-0.46x) on ~250 variants,
# direct 10.5 s against cached 22.4 s. At that size per-variant cost is
# irrelevant; the loss is the N-independent term. The fused engine halved the
# per-variant slope and left start-up flat, so it cannot help there.
#
# A head-node pre-flight already showed what the term is made of. 250 variants,
# zero misses, fused: 76.9 s total, of which the `route` stage was 76.5 s. Route
# streams the cache overlay, so the dominant fixed cost looks like a CACHE SCAN
# -- a function of how much of the 18.87 M-record cache must be touched, not of
# how many variants the sample has. That pre-flight input was also accidentally
# single-contig (the first 250 cache records, all chr1), which additionally
# disables the fused engine's contig parallelism.
#
# This job therefore separates three things that the pre-flight confounded:
#
#   ARM A  size at fixed genome-wide span. Inputs are stride-sampled across the
#          WHOLE cache, so 1 / 250 / 2000 / 75000 variants all span every contig,
#          as a real panel or exome does. Drawn FROM the cache, so every variant
#          is a hit, the recipe is never invoked, and the measurement is
#          VCFcache's own cost alone: T = s_vc + c_vc * N. N=1 observes s_vc
#          instead of extrapolating to it, which is the point.
#
#   ARM C  span at fixed size. 250 variants confined to 1, 4 and 22 contigs. If
#          the fixed term is a cache scan, this moves it while N is constant --
#          the decisive test, and the one that says whether region-restricted
#          cache access is the lever worth building.
#
#   ARM B  the real assays. Actual ACMG SF v3.3 panel and Twist core exome
#          samples, which have misses and so do invoke VEP. Ties the isolated
#          numbers to the campaign's observed loss and splits cached wall time
#          into VCFcache overhead versus annotator start-up. Guarded: a failure
#          is recorded and skipped, never fatal, because Arms A and C are the
#          deliverable.
#
# Every cell gets one UNTIMED warm-up then three timed reps. Cold state
# dominates a sub-minute cell, and ignoring that invalidated a whole Panel/WES
# campaign once already (GUIDELINES.md).
#
# `vcfcache --version` prints 0.5.3 here: that is the stale dist-info of the
# editable install this venv was built from. The CODE is feature/fused-routing
# at c878d86c3, shadowed by PYTHONPATH and asserted below via --engine.

set -euo pipefail

DEP=/results/deployment/vcfcache-fused-c878d86c
VENV=/results/deployment/vcfcache-recovery-95a01bf3e305/.venv
ENGINE=$DEP/native/target/release/vcfcache-engine
CACHE=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything
CBCF=$CACHE/vcfcache_annotated.bcf
PANEL=/mnt/data/vcfcache_benchmarks/samples/GRCh38/panel_acmg_sf_v3.3
WES=/mnt/data/vcfcache_benchmarks/samples/GRCh38/wes_twist_core
J=${SLURM_JOB_ID:-manual}
OUT=/results/veprs/results/overhead_$J
RAW=$OUT/raw
W=/mnt/data/tmp/overhead-$J
mkdir -p "$OUT" "$RAW" "$W/in" "$W/run"
chmod -R 777 "$W"
cd "$W"          # vcfcache resolves .env from the CWD, which must be readable

log(){ echo "[$(date -Is)] $*"; }
for f in "$ENGINE" "$VENV/bin/vcfcache" "$CBCF" "$PANEL" "$WES"; do
  [ -e "$f" ] || { echo "FATAL: missing $f" >&2; exit 1; }
done
VC=( env HOME=/home/appuser PYTHONPATH="$DEP" VCFCACHE_ENGINE="$ENGINE" "$VENV/bin/vcfcache" )
"${VC[@]}" annotate --help 2>&1 | grep -q -- "--engine" \
  || { echo "FATAL: CLI lacks --engine; the PYTHONPATH shadow is not in effect" >&2; exit 1; }
"${VC[@]}" annotate --help > "$RAW/cli_help.txt" 2>&1
"$ENGINE" --version > "$RAW/engine_version.txt" 2>&1
NCACHE=$(bcftools index -n "$CBCF")
log "engine: $(cat "$RAW/engine_version.txt"); cache: $NCACHE records"
bcftools index -s "$CBCF" > "$RAW/cache.index_stats.txt"

#### inputs ####
# One pass over the cache body, stride-sampled so the 75000-record sample spans
# every contig. Smaller sizes are nested subsets of it, so size varies while
# span stays genome-wide. Stripping QUAL/FILTER/INFO makes the rows look like a
# sample VCF rather than an already-annotated one.
HDR=$W/hdr_min.txt
{ bcftools view -h "$CBCF" | grep -E '^##(fileformat|contig|reference)'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'; } > "$HDR"
# A genome-wide pool, built from one indexed window per contig rather than by
# scanning the whole 18.87 M-record cache. 22 index seeks over ~110 Mb instead
# of a 3.1 Gb pass, and the result still spans every autosome.
#
# No `exit` in these awk programs. An awk that exits early SIGPIPEs bcftools,
# `pipefail` then reports the pipeline as failed, and `set -e` kills the script
# -- which is exactly how job 1252 died, silently, after logging this very step.
# GUIDELINES.md warns about this twice. The windows are small, so draining them
# is cheap and the guard is unnecessary; `set +o pipefail` is belt and braces.
: > "$W/pool.txt"
PER=$(( 75000 / 22 + 1 ))
set +o pipefail
for i in $(seq 1 22); do
  bcftools view -H "$CBCF" -r "chr$i:20000000-25000000" 2>/dev/null \
    | awk -v n="$PER" 'BEGIN{OFS="\t"} c<n {$6=".";$7=".";$8=".";NF=8; print; c++}' \
    >> "$W/pool.txt"
done
set -o pipefail
POOLN=$(wc -l < "$W/pool.txt")
POOLK=$(cut -f1 "$W/pool.txt" | sort -u | wc -l)
[ "$POOLN" -ge 75000 ] || { echo "FATAL: pool has only $POOLN rows, wanted >=75000" >&2; exit 1; }
[ "$POOLK" -eq 22 ] || { echo "FATAL: pool spans $POOLK contigs, wanted 22" >&2; exit 1; }
log "pool: $POOLN rows across $POOLK contigs" 

mk_bcf(){   # name  body-file
  local name=$1 body=$2
  cat "$HDR" "$body" > "$W/in/$name.vcf"
  bcftools view "$W/in/$name.vcf" -Ob -o "$W/in/$name.bcf" -W 2>/dev/null
  local n=$(bcftools index -n "$W/in/$name.bcf")
  local k=$(cut -f1 "$body" | sort -u | wc -l)
  echo -e "$name\t$n\t$k" >> "$RAW/input_inventory.tsv"
  log "input $name: $n records across $k contigs"
}
printf 'input\trecords\tcontigs\n' > "$RAW/input_inventory.tsv"
POOL=$(wc -l < "$W/pool.txt")
for N in 1 250 2000 75000; do
  awk -v n="$N" -v p="$POOL" 'BEGIN{step=int(p/n); if(step<1)step=1} NR%step==0 && c<n {print; c++}' \
    "$W/pool.txt" > "$W/body.txt"
  mk_bcf "gw$N" "$W/body.txt"
done
# Arm C: 250 records at three different contig spans.
for K in 1 4 22; do
  awk -v k="$K" 'BEGIN{for(i=1;i<=k;i++) want["chr"i]=1} $1 in want' "$W/pool.txt" \
    | awk 'c<250 {print; c++}' > "$W/body.txt"
  mk_bcf "span${K}c" "$W/body.txt"
done

TS=$OUT/overhead_$J.tsv
printf 'arm\tcell\tengine\trecords\tcontigs\trep\twall_seconds\trc\n' > "$TS"
printf 'arm\tcell\tengine\tstage\twall_seconds\n' > "$RAW/stage_timings.tsv"

ncontig_of(){ awk -F'\t' -v n="$1" '$1==n{print $3}' "$RAW/input_inventory.tsv"; }

run_cell(){   # arm cell engine input nrec ncontig reps fatal
  local arm=$1 cell=$2 eng=$3 inp=$4 nrec=$5 nk=$6 reps=$7 fatal=$8
  local sd rc st en wall i tag=${arm}_${cell}_${eng}
  sd=$W/run/warm_$tag; rm -rf "$sd"; mkdir -p "$sd"; chmod 777 "$sd"
  set +e
  "${VC[@]}" annotate -a "$CACHE" -i "$inp" -o "$sd/out.bcf" --engine "$eng" --stats-dir "$sd" \
      > "$sd/run.log" 2>&1
  rc=$?
  set -e
  if [ "$rc" -ne 0 ]; then
    log "$tag warm-up rc=$rc"
    cp "$sd/run.log" "$RAW/${tag}.warmup_fail.log" 2>/dev/null || true
    [ "$fatal" -eq 1 ] && { echo "FATAL: $tag failed its warm-up" >&2; tail -30 "$sd/run.log" >&2; exit 1; }
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$arm" "$cell" "$eng" "$nrec" "$nk" 0 "" "$rc" >> "$TS"
    return 0
  fi
  for i in $(seq 1 "$reps"); do
    sd=$W/run/${tag}_$i; rm -rf "$sd"; mkdir -p "$sd"; chmod 777 "$sd"
    st=$(date +%s.%N)
    set +e
    "${VC[@]}" annotate -a "$CACHE" -i "$inp" -o "$sd/out.bcf" --engine "$eng" --stats-dir "$sd" \
        > "$sd/run.log" 2>&1
    rc=$?
    set -e
    en=$(date +%s.%N)
    wall=$(awk -v a="$st" -v b="$en" 'BEGIN{printf "%.6f", b-a}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$arm" "$cell" "$eng" "$nrec" "$nk" "$i" "$wall" "$rc" >> "$TS"
    log "$tag rep$i rc=$rc wall=${wall}s"
    # Stage lines the workflow log emits, e.g. "Fused route completed in 76.497s".
    if [ "$i" -eq 1 ]; then
      grep -hoE "[A-Za-z-]+ (stage )?completed in [0-9.]+s" "$sd"/*_vcstats/workflow.log "$sd/run.log" 2>/dev/null \
        | sed -E "s/^/${arm}\t${cell}\t${eng}\t/; s/ completed in /\t/; s/s$//" >> "$RAW/stage_timings.tsv" || true
      mkdir -p "$RAW/stats/$tag"
      cp -r "$sd"/*_vcstats "$sd/run.log" "$RAW/stats/$tag/" 2>/dev/null || true
    fi
    [ "$rc" -eq 0 ] || { [ "$fatal" -eq 1 ] && { echo "FATAL: $tag rc=$rc" >&2; tail -30 "$sd/run.log" >&2; exit 1; }; }
  done
}

log "=== ARM A: size at genome-wide span, zero misses ==="
for N in 1 250 2000 75000; do
  for eng in legacy fused; do run_cell A "gw$N" "$eng" "$W/in/gw$N.bcf" "$N" "$(ncontig_of gw$N)" 3 1; done
done

log "=== ARM C: fixed size (250), varying contig span ==="
for K in 1 4 22; do
  for eng in legacy fused; do run_cell C "span${K}c" "$eng" "$W/in/span${K}c.bcf" 250 "$K" 3 1; done
done

log "=== ARM B: real panel and WES samples (VEP runs on their misses) ==="
for f in $(ls "$PANEL"/*.vcf.gz 2>/dev/null | head -2); do
  s=$(basename "$f" .vcf.gz); n=$(bcftools index -n "$f" 2>/dev/null || echo 0)
  for eng in legacy fused; do run_cell B "panel_$s" "$eng" "$f" "$n" "" 2 0; done
done
for f in $(ls "$WES"/*.vcf.gz 2>/dev/null | head -2); do
  s=$(basename "$f" .vcf.gz); n=$(bcftools index -n "$f" 2>/dev/null || echo 0)
  for eng in legacy fused; do run_cell B "wes_$s" "$eng" "$f" "$n" "" 2 0; done
done

cp "$TS" "$RAW/" 2>/dev/null || true
log "--- stage timings ---"; cat "$RAW/stage_timings.tsv"
log "--- results ---";       cat "$TS"
log "done; outputs in $OUT"
