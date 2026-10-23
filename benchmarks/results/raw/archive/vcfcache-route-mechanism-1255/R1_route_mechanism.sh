#!/bin/bash
#SBATCH --job-name=R1-mech
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=6:00:00
#SBATCH --output=/results/veprs/logs/R1-mech-%j.out
#
# What drives the cost of the fused engine's `route` stage: the VOLUME of cache
# it must read, its POSITION in the file, or THREADING?
#
# Jobs 1253/1254 left this open. Per-worker route times at 22 contigs rose almost
# arithmetically (9.9 / 21.3 / 32.6 / 45.6 s) and three explanations fit equally
# well: grouping balances `bcftools index --stats` of the INPUT (variant counts)
# while cost is driven by cache records; repeated cache-prefix traversal; and the
# arms were not thread-controlled (legacy's overlay gets --threads 8, the
# bcftools annotate inside route gets none, and compression_threads was 8//4-2=0).
#
# The cache is coordinate-sorted, which makes volume and position separable:
#
#   contig   own cache records   cumulative (file position)
#   chr1          1,435,577            1,435,577   <- largest, at the START
#   chr21           303,065           17,845,933   <- small, near the END
#
# Volume predicts chr1 ~4.7x SLOWER than chr21. Prefix traversal predicts chr21
# ~12.4x slower than chr1. Opposite directions, so one measurement decides it.
# All 22 autosomes plus chrX are run, giving a regression rather than one pair.
#
# This times the route stage DIRECTLY, reproducing the exact command from
# fused.py rather than driving the CLI. That also sidesteps the workers==1 branch
# (which sets contigs=[None] and drops the -r restriction entirely), so region
# restriction is always active here -- the condition the open question is about.
#
# Reference arm: a pure indexed read of the same contig's cache
# (`bcftools view -H -r chrN cache | wc -l`). If route tracks THIS, cost is
# volume and the index seek works. If route instead tracks cumulative position
# while the reference does not, the traversal is inside annotate.

set -euo pipefail

DEP=/results/deployment/vcfcache-fused-c878d86c
NATIVE=$DEP/native/target/release/vcfcache-engine
CACHE=/results/../vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything
CACHE=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything
CBCF=$CACHE/vcfcache_annotated.bcf
TAG=CSQ
BCFTOOLS=$(command -v bcftools)
J=${SLURM_JOB_ID:-manual}
OUT=/results/veprs/results/mech_$J
RAW=$OUT/raw
W=/mnt/data/tmp/mech-$J
mkdir -p "$OUT" "$RAW" "$W/in" "$W/run"; chmod -R 777 "$W"
cd "$W"

log(){ echo "[$(date -Is)] $*"; }
for f in "$NATIVE" "$CBCF" "$BCFTOOLS"; do [ -e "$f" ] || { echo "FATAL: missing $f" >&2; exit 1; }; done
"$NATIVE" --version > "$RAW/engine_version.txt" 2>&1
log "engine: $(cat "$RAW/engine_version.txt"); bcftools: $("$BCFTOOLS" --version | head -1)"

CONTIGS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"

# Per-contig cache volume and cumulative file position, straight from the index.
"$BCFTOOLS" index -s "$CBCF" > "$RAW/cache.index_stats.txt"
awk -F'\t' 'BEGIN{OFS="\t"; print "contig","cache_records","cumulative_records"}
            {cum+=$3; print $1,$3,cum}' "$RAW/cache.index_stats.txt" > "$RAW/cache_geometry.tsv"
cat "$RAW/cache_geometry.tsv"

# Inputs: ~250 variants taken from each contig's own cache slice, so every
# variant is a hit and no annotator is ever involved. No awk `exit` anywhere --
# an early exit SIGPIPEs bcftools and pipefail then fails the script (killed job
# 1252). See GUIDELINES.md.
{ "$BCFTOOLS" view -h "$CBCF" | grep -E '^##(fileformat|contig|reference)'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'; } > "$W/hdr.txt"
set +o pipefail
for c in $CONTIGS; do
  "$BCFTOOLS" view -H "$CBCF" -r "$c" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"} c<250 {$6=".";$7=".";$8=".";NF=8; print; c++}' > "$W/body.txt"
  cat "$W/hdr.txt" "$W/body.txt" > "$W/in/$c.vcf"
  "$BCFTOOLS" view "$W/in/$c.vcf" -Ob0 -o "$W/in/${c}_norm.bcf" -W 2>/dev/null
  n=$("$BCFTOOLS" index -n "$W/in/${c}_norm.bcf")
  [ "$n" -gt 0 ] || { echo "FATAL: $c input empty" >&2; exit 1; }
done
set -o pipefail
log "inputs built for $(echo $CONTIGS | wc -w) contigs"

TS=$OUT/mech_$J.tsv
printf 'contig\tcache_records\tcumulative_records\tarm\tthreads\trep\twall_seconds\trc\thits\n' > "$TS"

geom(){
  local g=$(awk -F'\t' -v c="$1" '$1==c{print $2"\t"$3}' "$RAW/cache_geometry.tsv")
  if [ -z "$g" ]; then
    g=$(awk -F'\t' -v c="$1" '$1==c{print $3"\t"}' "$RAW/group_geometry.tsv")
    [ -n "$g" ] || g=$'\t'
  fi
  printf '%s' "$g"
}

time_it(){  # contig arm threads reps cmd...
  local c=$1 arm=$2 th=$3 reps=$4; shift 4
  local g=$(geom "$c") st en wall rc i hits sd
  for i in $(seq 0 "$reps"); do       # i=0 is an untimed warm-up
    sd=$W/run/${c}_${arm}_t${th}; rm -rf "$sd"; mkdir -p "$sd"; cd "$sd"
    st=$(date +%s.%N)
    set +e
    ( "$@" ) > "$sd/out.log" 2> "$sd/err.log"
    rc=$?
    set -e
    en=$(date +%s.%N)
    cd "$W"
    [ "$i" -eq 0 ] && continue
    wall=$(awk -v a="$st" -v b="$en" 'BEGIN{printf "%.4f", b-a}')
    hits=""
    [ -s "$sd/route.json" ] && hits=$(python3 -c "
import json; d=json.load(open('$sd/route.json'))
print('%s/%s' % (d.get('input_variants'), d.get('missing_variants')))" 2>/dev/null || true)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$c" "$g" "$arm" "$th" "$i" "$wall" "$rc" "$hits" >> "$TS"
    log "$c/$arm/threads=$th rep$i rc=$rc wall=${wall}s hits=$hits"
    [ "$rc" -eq 0 ] || { echo "FATAL: $c/$arm rc=$rc" >&2; tail -20 "$sd/err.log" >&2; exit 1; }
  done
}

# Reference: a pure indexed read of this contig's cache slice. Volume only.
ref_cmd(){ bash -c "$BCFTOOLS view -H -r $1 $CBCF | wc -l"; }

# The real route stage, verbatim from fused.py line 243-246, with -r always set.
route_cmd(){
  local c=$1 th=$2 sel="-r $c --regions-overlap 0" thr=""
  [ "$th" -gt 0 ] && thr="--threads $th"
  bash -c "$BCFTOOLS annotate $sel $thr -a $CBCF $W/in/${c}_norm.bcf -c INFO -Ou \
    | $NATIVE route --hits hits.bcf --misses misses.raw.bcf --order order.bin \
      --stats route.json --tag $TAG"
}

route_grp(){
  local gi=$1 g=$2
  bash -c "$BCFTOOLS annotate -r '$g' --regions-overlap 0 -a $CBCF $W/in/grp${gi}_norm.bcf -c INFO -Ou \
    | $NATIVE route --hits hits.bcf --misses misses.raw.bcf --order order.bin \
      --stats route.json --tag $TAG"
}
printf 'group\tcontigs\tcache_records\n' > "$RAW/group_geometry.tsv"

log "=== ARM ref: indexed read of each contig's cache slice ==="
for c in $CONTIGS; do time_it "$c" ref 0 2 ref_cmd "$c"; done

log "=== ARM route: the fused route stage, threads as shipped (none) ==="
for c in $CONTIGS; do time_it "$c" route 0 2 route_cmd "$c" 0; done

log "=== ARM route8: same, with --threads 8 on the annotate, to size threading ==="
for c in chr1 chr11 chr21 chr22; do time_it "$c" route8 8 2 route_cmd "$c" 8; done

# The four groups fused.py builds for a 22-contig input at 4 workers. Contigs
# arrive in LEXICOGRAPHIC order and are grouped by INPUT variant count, which for
# ~11 variants per contig means roughly equal contig counts and therefore very
# unequal cache volume. Run serially, so each group's own cost is visible without
# the concurrency contention that inflated job 1254's per-worker numbers.
log "=== ARM group: the real lexicographic worker groups, run serially ==="
G1=chr1,chr10,chr11,chr12,chr13
G2=chr14,chr15,chr16,chr17,chr18
G3=chr19,chr2,chr20,chr21,chr22
G4=chr3,chr4,chr5,chr6,chr7,chr8,chr9
gi=0
for g in "$G1" "$G2" "$G3" "$G4"; do
  gi=$((gi+1))
  set +o pipefail
  "$BCFTOOLS" view -H "$CBCF" -r "$g" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"} c<250 {$6=".";$7=".";$8=".";NF=8; print; c++}' > "$W/body.txt"
  set -o pipefail
  cat "$W/hdr.txt" "$W/body.txt" > "$W/in/grp$gi.vcf"
  "$BCFTOOLS" view "$W/in/grp$gi.vcf" -Ob0 -o "$W/in/grp${gi}_norm.bcf" -W 2>/dev/null
  vol=$(awk -F'\t' -v g="$g" 'BEGIN{n=split(g,a,",");for(i=1;i<=n;i++)want[a[i]]=1}
                               $1 in want {s+=$3} END{print s+0}' "$RAW/cache.index_stats.txt")
  log "group $gi = $g : $vol cache records"
  echo -e "grp$gi	$g	$vol" >> "$RAW/group_geometry.tsv"
  time_it "grp$gi" group 0 2 route_grp "$gi" "$g"
done

cp "$TS" "$RAW/" 2>/dev/null || true
log "--- results ---"; cat "$TS"
log "done; outputs in $OUT"
