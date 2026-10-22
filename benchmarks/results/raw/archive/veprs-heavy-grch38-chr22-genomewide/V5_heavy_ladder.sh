#!/bin/bash
#SBATCH --job-name=V5-heavy
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=6:00:00
#SBATCH --output=/results/veprs/logs/V5-heavy-%j.out
#
# vep-rs 0.4.0 in the annotation-heavy regime, over GENOME-WIDE gnomAD.
#
# Job 1242 put vep-rs at 5.68 us/variant on transcripts alone, far under the
# 29.7 us/variant caching needs to pay for itself. That is the tool's cheapest
# state. This job measures the other end: per-variant cost when vep-rs must
# traverse supplementary databases, and how many sources a recipe needs before
# caching wins.
#
# Nothing is downloaded. Sources already on the cluster:
#   gnomad38.vcf.gz / gnomad_v4.1_GRCh38_joint_master_af0001.bcf
#       99,942,855 records, GENOME-WIDE joint genomes+exomes, AF>=0.01%,
#       AF populated on every record. chr22: 1,626,030. The primary AF source
#       and the one the blueprints derive from.
#   gnomad_exomes_GRCh38.vcf.gz
#       183,717,261 records but exome-restricted, AF absent off-target.
#       Kept only for the exome-vs-genome-wide contrast.
#   clinvar_GRCh38.vcf.gz   4,465,932 records.
#
# Supersedes two failed attempts, both recorded rather than hidden:
#   1243 died building its first TSV -- the container was bound -B /mnt/data
#        only while the TSV path was under /results, so bgzip could not see it.
#   1245 died on the opposite mistake: moving everything to
#        /mnt/data/slurm-results, which is the HEAD node's name for the share.
#        On a worker that path exists but is an empty mountpoint. /results is
#        the real one, and both paths are bound below.
#   1247 died because the genome-wide gnomAD is under /mnt/data/vcfcache_data/raw,
#        which is head-node only. It is now staged into the shared tree.
#   V4   was never submitted: pre-flight found that vep-rs's dbNSFP plugin
#        SILENTLY DROPS any requested column named `AF` or ending in `_AF`.
#        Proven with identical values under four names: X_FREQ emitted on 60752
#        rows, X_AF / AF / gnomAD_joint_AF emitted on none, no warning. Hence
#        the frequency column below is GJ_FREQ.
#
# Division of labour between the two plugin shapes, which is not arbitrary:
#   gnomADc <whole VCF> annotates EVERY variant and is the correct vehicle for
#       variant-level allele frequency. It also does the realistic heavy I/O:
#       random access into the complete 1.1 GB / 99.9 M record genome-wide file.
#   dbNSFP <tabix TSV> annotates only TRANSCRIPT-OVERLAPPING consequences, as
#       the real dbNSFP does. Right for transcript-centric scores, and it is why
#       a dbNSFP-delivered field covers a fraction of a WGS input.
#
# The W* sweep points several plugin instances at the SAME whole file. The
# payload is then redundant, so the sweep is a COST model of a multi-source
# recipe and not an annotation anyone would run -- only the per-instance I/O and
# decode work is real, which is the single quantity it claims. It is the I/O
# analogue of benchmarks/vep_plugins/SyntheticDelay.pm.

set -euo pipefail

# Worker-side paths. /results is the real network mount on the compute nodes;
# /mnt/data/slurm-results exists there too but is a different, empty mountpoint
# -- that mix-up is what killed job 1245. Both /results and /mnt/data must be
# bound into the container, and both exist on the workers (verified by srun).
D=/results/veprs
SIF=$D/dist/ensembl-vep_116.2.sif
V=$D/dist/vep
JC=$D/json_cache/GRCh38_chr22
SRC=/results/clinical_sources
CV=$SRC/clinvar_GRCh38.vcf.gz
GE=$SRC/gnomad_exomes_GRCh38.vcf.gz
# The genome-wide joint gnomAD lives under /mnt/data/vcfcache_data/raw on the
# HEAD node only -- a worker's /mnt/data/vcfcache_data holds just `caches`, which
# is what killed job 1247. It is therefore staged into the shared results tree,
# with its .tbi built once outside this job (vep-rs refuses the shipped .csi).
GW=/results/veprs/data/gnomad_joint_gw_GRCh38.vcf.gz
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
J=${SLURM_JOB_ID:-manual}
OUT=$D/results/heavy_$J
RAW=$OUT/raw
TSV=$OUT/tsv
W=/mnt/data/tmp/veprs-heavy-$J
mkdir -p "$OUT" "$RAW" "$TSV" "$W"
APP=(apptainer exec -B /mnt/data -B /results "$SIF")

log(){ echo "[$(date -Is)] $*"; }
log "vep-rs $("$V" --version 2>&1 | tail -1) on $(hostname), $(nproc) CPUs"
for f in "$SIF" "$V" "$JC" "$CV" "$GE" "$GW" "$GW.tbi" "$BP"; do
  [ -e "$f" ] || { echo "FATAL: missing $f" >&2; exit 1; }
done
"${APP[@]}" bash -c 'command -v bgzip tabix' > "$RAW/container_htslib.txt" \
  || { echo "FATAL: bgzip/tabix not reachable in the container" >&2; exit 1; }

#### 1. input slices -- the same chr22 cuts job 1242 and the SnpEff gate used ####

bcftools view -r chr22                   "$BP" -Ob -o "$W/big.bcf"   -W --threads 8
bcftools view -r chr22:16000000-17500000 "$BP" -Ob -o "$W/small.bcf" -W --threads 8
declare -A NREC
for s in small big; do
  bcftools index -s "$W/$s.bcf" > "$RAW/$s.index_stats.txt"
  n=$(awk '{s+=$3} END {print s+0}' "$RAW/$s.index_stats.txt")
  case "$n" in ''|*[!0-9]*) echo "FATAL: $s record count not numeric: '$n'" >&2; exit 1;; esac
  [ "$n" -gt 0 ] || { echo "FATAL: $s is empty" >&2; exit 1; }
  NREC[$s]=$n
  bcftools view "$W/$s.bcf" -Ov -o "$W/$s.vcf"
  log "$s slice: $n records"
done

#### 2. genome-wide source check ####

log "genome-wide joint gnomAD: $(du -h "$GW" | cut -f1), tbi $(du -h "$GW.tbi" | cut -f1)"
nrec=$(bcftools index -n "$GW")
[ "$nrec" -gt 90000000 ] || { echo "FATAL: genome-wide source has only $nrec records" >&2; exit 1; }
log "genome-wide joint gnomAD: $nrec records"

#### 3. dbNSFP-format supplementary TSVs, built from what we hold ####

build_tsv(){            # name  source  "bcftools format"  "col names"
  # Two statements, not one: `local` expands all its words before it runs, so
  # $name inside the same declaration is still unset under `set -u`.
  local name=$1 src=$2 fmt=$3 cols=$4
  local out=$TSV/$name.tsv
  # %b, not %s: %s would leave the \t separators in $cols as literal text.
  printf '#chr\tpos(1-based)\tref\talt\t%b\n' "$cols" > "$out"
  bcftools query -r chr22 -f "$fmt" "$src" >> "$out"
  local n=$(( $(wc -l < "$out") - 1 ))
  [ "$n" -gt 0 ] || { echo "FATAL: $name TSV is empty" >&2; exit 1; }
  "${APP[@]}" bgzip -f "$out"
  "${APP[@]}" tabix -s1 -b2 -e2 -f "$out.gz"
  [ -s "$out.gz.tbi" ] || { echo "FATAL: no tabix index for $out.gz" >&2; exit 1; }
  log "built $name.tsv.gz: $n rows, $(du -h "$out.gz" | cut -f1)"
}

# GJ_FREQ, not GJ_AF: a column named *_AF is dropped without warning.
build_tsv gnomad_joint "$GW" '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n' 'GJ_FREQ\tGJ_AC'
build_tsv clinvar      "$CV" '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\n' 'CV_CLNSIG\tCV_REVSTAT\tCV_CLNDN'
ls -la "$TSV" > "$RAW/tsv_listing.txt"

GA="dbNSFP,$TSV/gnomad_joint.tsv.gz,GJ_FREQ,GJ_AC"
CA="dbNSFP,$TSV/clinvar.tsv.gz,CV_CLNSIG,CV_REVSTAT,CV_CLNDN"
HJ="gnomADc,$GW"    # whole-file random access, genome-wide, 99.9 M records
HE="gnomADc,$GE"    # whole-file random access, exome-restricted, 183.7 M records

#### 4. the recipe ladder ####

declare -a LADDER=(
  "L0_core|"
  "L1_gw_af|--plugin $HJ"
  "L2_clinical_gw|--plugin $HJ --plugin $CA"
  "L3_clinical_hgvs|--hgvs --symbol --canonical --plugin $HJ --plugin $CA"
  "L4_clinical_plus_scores|--hgvs --symbol --canonical --plugin $HJ --plugin $CA --plugin $GA"
  "W2_gw_wholefile_x2|--plugin $HJ --plugin $HJ"
  "W4_gw_wholefile_x4|--plugin $HJ --plugin $HJ --plugin $HJ --plugin $HJ"
  "W8_gw_wholefile_x8|--plugin $HJ --plugin $HJ --plugin $HJ --plugin $HJ --plugin $HJ --plugin $HJ --plugin $HJ --plugin $HJ"
  "X1_exome_wholefile|--plugin $HE"
  "E1_everything_clinical_gw|--everything --plugin $HJ --plugin $CA"
  "E2_everything_heavy_gw|--everything --plugin $HJ --plugin $CA --plugin $GA --plugin $HJ"
)

TS=$OUT/heavy_$J.tsv
printf 'config\tsize\trecords\trep\twall_seconds\trc\toutput_lines\n' > "$TS"
printf 'config\tsize\tcolumn\trows_with_value\n' > "$RAW/field_assertions.tsv"

for entry in "${LADDER[@]}"; do
  cfg=${entry%%|*}; args=${entry#*|}
  for size in small big; do
    n=${NREC[$size]}
    # rep 1 is first touch for this config. Page-cache state is part of the
    # measurement and is reported per rep, never averaged away.
    for rep in 1 2 3; do
      o=$W/${cfg}_${size}.txt
      st=$(date +%s.%N)
      set +e
      # shellcheck disable=SC2086
      "$V" --assembly GRCh38 --force_overwrite --fork 8 --json_cache "$JC" \
           --input_file "$W/$size.vcf" --output_file "$o" $args \
           > "$W/stdout.txt" 2> "$W/stderr.txt"
      rc=$?
      set -e
      en=$(date +%s.%N)
      wall=$(awk -v a="$st" -v b="$en" 'BEGIN{printf "%.6f", b-a}')
      lines=$(grep -vc '^#' "$o" 2>/dev/null || echo 0)
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$cfg" "$size" "$n" "$rep" "$wall" "$rc" "$lines" >> "$TS"
      [ "$rep" = 1 ] && cp "$W/stderr.txt" "$RAW/${cfg}_${size}.stderr"
      log "$cfg/$size rep$rep rc=$rc wall=${wall}s lines=$lines"
      # Before any artifact capture: a failed run has no output file, and reaching
      # for it would report a missing file instead of the error that caused it.
      [ "$rc" -eq 0 ] || { echo "FATAL: $cfg/$size returned $rc" >&2; sed -n '1,40p' "$W/stderr.txt" >&2; exit 1; }
      if [ "$rep" = 1 ] && [ "$size" = small ]; then
        head -70 "$o" > "$RAW/${cfg}_small.head.txt"
      fi
    done
    # Assert that every field the level is defined by is actually populated.
    # This is the check that caught the silent *_AF drop; it is not ceremony.
    case "$cfg" in
      L0_*)                     want= ;;
      L1_*|W*|X1_*)             want=gnomADc_INFO ;;
      L2_*|L3_*|E1_*)           want=dbNSFP_CV_CLNSIG ;;
      L4_*|E2_*)                want=dbNSFP_GJ_FREQ ;;
      *)                        want= ;;
    esac
    if [ -n "$want" ]; then
      set +o pipefail   # a matching grep would SIGPIPE the producer; see GUIDELINES.md
      hits=$(grep -v '^#' "$W/${cfg}_${size}.txt" | grep -c "$want" || true)
      set -o pipefail
      printf '%s\t%s\t%s\t%s\n' "$cfg" "$size" "$want" "$hits" >> "$RAW/field_assertions.tsv"
      [ "$hits" -gt 0 ] || { echo "FATAL: $cfg/$size claims $want but emitted none" >&2; exit 1; }
      log "$cfg/$size asserted $want on $hits rows"
    fi
  done
done

#### 5. how much of each supplementary column is actually populated ####
# Genome-wide joint versus exome-restricted is the point of this job, and
# dbNSFP's transcript-only reach is the other thing a reader needs. Both are
# recorded here rather than inferred later from a figure.
{
  printf 'config\tcolumn\trows_with_value\ttotal_rows\n'
  for cfg in L1_gw_af L2_clinical_gw L4_clinical_plus_scores X1_exome_wholefile E1_everything_clinical_gw; do
    f=$W/${cfg}_big.txt
    [ -s "$f" ] || continue
    set +o pipefail
    tot=$(grep -vc '^#' "$f" || true)
    for col in gnomADc_INFO dbNSFP_CV_CLNSIG dbNSFP_GJ_FREQ dbNSFP_GJ_AC; do
      c=$(grep -v '^#' "$f" | grep -c "$col=" || true)
      printf '%s\t%s\t%s\t%s\n' "$cfg" "$col" "$c" "$tot"
    done
    set -o pipefail
  done
} > "$RAW/supplementary_coverage.tsv"

# Keep one real annotated row per level, so the output can be inspected later
# without rerunning anything.
for cfg in L1_gw_af L2_clinical_gw L4_clinical_plus_scores E2_everything_heavy_gw; do
  f=$W/${cfg}_small.txt
  [ -s "$f" ] || continue
  set +o pipefail
  grep -v '^#' "$f" | grep -m3 'gnomADc_INFO=' > "$RAW/${cfg}.sample_rows.txt" || true
  set -o pipefail
done

cp "$TS" "$RAW/" 2>/dev/null || true
log "--- supplementary coverage ---"; cat "$RAW/supplementary_coverage.tsv"
log "--- field assertions ---";       cat "$RAW/field_assertions.tsv"
log "--- results ---";                cat "$TS"
log "done; outputs in $OUT"
