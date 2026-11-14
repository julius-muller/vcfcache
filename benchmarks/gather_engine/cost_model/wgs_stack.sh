#!/bin/bash
# Build the per-chromosome sources of fastVEP's ACMG stack for one chromosome
# (gnomAD exomes, and SpliceAI and PhyloP distilled from it, and REVEL), as
# acmg_probe/build.sh does for chr1. Usage: wgs_stack.sh <chromosome name
# without "chr">; run as one Slurm task per chromosome.
set -euo pipefail
C=$1
D=/results/fastvep05/acmg; S=$D/sources; OUT=$D/sa/wgs; mkdir -p $OUT
FV=/mnt/data/apps/fastvep/0.5.0/bin/fastvep
[ -e $OUT/chr$C.DONE ] && exit 0
W=$(mktemp -d /mnt/data/tmp/wgs-stack-$C-XXXX); trap 'rm -rf $W' EXIT
log() { echo "[$(date -Is)] chr$C $*"; }
b() { local s=$(date +%s); $FV sa-build --no-progress --assembly GRCh38 "$@"; log "sa-build $1 $2 $(( $(date +%s) - s ))s"; }
G=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr$C.vcf.bgz
if [ $C = 1 ]; then V=$S/gnomad.exomes.v4.1.sites.chr1.vcf.bgz; else
  V=$W/gnomad.vcf.bgz
  curl -fsSL --retry 10 --retry-delay 30 -o $V $G; curl -fsSL --retry 10 -o $V.tbi $G.tbi
  log "downloaded $(du -h $V | cut -f1)"
fi
b --source gnomad -i $V -o $W/gnomad_chr$C --chunk-bits 16
python3 $D/scripts/extract_gnomad_scores.py --chrom $C --url $V --out-dir $W
b --source spliceai -i $W/spliceai_chr$C.vcf -o $W/spliceai_chr$C --chunk-bits 16
b --source phylop -i $W/phylop_chr$C.tsv -o $W/phylop_chr$C --chunk-bits 16
python3 -c 'import shutil,sys,zipfile; shutil.copyfileobj(zipfile.ZipFile(sys.argv[1]).open("revel_with_transcript_ids"), sys.stdout.buffer, 1 << 24)' \
  $S/revel-v1.3_all_chromosomes.zip | awk -F, -v c=$C 'NR==1 || $1==c' > $W/revel_chr$C.csv
if [ $(wc -l < $W/revel_chr$C.csv) -gt 1 ]; then b --source revel -i $W/revel_chr$C.csv -o $W/revel_chr$C; fi
cp $W/*_chr$C.osa2 $OUT/
touch $OUT/chr$C.DONE; log "done"; ls -la $OUT | grep "_chr$C\."
