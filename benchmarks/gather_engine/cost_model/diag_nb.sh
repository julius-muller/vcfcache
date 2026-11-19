#!/bin/bash
# How narrow should the miss-step sources be? chr1 of one sample, fastVEP alone (16 threads):
# all variants at fastVEP's advised 65 kb blocks (the uncached reference), and the misses of the
# 0.1% and the cohort cache at 65 kb, 4 kb, 1 kb and 256 bp blocks; for the 11- and 15-source stack.
# Needs the node-local caches of wgs_mission.sh STACK2=1 (mission-wgs2f, here KPGP-00128 on sl-w1)
# and the chr1 builds with NB=10 and NB=8 (wgs_mission_stack.sh, wgs_mission_stack2.sh, build_tracks1.sh).
# Usage: diag_nb.sh [sample] [chromosome]
set -uo pipefail
SAMPLE=${1:-KPGP-00128}; C=${2:-1}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
R=/mnt/data/tmp/mission-wgs2f; F=/results/fastvep05
OUT=$F/cost_model/diag_nb; ls $F $F/mission $F/cost_model > /dev/null 2>&1; mkdir -p $OUT
W=$(mktemp -d /mnt/data/tmp/diagnb-XXXX); trap 'rm -rf $W' EXIT; cd $W
log() { echo "[$(date -Is)] $*"; }
# ---- inputs: this chromosome, and its misses against each cache (exact CHROM POS REF ALT)
bcftools view -r chr$C $R/$SAMPLE/in.bcf -Ob -W -o in.bcf || exit 1
for k in af001 cohort; do
  bcftools view -r chr$C $R/$SAMPLE/c_$k.bcf -G -Ob -W -o c.bcf
  bcftools isec -C -w1 -c none in.bcf c.bcf -Ov -o miss_$k.vcf 2>/dev/null || exit 1
done
bcftools view in.bcf -Ov -o all.vcf
echo "chr$C of $SAMPLE: $(grep -vc '^#' all.vcf) variants; misses af001 $(grep -vc '^#' miss_af001.vcf), cohort $(grep -vc '^#' miss_cohort.vcf)" | tee $OUT/info.txt
# ---- sources, local copies: this chromosome's files plus the genome-wide ones; no per-chromosome tracks
new='^(onekg|gerp|dbnsfp)_|^encode_ccre'
copy() {  # dest dir...
  local d=$1; shift; mkdir -p $d
  for s in "$@"; do ls $s > /dev/null 2>&1
    for f in $s/*; do b=$(basename $f)
      case $b in tracka*|trackb*|gnomad_genomes_*) continue;; esac
      [[ $b == *_chr* && $b != *_chr${C}.* ]] && continue
      cp $f $d/ || exit 1
    done
  done
}
copy s16 $R/sa16
copy s12 $R/sa12
copy s10 $F/mission10/sa12 $F/mission210/sa12
copy s8 $F/mission8/sa12 $F/mission28/sa12
for b in 16 12 10 8; do cp $F/mission/tracks1/sa$b/track[ab].* s$b/ || exit 1
  mkdir -p e$b; for f in s$b/*; do [[ $(basename $f) =~ $new ]] || ln -s $W/$f e$b/; done   # the 11-source stack
done
for b in 16 12 10 8; do echo "s$b: $(ls s$b | wc -l) files, $(du -sh s$b | cut -f1); e$b: $(ls e$b | wc -l) files"; done | tee -a $OUT/info.txt
run() {  # label input sa-dir
  /usr/bin/time -f "%e %U" -o t.txt env RAYON_NUM_THREADS=16 $FV annotate --input $2 --output o.vcf --output-format vcf \
    --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $3 2> err.txt || { echo "FAILED $1 $3"; tail -3 err.txt; return; }
  printf "%s\t%s\t%s\t%s\n" $1 $3 "$(tail -1 t.txt)" $(grep -v '^#' o.vcf | md5sum | cut -c1-12) | tee -a $OUT/times.tsv
}
printf "input\tsa\twall user\tmd5\n" > $OUT/times.tsv
run warmup all.vcf s16
for rep in 1 2; do
  for st in s e; do
    run all all.vcf ${st}16
    for k in af001 cohort; do for b in 16 12 10 8; do run $k miss_$k.vcf $st$b; done; done
  done
done
echo DONE
