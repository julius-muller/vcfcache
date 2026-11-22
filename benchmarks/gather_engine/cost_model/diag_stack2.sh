#!/bin/bash
# Which of the five further sources makes fastVEP ~11x slower? chr22 of one sample,
# fastVEP alone (16 threads), the mission stack plus each new source on its own.
# Needs the node-local source copies of wgs_mission.sh STACK2=1 (sl-w1/3/5).
# Usage: diag_stack2.sh [sample] [chromosome]
set -uo pipefail
SAMPLE=${1:-KPGP-00128}; C=${2:-22}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
R=/mnt/data/tmp/mission-wgs2; M16=$R/sa16
OUT=/results/fastvep05/cost_model/diag_stack2; ls /results/fastvep05/cost_model > /dev/null 2>&1; mkdir -p $OUT
W=$(mktemp -d /mnt/data/tmp/diag2-XXXX); trap 'rm -rf $W' EXIT; cd $W
bcftools view -r chr$C $R/$SAMPLE/in.bcf -Ov -o in.vcf || exit 1
echo "chr$C of $SAMPLE: $(grep -vc '^#' in.vcf) variants, node $(hostname)" | tee $OUT/info.txt
genome=$(ls $M16 | grep -v '_chr')                         # gene tables and BED indexes, genome-wide
base="alphamissense dbsnp gnomad phylop revel spliceai tracka trackb"
new="gnomad_genomes onekg gerp dbnsfp"
mk() {  # dir sources...: links to this chromosome's files of the given sources
  local d=$1; shift; mkdir -p $d
  for g in $genome; do case $g in encode_ccre*) ;; *) ln -s $M16/$g $d/;; esac; done
  for s in "$@"; do
    case $s in encode) ln -s $M16/encode_ccre.osi $d/;;
      *) for f in $M16/${s}_chr${C}.* $M16/${s}_chr${C}_p*; do [ -e $f ] && ln -s $f $d/; done;; esac
  done
}
mk base $base
for s in $new encode; do mk plus_$s $base $s; done
mk all $base $new encode
run() {  # config threads
  local s=$(date +%s.%N)
  /usr/bin/time -f "%e %U %M" -o t.txt env RAYON_NUM_THREADS=$2 $FV annotate --input in.vcf --output o.vcf --output-format vcf \
    --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $1 2> err.txt || { echo "FAILED $1"; tail -3 err.txt; return; }
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" $1 $2 $(cat t.txt) $(stat -c %s o.vcf) $(ls $1 | wc -l) | tee -a $OUT/times.tsv
}
printf "config\tthreads\twall_s\tuser_s\tmaxrss_kb\tout_bytes\tfiles\n" > $OUT/times.tsv
run base 16   # warm-up of the transcript cache
for rep in 1 2; do for cfg in base plus_gnomad_genomes plus_onekg plus_gerp plus_dbnsfp plus_encode all; do run $cfg 16; done; done
# the bcftools pipe of the benchmark, on the full stack
s=$(date +%s.%N); env RAYON_NUM_THREADS=16 $FV annotate --input in.vcf --output - --output-format vcf --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir all 2> /dev/null | bcftools view -Ob -o p.bcf -W --threads 8
echo "all_piped_to_bcftools $(echo "$(date +%s.%N) - $s" | bc)" | tee -a $OUT/times.tsv
# per-source fields added, on one record that has dbNSFP values
grep -v '^#' o.vcf | awk -F'\t' 'length($8) > 2000' | head -1 | cut -f8 | tr ';' '\n' | cut -c1-120 > $OUT/long_record_fields.txt
grep -v '^#' o.vcf | awk -F'\t' '{n+=length($8)} END {print "mean INFO length " n/NR}' | tee -a $OUT/info.txt
echo DONE
