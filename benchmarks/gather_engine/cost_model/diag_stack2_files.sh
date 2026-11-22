#!/bin/bash
# Does fastVEP's cost grow with the number of files a source has, even for files of other
# chromosomes? chr22 input; the mission stack plus gnomAD genomes as chr22's 3 pieces
# or as all 166 pieces of the genome; likewise onekg (1 vs 23 files) and gerp (1 vs 24).
# Usage: diag_stack2_files.sh [sample] [chromosome]
set -uo pipefail
SAMPLE=${1:-KPGP-00006}; C=${2:-22}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
R=/mnt/data/tmp/mission-wgs2; M16=$R/sa16
OUT=/results/fastvep05/cost_model/diag_stack2; ls /results/fastvep05/cost_model > /dev/null 2>&1; mkdir -p $OUT
W=$(mktemp -d /mnt/data/tmp/diag2f-XXXX); trap 'rm -rf $W' EXIT; cd $W
bcftools view -r chr$C $R/$SAMPLE/in.bcf -Ov -o in.vcf || exit 1
genome=$(ls $M16 | grep -v '_chr' | grep -v encode)
base="alphamissense dbsnp gnomad phylop revel spliceai tracka trackb"
mk() {  # dir pattern-of-extra-files: base sources of this chromosome plus the extra files
  local d=$1 x=$2; mkdir -p $d
  for g in $genome; do ln -s $M16/$g $d/; done
  for s in $base; do ln -s $M16/${s}_chr${C}.* $d/; done
  for f in $(ls $M16 | grep -E "$x"); do ln -s $M16/$f $d/; done
}
mk gg_here "^gnomad_genomes_chr${C}_"
mk gg_all "^gnomad_genomes_"
mk kg_here "^onekg_chr${C}\."
mk kg_all "^onekg_"
mk gerp_here "^gerp_chr${C}\."
mk gerp_all "^gerp_"
run() {
  /usr/bin/time -f "%e %U %M" -o t.txt env RAYON_NUM_THREADS=16 $FV annotate --input in.vcf --output o.vcf --output-format vcf \
    --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $1 2> err.txt || { echo "FAILED $1"; tail -3 err.txt; return; }
  printf "%s\t%s\t%s\t%s\n" $1 "$(cat t.txt | tail -1)" $(md5sum < o.vcf | cut -c1-12) $(ls $1 | wc -l) | tee -a $OUT/files_times.tsv
}
printf "config\twall user maxrss_kb\tmd5\tfiles\n" > $OUT/files_times.tsv
run gg_here
for rep in 1 2; do for c in gg_here gg_all kg_here kg_all gerp_here gerp_all; do run $c; done; done
echo DONE
