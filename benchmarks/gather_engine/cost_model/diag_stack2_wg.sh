#!/bin/bash
# Is the larger stack slow because its sources (35 GB at 65 kb blocks) no longer fit
# into a 30 GB node's page cache? Whole genome, fastVEP alone, 16 threads:
# the full stack (stopped after 40 min) and the stack without gnomAD genomes (19 GB).
# Records wall, user, major page faults, and disk reads per minute (vmstat).
# Usage: diag_stack2_wg.sh [sample]
set -uo pipefail
SAMPLE=${1:-KPGP-00128}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
R=/mnt/data/tmp/mission-wgs2; M16=$R/sa16
OUT=/results/fastvep05/cost_model/diag_stack2; ls /results/fastvep05/cost_model > /dev/null 2>&1; mkdir -p $OUT
W=$(mktemp -d /mnt/data/tmp/diag2wg-XXXX); trap 'rm -rf $W' EXIT; cd $W
free -g | tee $OUT/wg_mem.txt
bcftools view $R/$SAMPLE/in.bcf -Ov -o in.vcf || exit 1
mkdir nogg; for f in $M16/*; do case $(basename $f) in gnomad_genomes_*) ;; *) ln -s $f nogg/;; esac; done
run() {  # name sa-dir timeout
  vmstat -S M 60 > $OUT/wg_vmstat_$1.txt & local vm=$!
  /usr/bin/time -f "%e %U %S %F %M" -o t.txt timeout $3 env RAYON_NUM_THREADS=16 $FV annotate --input in.vcf --output /dev/null \
    --output-format vcf --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $2 2> err.txt; local rc=$?
  kill $vm
  printf "%s\trc=%s\t%s\n" $1 $rc "$(cat t.txt | tail -1)" | tee -a $OUT/wg_times.tsv
}
printf "config\trc\twall user sys majflt maxrss_kb\n" > $OUT/wg_times.tsv
run nogg nogg 3000
run all $M16 2400
run nogg_again nogg 3000
echo DONE
