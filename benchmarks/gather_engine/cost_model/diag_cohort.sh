#!/bin/bash
# Why is the cohort cache (fewer misses) slower than the 1% blueprint cache
# when the misses are annotated from the 65 kb sources? Reuses the KPGP-00128
# caches left on sl-w1 by wgs_mission.sh; keeps vcfcache's work directory and
# times each step, then times fastVEP alone on each arm's misses.
# Usage: diag_cohort.sh <deployment>
set -uo pipefail
DEP=$1
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
R=/mnt/data/tmp/mission-wgs; cd $R || exit 1
ls /results/fastvep05 > /dev/null 2>&1
OUT=/results/fastvep05/cost_model/diag_cohort; mkdir -p $OUT
D=$R/diag; rm -rf $D; mkdir -p $D
T="/usr/bin/time -f %e\t%U\t%S\t%P\t%M"
echo -e "what\tarm\trep\twall\tuser\tsys\tcpu\tmaxrss_kb" > $OUT/times.tsv
for rep in 1 2; do
  for k in af001 cohort af001_narrow cohort_narrow; do
    C=$R/caches/$k/cache/$k; S=$D/st.$k.$rep; rm -rf $S out.bcf*
    $T -o $D/time.$k.$rep env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $C -i in.bcf -o $D/out.$k.bcf \
      -y $C/params.snapshot.yaml --statistics light --stats-dir $S -vv --debug > $D/log.$k.$rep 2>&1 || { echo "FAILED $k"; tail $D/log.$k.$rep; }
    echo -e "vcfcache\t$k\t$rep\t$(tail -1 $D/time.$k.$rep)" | tee -a $OUT/times.tsv
    cp $D/log.$k.$rep $OUT/
    if [ $rep = 1 ]; then
      bcftools concat --naive $(ls $S/*/work/annotate/shard-*/misses.bcf | sort -V) -Ob -o $D/miss.$k.bcf
      echo "$k misses $(bcftools view -H $D/miss.$k.bcf | wc -l)" | tee -a $OUT/notes.txt
    fi
    rm -rf $S $D/out.$k.bcf*
  done
done
# fastVEP alone on each arm's misses, from the advised 65 kb sources and from the narrow copies
for rep in 1 2; do
  for k in af001 cohort; do
    for sa in sa16 sa12; do
      bcftools view $D/miss.$k.bcf -Ov > $D/in.vcf
      $T -o $D/t env RAYON_NUM_THREADS=8 $FV annotate --input $D/in.vcf --output $D/o.vcf --output-format vcf \
        --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $R/$sa 2> $D/fv.err
      echo -e "fastvep_misses\t$k.$sa\t$rep\t$(tail -1 $D/t)" | tee -a $OUT/times.tsv
    done
  done
done
# the hit lookup alone: the overlay of each cache over the whole input, no annotator
for k in af001 cohort; do
  C=$R/caches/$k/cache/$k/vcfcache_annotated.bcf
  $T -o $D/t bcftools annotate --no-version --pair-logic exact -i 'N_ALT=1' -k -a $C -c INFO in.bcf -Ob -o $D/ov.bcf --threads 8
  echo -e "overlay_1thread_job\t$k\t1\t$(tail -1 $D/t)" | tee -a $OUT/times.tsv
  echo "$k cache records $(bcftools index -n $C)" | tee -a $OUT/notes.txt
done
echo DONE
