#!/bin/bash
# fastVEP 0.5.0 everywhere, its two caches, all exome caches, then the cohort launch.
# Usage: submit_fastvep050.sh <deployment> <commit>
set -euo pipefail
DEP=$1; COMMIT=$2; V=0.5.0
P=$DEP/benchmarks/gather_engine/prep; L=/results/apps-staging/logs
sb() { sbatch --parsable --partition=benchmark --chdir=/tmp --output=$L/%x-%j.out "$@" | cut -d';' -f1; }
dl=$(sb --job-name=dl-fastvep-$V -c1 --mem=2000M --time=01:00:00 --wrap "bash $P/dl_fastvep.sh $V")
inst=()
for n in sl-w1 sl-w2 sl-w3 sl-w4 sl-w5 sl-w6; do
  inst+=($(sb --job-name=inst-fastvep-$V-$n --nodelist=$n -c8 --mem=16000M --time=02:00:00 --dependency=afterok:$dl --wrap "bash $P/install_fastvep.sh $V"))
done
all=$(IFS=:; echo "${inst[*]}")
c38=$(sb --job-name=fastvep-$V-cache-GRCh38 -c16 --mem=28000M --time=06:00:00 --dependency=afterok:$all --wrap "bash $P/build_fastvep_cache.sh GRCh38 $V")
c37=$(sb --job-name=fastvep-$V-cache-GRCh37 -c16 --mem=28000M --time=06:00:00 --dependency=afterok:$all --wrap "bash $P/build_fastvep_cache.sh GRCh37 $V")
B38=/results/richness/regions/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed
B37=/results/richness/regions/twist_human_core_exome_hg19.chr1-22-X.pad125.merged_GRCh37.bed
ex=()
ex+=($(sb --job-name=exome-fastvep38 -c8 --mem=8000M --time=01:00:00 --dependency=afterok:$c38 --wrap "bash $P/build_exome_cache.sh /results/fastvep05/caches/fastvep050_clinical_grch38_af001/cache/fastvep050_clinical_grch38_af001 $B38"))
ex+=($(sb --job-name=exome-fastvep37 -c8 --mem=8000M --time=01:00:00 --dependency=afterok:$c37 --wrap "bash $P/build_exome_cache.sh /results/fastvep05/caches/fastvep050_clinical_grch37_af001/cache/fastvep050_clinical_grch37_af001 $B37"))
for spec in \
  /results/veprs05/caches/veprs05_core/cache/veprs05_core:$B38 \
  /results/veprs05/caches/veprs05_clinical/cache/veprs05_clinical:$B38 \
  /results/vep116/caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep116.3_everything:$B38 \
  /results/snpeff/caches/snpeff_clinical_grch38_af001/cache/snpeff_clinical_grch38_af001:$B38 \
  /results/snpeff/caches/snpeff_clinical_grch37_af001/cache/snpeff_clinical_grch37_af001:$B37 \
  /results/snpeff/caches/snpeff_core_grch37_af001/cache/snpeff_core_grch37_af001:$B37; do
  ex+=($(sb --job-name=exome-$(basename ${spec%%:*}) -c8 --mem=8000M --time=01:00:00 --wrap "bash $P/build_exome_cache.sh ${spec%%:*} ${spec##*:}"))
done
deps=$(IFS=:; echo "${ex[*]}")
launch=$(sb --job-name=cohort-launch -c1 --mem=2000M --time=00:30:00 --dependency=afterok:$deps --wrap "bash $P/launch_cohort.sh $DEP $COMMIT /results/gather-engine/cohort-v116-fv050-20261111")
echo "download=$dl installs=$all caches=$c38,$c37 exome=$deps launch=$launch"
