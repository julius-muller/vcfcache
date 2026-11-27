#!/bin/bash
# Submits the 40-genome campaign (cohort40.sh) after cohort40_prep.sh: the cache parts of each
# annotator (VEP 116.3 and SnpEff reuse their mission's bp and new parts), then one job per
# annotator and evaluation genome, each waiting for its annotator's parts. The genomes go in order of
# record count, largest first, so that a genome that does not fit (SnpEff's memory grows with the
# input) shows up early. SnpEff waits for its sources too (job SRC, snpeff_sources_rebuild.sh).
# Usage (on the head node, from the cost_model folder): cohort40_submit.sh <deployment> [SRC job id] [annotators]
set -uo pipefail
DEP=$1; SRC=${2:-}; ANS=${3:-"fastvep veprs vep116 snpeff"}
PR=/mnt/data/slurm-results/cohort40/prep; L=/results/cohort40/logs
test -e $PR/READY || { echo "prep missing"; exit 1; }
NX=$(ls $PR/extra-*.regions | wc -l)
# a finished job can no longer be named as a dependency: then the sources must be there
if [ -n "$SRC" ] && [ -z "$(squeue -h -j $SRC 2> /dev/null)" ]; then
  test -e /mnt/data/slurm-results/snpeff/mission/src/READY || { echo "SnpEff sources missing"; exit 1; }; SRC=""; fi
sb() { sbatch --parsable -p benchmark --exclusive "$@"; }
declare -A PARTS LIM PLIM
LIM=([fastvep]=4:00:00 [veprs]=5:00:00 [vep116]=12:00:00 [snpeff]=8:00:00)
PLIM=([fastvep]=1:00:00 [veprs]=2:00:00 [vep116]=3:00:00 [snpeff]=3:00:00)
for an in $ANS; do
  sets="extra"; case $an in fastvep|veprs) sets="bp new extra";; esac
  dep=""; [ $an = snpeff ] && [ -n "$SRC" ] && dep="--dependency=afterok:$SRC"
  ids=""
  for s in $sets; do for f in $PR/$s-*.regions; do k=$(basename $f .regions); k=${k#$s-}
    id=$(sb -t ${PLIM[$an]} $dep -J c40p-$an-$s$k -o $L/part-$an-$s-$k-%j.out cohort40.sh $DEP $an part $s $k) || exit 1
    ids="$ids:$id"; done; done
  PARTS[$an]=$ids; echo "$an parts: $(echo $ids | tr ':' ' ' | wc -w) jobs ($ids)"
done
for s in $(awk -F'\t' '$5 == "GRCh38" && $3 == "evaluation" {print $7, $2}' $PR/samples.tsv | sort -k1,1nr | cut -d' ' -f2); do
  for an in $ANS; do
    d="--dependency=afterok${PARTS[$an]}"; [ $an = snpeff ] && [ -n "$SRC" ] && d="$d:$SRC"
    id=$(sb -t ${LIM[$an]} $d -J c40-$an-$s -o $L/measure-$an-$s-%j.out cohort40.sh $DEP $an measure $s) || exit 1
    echo -e "$an\t$s\t$id"
  done
done | tee -a /mnt/data/slurm-results/cohort40/jobs.tsv
echo "extra parts per annotator: $NX"
