#!/bin/bash
# The public Ensembl VEP 116.3 cache (see build_public_cache.sh) in parts: one VEP process writes
# about 400 records/s whatever --fork is, so the 18.9 M-record blueprint takes about 13 h, longer
# than the 12 h job limit. Each part runs vcfcache cache-build with the published recipe on a slice
# of the blueprint (the 12 slices of wgs_mission_vep116.sh, cut between positions); join
# concatenates the slices in order into one cache and records the parts in build_provenance.json.
# Usage: build_public_vep_parts.sh <deployment> part <k> | join
set -uo pipefail
DEP=$1; MODE=$2; K=${3:-}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
BPROOT=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001
SL=/results/vep116/mission/parts   # bp-<k>.regions
OUT=/results/publish/vep; RC=$OUT/recipe; NAME=vep116.3_e; PP=$OUT/parts; ROOT=$OUT/gnomad_v4.1_GRCh38_joint_af001
NP=$(ls $SL/bp-*.regions | wc -l)
ls /results /results/publish /results/vep116/mission > /dev/null 2>&1
vc() { env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" "$@"; }
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname), mode $MODE $K"
test -e $RC/annotation.yaml || { echo "recipe missing"; exit 1; }
mkdir -p $PP

if [ $MODE = part ]; then
  PR=$PP/root-$K
  [ -e $PR/cache/$NAME/.vcfcache_complete ] && { echo "part $K exists"; exit 0; }
  W=$(mktemp -d /mnt/data/tmp/pubvep-$K-XXXX); trap 'rm -rf $W' EXIT
  rm -rf $PR; mkdir -p $PR && cp -a $BPROOT/. $PR/ && rm -f $PR/zenodo_blueprint_provenance.json $PR/blueprint/vcfcache.bcf* || { echo FAILED copy; exit 1; }
  bcftools view -R $SL/bp-$K.regions --regions-overlap pos $BPROOT/blueprint/vcfcache.bcf -Ob -W -o $PR/blueprint/vcfcache.bcf || { echo FAILED slice; exit 1; }
  sed "s#^temp_dir: .*#temp_dir: \"$W\"#" $RC/params.yaml > $W/params.yaml
  s=$(date +%s)
  vc cache-build -d $PR -n $NAME -a $RC/annotation.yaml -y $W/params.yaml > $PP/build-$K.log 2>&1 || { echo "FAILED part $K"; tail -20 $PP/build-$K.log; exit 1; }
  n=$(bcftools index -n $PR/cache/$NAME/vcfcache_annotated.bcf); want=$(bcftools index -n $PR/blueprint/vcfcache.bcf); e=$(( $(date +%s) - s ))
  [ "$n" = "$want" ] || { echo "FAILED part $K: $n of $want records"; exit 1; }
  echo -e "$K\t$n\t$e\t$(hostname)" >> $PP/parts.tsv; log "part $K: $n records, ${e}s"; echo DONE; exit 0
fi

[ $MODE = join ] || { echo "unknown mode $MODE"; exit 1; }
for k in $(seq 0 $((NP-1))); do test -e $PP/root-$k/cache/$NAME/.vcfcache_complete || { echo "part $k missing"; exit 1; }; done
rm -rf $ROOT; mkdir -p $ROOT && cp -a $BPROOT/. $ROOT/ && rm -f $ROOT/zenodo_blueprint_provenance.json || { echo FAILED copy; exit 1; }
C=$ROOT/cache/$NAME; P0=$PP/root-0/cache/$NAME; mkdir -p $C/auxiliary
cp $P0/annotation.yaml $P0/params.snapshot.yaml $P0/blueprint_snapshot.info $C/
for k in $(seq 0 $((NP-1))); do
  mkdir -p $C/auxiliary/part-$k && cp -r $PP/root-$k/cache/$NAME/auxiliary/. $C/auxiliary/part-$k/ && cp $PP/root-$k/cache/$NAME/workflow.log $C/auxiliary/part-$k/
done
bcftools concat $(for k in $(seq 0 $((NP-1))); do echo $PP/root-$k/cache/$NAME/vcfcache_annotated.bcf; done) -Ob -W -o $C/vcfcache_annotated.bcf --threads 8 || { echo FAILED concat; exit 1; }
n=$(bcftools index -n $C/vcfcache_annotated.bcf); want=$(bcftools index -n $ROOT/blueprint/vcfcache.bcf)
[ "$n" = "$want" ] || { echo "FAILED join: $n of $want records"; exit 1; }
# every blueprint record, in blueprint order
cmp <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $C/vcfcache_annotated.bcf) <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $ROOT/blueprint/vcfcache.bcf) || { echo "FAILED join: records differ from the blueprint"; exit 1; }
tot=$(awk '{s+=$3} END {print s}' $PP/parts.tsv)
printf "task_id\tname\tstatus\texit\tduration\n1\tcache-build\tCOMPLETED\t0\t%ss (sum of %s parts)\n" $tot $NP > $C/${NAME}_trace.txt
$PY - $C $PP/parts.tsv $NP <<'P'
import json, sys
c, parts, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = [l.split("\t") for l in open(parts).read().split("\n") if l]
json.dump({"kind": "vcfcache cache-build of the published recipe, run on blueprint slices in parallel",
           "reason": "one Ensembl VEP process writes about 400 records/s; the whole blueprint exceeds the 12 h job limit",
           "slices": n, "slicing": "equal record counts, cut only between positions; concatenated in order",
           "parts": [dict(part=int(r[0]), records=int(r[1]), seconds=int(r[2]), node=r[3].strip()) for r in sorted(rows, key=lambda r: int(r[0]))]},
          open(f"{c}/build_provenance.json", "w"), indent=1)
P
env PYTHONPATH=$DEP VCFCACHE_GIT_COMMIT=${GIT_COMMIT:-unknown} $PY -c "
from pathlib import Path; from vcfcache.utils.completion import write_completion_flag
write_completion_flag(Path('$C'), 'cache-build', mode='cache-build')"
cat $C/.vcfcache_complete; du -sh $ROOT; log "joined $n records"; echo DONE
