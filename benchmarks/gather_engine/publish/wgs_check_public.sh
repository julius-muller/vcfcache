#!/bin/bash
# Whole-genome check of a public cache built by build_public_cache.sh: one KPGP genome, cached
# and --uncached, 2 repeats each, times and a record-by-record comparison of the last pair.
# Usage: wgs_check_public.sh <deployment> vep|fastvep [sample]
set -uo pipefail
DEP=$1; TOOL=$2; SAMPLE=${3:-KPGP-00128}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp/$SAMPLE.GRCh38.small_variants.vcf.gz
OUT=/results/publish/$TOOL; ROOT=$OUT/gnomad_v4.1_GRCh38_joint_af001; RC=$OUT/recipe
NAME=$([ $TOOL = vep ] && echo vep116.3_e || echo fastvep0.5.0_core); C=$ROOT/cache/$NAME
ls /results /results/publish > /dev/null 2>&1
vc() { env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" "$@"; }
W=$(mktemp -d /mnt/data/tmp/pubwgs-$TOOL-XXXX); trap 'rm -rf $W' EXIT; cd $W
sed "s#^temp_dir: .*#temp_dir: \"$W\"#" $RC/params.yaml > params.yaml
bcftools view $KP -Ob -W -o in.bcf; N=$(bcftools index -n in.bcf)
echo "node $(hostname), $TOOL, $SAMPLE, $N variants"
for rep in $(seq ${REPS:-2}); do for arm in cached uncached; do
  rm -rf $arm.bcf* $arm.stats; s=$(date +%s.%N)
  vc annotate -a $C -i in.bcf -o $arm.bcf -y params.yaml --stats-dir $arm.stats $([ $arm = uncached ] && echo --uncached) > $arm.log 2>&1 || { echo "FAILED $arm"; tail -10 $arm.log; exit 1; }
  miss=$(grep -ho "missing_variants: [0-9]*" -r $arm.stats | head -1 | cut -d' ' -f2)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" $arm $N "${miss:-NA}" "$(echo "$(date +%s.%N) - $s" | bc)" $SAMPLE $(hostname) | tee -a $OUT/wgs_check.tsv
done; done
$PY -I /results/fastvep05/cost_model/vep116_compare.py uncached.bcf cached.bcf | sed "s/^/$SAMPLE\t/" | tee -a $OUT/wgs_compare.tsv
echo DONE
