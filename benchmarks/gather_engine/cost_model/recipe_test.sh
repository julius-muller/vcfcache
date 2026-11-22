#!/bin/bash
# Functional test of vcfcache/recipes/grch38_fastvep0.5 on chr22: build a cache from chr22 of the
# gnomAD AF >= 0.1% blueprint with params.yaml (65 kb sources), annotate chr22 of one sample with
# params.annotate.yaml (1 kb sources), and compare with direct annotation (--uncached, 65 kb).
# Uses the node-local sources of wgs_mission.sh N10=1 (sa16, sa10).
# Usage: [STACK2=1] recipe_test.sh <deployment> <recipe dir> [sample]
set -uo pipefail
DEP=$1; RC=$2; SAMPLE=${3:-KPGP-00128}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
FVD=/mnt/data/apps/fastvep/0.5.0
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
ls /results/fastvep05 /results/fastvep05/cost_model > /dev/null 2>&1
R=/mnt/data/tmp/mission-wgs${STACK2:+2}f
for f in $R/SOURCES10 $R/sa16 $R/sa10 $RC/annotation.yaml $RC/params.yaml $RC/params.annotate.yaml; do test -e $f || { echo "missing $f"; exit 1; }; done
W=$(mktemp -d /mnt/data/tmp/recipe-XXXX); trap 'rm -rf $W' EXIT; cd $W; mkdir tmp
vc() { env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" "$@"; }
digest() { bcftools view -H $1 | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12; }
t() { local s=$(date +%s.%N); "$@" > last.log 2>&1 || { echo "FAILED: $*"; tail -20 last.log; exit 1; }; echo "$(date +%s.%N) - $s" | bc; }
for p in params params.annotate; do
  sed -e "s#^annotation_tool_cmd: .*#annotation_tool_cmd: \"$FVD/bin/fastvep\"#" -e "s#^tool_version_command: .*#tool_version_command: \"$FVD/bin/fastvep --version\"#" \
      -e "s#^transcript_cache: .*#transcript_cache: \"$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache\"#" \
      -e "s#^fasta: .*#fasta: \"$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa\"#" -e "s#^temp_dir: .*#temp_dir: \"$W/tmp\"#" \
      -e "s#/data/fastvep/sa/GRCh38_65kb#$R/sa16#" -e "s#/data/fastvep/sa/GRCh38_1kb#$R/sa10#" $RC/$p.yaml > $p.yaml
done
grep -H "sa_dir" params.yaml params.annotate.yaml
bcftools view -r chr22 $BP -Ob -W -o bp22.bcf
bcftools view -r chr22 $KP/$SAMPLE.GRCh38.small_variants.vcf.gz -Ob -W -o in.bcf
echo "blueprint chr22 $(bcftools index -n bp22.bcf) records, sample chr22 $(bcftools index -n in.bcf)"
echo "blueprint-init $(t vc blueprint-init -i bp22.bcf -o $W/bp -y params.yaml) s"
echo "cache-build $(t vc cache-build --db $W/bp -n fv -a $RC/annotation.yaml -y params.yaml) s"
C=$(dirname $(find $W/bp -name vcfcache_annotated.bcf | head -1)); echo "cache $C"
grep -c "sa_dir" $C/annotation.yaml > /dev/null && grep "sa-dir" $C/annotation.yaml
# twice, so that the second round reads every source folder warm
for rep in 1 2; do rm -rf st1 st2 st3 direct.bcf* c65.bcf* c1k.bcf*
  echo "round $rep: direct (uncached, 65 kb) $(t vc annotate -a $C -i in.bcf -o direct.bcf -y params.yaml --uncached --stats-dir st1) s"
  echo "round $rep: cached, misses at 65 kb  $(t vc annotate -a $C -i in.bcf -o c65.bcf -y params.yaml --stats-dir st2) s"
  echo "round $rep: cached, misses at 1 kb   $(t vc annotate -a $C -i in.bcf -o c1k.bcf -y params.annotate.yaml --stats-dir st3) s"
done
grep -rho "missing_variants: [0-9]*" st3 | head -1
grep -rh "^sa_dir" st2 st3
for f in direct c65 c1k; do echo "$f $(bcftools index -n $f.bcf) $(digest $f.bcf)"; done
[ $(digest direct.bcf) = $(digest c1k.bcf) ] && [ $(digest direct.bcf) = $(digest c65.bcf) ] && echo "IDENTICAL" || echo "DIFFERENT"
echo DONE
