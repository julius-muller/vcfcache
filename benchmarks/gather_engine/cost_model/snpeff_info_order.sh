#!/bin/bash
# Why cached SnpEff output differs from direct output in INFO tag order (snpeff_compare.py: "identical
# (INFO tag order only)", 6.5% of genome records, 92-99% of exome records): the ACMG panel of a sample
# direct and from its rerun cache (caches as built by wgs_mission_snpeff.sh on this node); prints the
# INFO keys of the first records whose order differs, and the INFO header order of input and cache.
# Result (job 3727, KPGP-00006 panel, 284 of 284 records): headers identical; SnpSift writes the tags
# of a source in the order of the source's record (ClinVar CLNDN,CLNSIG,REVSTAT; gnomAD AC,AF,AF_afr,
# ...) but declares them in another order (CLNSIG,REVSTAT,CLNDN; AF_afr,AC,AF,...). vcfcache transfers
# the cached tags with bcftools annotate, which writes them in header order. Records with a multi-tag
# source (ClinVar, gnomAD exomes: most exome and panel records, 6-7% of a genome) therefore list the
# same tags and values in another order; the values are identical.
# Usage: snpeff_info_order.sh <deployment> <sample>
set -uo pipefail
DEP=$1; S=$2
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
W=/mnt/data/tmp/snpeff-mission/$S; test -e $W/CACHES || { echo "no caches for $S on $(hostname)"; exit 1; }
T=$(mktemp -d /mnt/data/tmp/seord-XXXX); trap 'rm -rf $T' EXIT; cd $T
bcftools view -R /results/fastvep05/cost_model/regions/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed $W/in.bcf -Ob -W -o panel.bcf
C=$W/caches/rerun/cache/rerun
for arm in direct cached; do
  env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $C -i panel.bcf -o $arm.bcf \
    -y $C/params.snapshot.yaml --statistics light --stats-dir st.$arm $([ $arm = direct ] && echo --uncached) > log.$arm 2>&1 || { tail log.$arm; exit 1; }
done
keys() { bcftools view -H $1 | awk -F'\t' '{n=split($8,a,";"); s=""; for(i=1;i<=n;i++) {sub(/=.*/,"",a[i]); s=s (i>1?",":"") a[i]}; print $1":"$2":"$4">"$5"\t"s}'; }
keys direct.bcf > d.keys; keys cached.bcf > c.keys
echo "records whose INFO order differs: $(paste d.keys c.keys | awk -F'\t' '$2!=$4' | wc -l) of $(wc -l < d.keys)"
paste d.keys c.keys | awk -F'\t' '$2!=$4 {print "direct\t" $1 "\t" $2; print "cached\t" $3 "\t" $4; if (++n==3) exit}'
for f in panel.bcf direct.bcf cached.bcf $C/vcfcache_annotated.bcf; do echo "INFO header $f:"; bcftools view -h $f | grep -o '^##INFO=<ID=[^,]*' | cut -d= -f3 | tr '\n' ',' | cut -c1-900; echo; done
