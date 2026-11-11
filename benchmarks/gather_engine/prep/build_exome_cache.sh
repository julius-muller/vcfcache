#!/bin/bash
# Restrict one vcfcache cache to an exome capture BED (same records, same recipe).
# Usage: build_exome_cache.sh <cache dir> <bed>
set -euo pipefail
SRC=$1; BED=$2
name=$(basename $SRC)_twist
ROOT=/results/exome-caches/$name; C=$ROOT/cache/$name
test -s $SRC/vcfcache_annotated.bcf && test -s $BED
rm -rf $ROOT.tmp; mkdir -p $ROOT.tmp/cache/$name $ROOT.tmp/blueprint
cp -r $(dirname $(dirname $SRC))/workflow $ROOT.tmp/workflow
for f in annotation.yaml params.snapshot.yaml build_provenance.json; do [ -e $SRC/$f ] && cp $SRC/$f $ROOT.tmp/cache/$name/; done
start=$(date +%s)
# Default -R overlap, as used to extract the WES inputs from the genomes.
bcftools view -R $BED $SRC/vcfcache_annotated.bcf -Ob -W --threads 8 -o $ROOT.tmp/cache/$name/vcfcache_annotated.bcf
n=$(bcftools index -n $ROOT.tmp/cache/$name/vcfcache_annotated.bcf); test "$n" -gt 0
cat > $ROOT.tmp/cache/$name/exome_provenance.json <<J
{"source_cache": "$SRC", "source_records": $(bcftools index -n $SRC/vcfcache_annotated.bcf), "bed": "$BED",
 "bed_sha256": "$(sha256sum $BED | cut -d' ' -f1)", "records": $n, "seconds": $(( $(date +%s) - start )),
 "cache_sha256": "$(sha256sum $ROOT.tmp/cache/$name/vcfcache_annotated.bcf | cut -d' ' -f1)"}
J
rm -rf $ROOT; mv $ROOT.tmp $ROOT
cat $C/exome_provenance.json
