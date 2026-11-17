#!/bin/bash
# Downloads for the larger mission stack: dbNSFP 5.4a (academic; its links are
# in private files next to the data, never in this repository), GERP (Ensembl
# 91 mammals), ENCODE cCREs, and UCSC's bigWigToBedGraph. Usage: dl_sources.sh
set -uo pipefail
ls /results > /dev/null 2>&1
D=/results/sources; mkdir -p $D/dbnsfp $D/gerp $D/encode /results/fastvep05/tools
log() { echo "[$(date -Is)] $*"; }
get() { curl -fsSL --retry 20 --retry-delay 30 -C - -o "$2" "$1"; }
U=$(cat $D/dbnsfp/URL_GRCH38)
( get $U.tbi $D/dbnsfp/dbNSFP5.4a_grch38.gz.tbi && get $U.md5 $D/dbnsfp/dbNSFP5.4a_grch38.gz.md5 && get $U $D/dbnsfp/dbNSFP5.4a_grch38.gz &&
  log "dbNSFP $(du -h $D/dbnsfp/dbNSFP5.4a_grch38.gz | cut -f1)" &&
  (cd $D/dbnsfp && a=$(md5sum dbNSFP5.4a_grch38.gz | cut -d' ' -f1) && grep -q $a dbNSFP5.4a_grch38.gz.md5 && log "dbNSFP md5 ok") ) & a=$!
( get $(cat $D/dbnsfp/URL_GENE) $D/dbnsfp/dbNSFP5.4_gene.gz && log "dbNSFP gene table" ) & b=$!
( get https://ftp.ensembl.org/pub/release-113/compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw $D/gerp/gerp_91_mammals.GRCh38.bw &&
  log "GERP $(du -h $D/gerp/gerp_91_mammals.GRCh38.bw | cut -f1)" ) & c=$!
( get https://downloads.wenglab.org/V3/GRCh38-cCREs.bed $D/encode/GRCh38-cCREs.bed && log "cCREs $(wc -l < $D/encode/GRCh38-cCREs.bed)" ) & d=$!
( get https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph /results/fastvep05/tools/bigWigToBedGraph &&
  chmod +x /results/fastvep05/tools/bigWigToBedGraph && log "bigWigToBedGraph" ) & e=$!
f=0; for p in $a $b $c $d $e; do wait $p || f=1; done
[ $f = 0 ] && touch $D/DOWNLOADED && echo DONE || echo "a download failed"
