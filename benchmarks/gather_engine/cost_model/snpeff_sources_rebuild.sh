#!/bin/bash
# Rebuilds /results/snpeff/mission/src (the plain SnpSift sources with their .sidx) after it was
# deleted to free disk (2026-11-27), from the files it was made from: the --custom VCFs of the
# Ensembl VEP 116.3 mission and the SpliceAI-style VCFs of the vep-rs mission, decompressed with
# bcftools view as before; dbSNP in the five chromosome groups of snpeff_dbsnp_split.sh. Each
# .sidx is built by SnpSift on a one-record input, as snpsift_memory_probe.sh and
# snpeff_dbsnp_split.sh did (16 GB heap, 26 GB for the dbSNP groups). The sizes of the plain VCFs
# and indexes are compared with the listing of the first build (memory_probe/src.txt).
# Usage: snpeff_sources_rebuild.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
M=/results/snpeff/mission; S=$M/src; VM=/results/vep116/mission
ls /results /results/snpeff $M $VM > /dev/null 2>&1   # clears a stale "no such directory" on some nodes
R=/mnt/data/tmp/snpeff-mission; SRC=$R/src; mkdir -p $SRC
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname)"
test -e $S/READY && { echo "shared sources exist"; exit 0; }
ERR=$R/rebuild.err; CU=$SRC; T=$SRC; . /results/fastvep05/cost_model/snpeff_stack.sh
s=$(date +%s); pids=()
for f in $VM/custom/{clinvar,gnomad,phylop,trackA,trackB}.vcf.gz /results/veprs05/mission/tables/spliceai_{snv,indel}.vcf.gz; do
  [ -s $SRC/$(basename $f .gz) ] || { bcftools view $f -Ov -o $SRC/$(basename $f .gz) & pids+=($!); }
done
DB=$VM/custom/dbsnp.vcf.gz
G1="chr1 chr2 chr3"; G2="chr4 chr5 chr6 chr7"; G3="chr8 chr9 chr10 chr11 chr12"; G4="chr13 chr14 chr15 chr16 chr17 chr18"
G5=$(tabix -l $DB | grep -v -x -F -f <(echo $G1 $G2 $G3 $G4 | tr ' ' '\n') | tr '\n' ' ')
for g in 1 2 3 4 5; do eval "regs=\$G$g"
  [ -s $SRC/dbsnp_g$g.vcf ] || { bcftools view $DB -r $(echo $regs | tr ' ' ',') -Ov -o $SRC/dbsnp_g$g.vcf & pids+=($!); }
done
for p in ${pids[@]}; do wait $p || { echo "FAILED decompress"; exit 1; }; done
log "decompressed: $(du -sh $SRC | cut -f1), $(( $(date +%s) - s )) s"
bcftools view -h $SRC/clinvar.vcf > $R/one.vcf; bcftools view -H $SRC/clinvar.vcf | head -1 >> $R/one.vcf
for f in clinvar spliceai_indel spliceai_snv trackA trackB phylop gnomad dbsnp_g1 dbsnp_g2 dbsnp_g3 dbsnp_g4 dbsnp_g5; do
  [ $SRC/$f.vcf.sidx -nt $SRC/$f.vcf ] && continue
  h=16g; [[ $f == dbsnp_* ]] && h=26g
  s=$(date +%s)
  $JAVA -Xmx$h -jar $SE/SnpSift.jar annotate -noLog -c $SE/snpEff/snpEff.config -noId $SRC/$f.vcf $R/one.vcf > /dev/null 2>> $ERR \
    || { echo "FAILED sidx $f"; tail -n 3 $ERR; exit 1; }
  log "$f.vcf.sidx: $(stat -c %s $SRC/$f.vcf.sidx) bytes, $(( $(date +%s) - s )) s"
done
# sizes against the first build (dbSNP groups were not in that listing)
awk 'NF == 9 && $9 ~ /\.vcf/ {print $9, $5}' $M/memory_probe/src.txt | sort -u | while read f n; do
  [ -e $SRC/$f ] || continue; m=$(stat -c %s $SRC/$f); echo -e "$f\t$n\t$m\t$([ $n = $m ] && echo same || echo DIFFERS)"
done | tee $R/sizes.tsv
grep -q DIFFERS $R/sizes.tsv && echo "WARNING: sizes differ from the first build"
mkdir -p $S; s=$(date +%s)
cp --preserve=timestamps $SRC/*.vcf $SRC/*.vcf.sidx $S/ || { echo "FAILED copy"; exit 1; }
for x in $S/*.sidx; do [ $x -nt ${x%.sidx} ] || touch $x; done
echo split > $S/DBSNP_MODE; cp $R/sizes.tsv $M/rebuild_sizes.tsv; touch $S/READY
log "shared sources: $(du -sh $S | cut -f1), copied in $(( $(date +%s) - s )) s"; echo DONE
