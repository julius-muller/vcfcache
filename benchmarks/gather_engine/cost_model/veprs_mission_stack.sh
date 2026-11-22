#!/bin/bash
# The fastVEP mission stack's per-chromosome content as vep-rs 0.5.0 plugin
# tables (vep-rs has no --custom; its built-in dbNSFP plugin reads any tabix
# table with #chr, pos(1-based), ref, alt columns). Same inputs as
# wgs_mission_stack.sh: gnomAD v4.1 exomes (with SpliceAI and PhyloP taken from
# it), dbSNP and the two mock tracks. Writes headerless parts; veprs_mission_join.sh
# joins them. Usage: veprs_mission_stack.sh <chromosome without chr>
set -uo pipefail
C=$1
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
FA=/mnt/data/apps/fastvep/0.5.0/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
A=/results/fastvep05/acmg; S=$A/sources
ls /results /results/veprs05 > /dev/null 2>&1   # clears a stale "no such directory" on some nodes
OUT=/results/veprs05/mission; mkdir -p $OUT/parts $OUT/done
[ -e $OUT/done/chr$C ] && exit 0
W=$(mktemp -d /mnt/data/tmp/veprs-stack-$C-XXXX); trap 'rm -rf $W' EXIT; cd $W
log() { echo "[$(date -Is)] chr$C $*"; }
if [ $C = 1 ]; then V=$S/gnomad.exomes.v4.1.sites.chr1.vcf.bgz; else
  V=$W/gnomad.vcf.bgz; G=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr$C.vcf.bgz
  curl -fsSL --retry 10 --retry-delay 30 -o $V $G && curl -fsSL --retry 10 -o $V.tbi $G.tbi || { echo "download failed"; exit 1; }
  log "downloaded $(du -h $V | cut -f1)"
fi
fail=0
# gnomAD: the fields fastVEP's gnomad source keeps (frequency, counts, ten groups); SpliceAI max and PhyloP
POPS="afr amr asj eas fin mid nfe remaining sas"
F='%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%INFO/nhomalt'; for p in $POPS; do F="$F\t%INFO/AF_$p"; done
( bcftools norm -m- $V -Ou --threads 2 | bcftools query -f "$F\n" | bgzip -@2 > $W/gnomad.part.gz ) & p1=$!
( bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/spliceai_ds_max\t%INFO/phylop\n' -i 'INFO/spliceai_ds_max!="." || INFO/phylop!="."' $V > $W/sp.tsv
  # SpliceAI as the plugin's VCF (max score in all four delta fields, as wgs_mission_stack.sh's extract does)
  awk -F'\t' 'BEGIN{OFS="\t"} $5!="." {print $1,$2,".",$3,$4,".",".","SpliceAI="$4"|.|"$5"|"$5"|"$5"|"$5"|0|0|0|0"}' $W/sp.tsv | bgzip > $W/spliceai.part.gz
  # PhyloP is positional in fastVEP: every SNV allele at a scored position
  awk -F'\t' 'BEGIN{OFS="\t"} $6!="." && !seen[$2]++ && length($3)==1 {for(i=1;i<=4;i++){a=substr("ACGT",i,1); if(a!=$3) print $1,$2,$3,a,$6}}' $W/sp.tsv | bgzip > $W/phylop.part.gz ) & p2=$!
( n=$(case $C in X) echo 23;; Y) echo 24;; *) echo $C;; esac)
  acc=$(bcftools view -h $S/GCF_000001405.40.gz | grep -o "NC_0000$(printf %02d $n)\.[0-9]*" | head -1)
  [ -n "$acc" ] || { echo "FAILED no dbSNP accession"; exit 1; }
  bcftools view -r $acc $S/GCF_000001405.40.gz -Ou | bcftools norm -m- -Ou | bcftools query -f "chr$C\t%POS\t%REF\t%ALT\t%ID\n" | bgzip -@2 > $W/dbsnp.part.gz ) & p3=$!
( $PY - $FA $C <<'PYE' | bgzip -@2 > $W/track.part.gz
import sys
fa, c = sys.argv[1:]
for line in open(fa + ".fai"):
    f = line.split("\t")
    if f[0] == c:
        n, off, bases, width = map(int, f[1:5]); break
out = sys.stdout
with open(fa, "rb") as s:
    for pos in range(1, n + 1, 93):
        s.seek(off + (pos - 1) // bases * width + (pos - 1) % bases)
        r = s.read(1).decode().upper()
        if r not in "ACGT":
            continue
        for k, a in enumerate("ACGT"):
            if a != r:
                h = (pos * 48271 + k * 40503) % 1000003
                out.write("chr%s\t%d\t%s\t%s\t%.3f\t%.2f\n" % (c, pos, r, a, (h % 991) / 991, (h % 97) / 97))
PYE
) & p4=$!
for p in $p1 $p2 $p3 $p4; do wait $p || fail=1; done
for f in gnomad spliceai phylop dbsnp track; do
  [ $(zcat $W/$f.part.gz | head -1 | wc -l) = 1 ] || { [ $f = spliceai ] && [ $C = Y ] && continue; echo "missing $f"; fail=1; }
done
[ $fail = 0 ] || { echo "chr$C incomplete"; exit 1; }
for f in gnomad spliceai phylop dbsnp track; do
  for try in 1 2 3 4 5; do ls $OUT/parts > /dev/null 2>&1; cp $W/$f.part.gz $OUT/parts/${f}_chr$C.gz && break; sleep 20; done
  [ -s $OUT/parts/${f}_chr$C.gz ] || { echo "chr$C not saved"; exit 1; }
done
touch $OUT/done/chr$C
log "done: $(for f in gnomad spliceai phylop dbsnp track; do echo -n "$f $(zcat $W/$f.part.gz | wc -l) "; done)"
