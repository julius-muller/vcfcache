#!/bin/bash
# Mission-stack sources for one chromosome, twice: as fastVEP advises (dense
# sources in 65 kb blocks) into sa16, and in narrow blocks into sa12 for the
# misses of a cached run. Usage: wgs_mission_stack.sh <chromosome without chr>
set -uo pipefail
C=$1
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
A=/results/fastvep05/acmg; S=$A/sources
OUT=/results/fastvep05/mission; mkdir -p $OUT/sa16 $OUT/sa12 $OUT/done
# A chromosome marked done by the first version of this script lacks dbSNP.
ONLY=""
if [ -e $OUT/done/chr$C ]; then
  [ -e $OUT/sa16/dbsnp_chr$C.osa2 ] && [ -e $OUT/sa12/dbsnp_chr$C.osa2 ] && exit 0
  ONLY=dbsnp
fi
W=$(mktemp -d /mnt/data/tmp/mission-stack-$C-XXXX); trap 'rm -rf $W' EXIT; cd $W; mkdir 16 12
log() { echo "[$(date -Is)] chr$C $*"; }
fail=0
b() { local s=$(date +%s); $FV sa-build --no-progress --assembly GRCh38 "$@" > $W/b.$BASHPID.log 2>&1 || { cat $W/b.$BASHPID.log; echo FAILED "$@"; return 1; }; log "sa-build $2 ${*: -1} $(( $(date +%s) - s ))s"; }
both() {  # source input name [extra...]: build at 16 and 12 bits
  local src=$1 in=$2 name=$3; shift 3
  b --source $src -i $in -o $W/16/$name "$@" --chunk-bits 16 && b --source $src -i $in -o $W/12/$name "$@" --chunk-bits 12
}
dbsnp() {
  local n=$(case $C in X) echo 23;; Y) echo 24;; *) echo $C;; esac)
  local acc=$(bcftools view -h $S/GCF_000001405.40.gz | grep -o "NC_0000$(printf %02d $n)\.[0-9]*" | head -1)
  [ -n "$acc" ] || { echo "FAILED no dbSNP accession for chr$C"; return 1; }
  bcftools view -r $acc $S/GCF_000001405.40.gz -Oz -o $W/dbsnp.vcf.gz --threads 2 && both dbsnp $W/dbsnp.vcf.gz dbsnp_chr$C
}
save() {  # the shared disk sometimes answers "no such directory" on first use
  for try in 1 2 3 4 5; do
    ls $OUT/sa16 $OUT/sa12 > /dev/null 2>&1; cp $W/16/* $OUT/sa16/ && cp $W/12/* $OUT/sa12/ && return 0; sleep 20
  done; return 1
}
if [ -n "$ONLY" ]; then
  dbsnp && save && log "dbSNP added" && exit 0
  echo "chr$C dbSNP incomplete"; exit 1
fi
if [ $C = 1 ]; then V=$S/gnomad.exomes.v4.1.sites.chr1.vcf.bgz; else
  V=$W/gnomad.vcf.bgz; G=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr$C.vcf.bgz
  curl -fsSL --retry 10 --retry-delay 30 -o $V $G && curl -fsSL --retry 10 -o $V.tbi $G.tbi || { echo "download failed"; exit 1; }
  log "downloaded $(du -h $V | cut -f1)"
fi
( b --source gnomad -i $V -o $W/16/gnomad_chr$C --chunk-bits 16 ) & p1=$!
( b --source gnomad -i $V -o $W/12/gnomad_chr$C --chunk-bits 12 ) & p2=$!
( python3 $A/scripts/extract_gnomad_scores.py --chrom $C --url $V --out-dir $W > $W/extract.log 2>&1 &&
  both spliceai $W/spliceai_chr$C.vcf spliceai_chr$C && both phylop $W/phylop_chr$C.tsv phylop_chr$C ) & p3=$!
( dbsnp || echo FAILED dbsnp
  python3 -c 'import shutil,sys,zipfile; shutil.copyfileobj(zipfile.ZipFile(sys.argv[1]).open("revel_with_transcript_ids"), sys.stdout.buffer, 1 << 24)' \
    $S/revel-v1.3_all_chromosomes.zip | awk -F, -v c=$C 'NR==1 || $1==c' > $W/revel.csv
  [ $(wc -l < $W/revel.csv) -le 1 ] || both revel $W/revel.csv revel_chr$C
  zcat /results/veprs05/data/AlphaMissense_hg38.tsv.gz | awk -F'\t' -v c=chr$C '/^#/ || $1==c' | gzip > $W/am.tsv.gz
  [ $(zcat $W/am.tsv.gz | grep -vc '^#') = 0 ] || both alphamissense $W/am.tsv.gz alphamissense_chr$C
  $PY - $FA $C <<'PYE' | bgzip -@2 -l1 > $W/mock.vcf.gz
import sys
fa, c = sys.argv[1:]
for line in open(fa + ".fai"):
    f = line.split("\t")
    if f[0] == c:
        n, off, bases, width = map(int, f[1:5]); break
out = sys.stdout
out.write("##fileformat=VCFv4.2\n##contig=<ID=chr%s,length=%d>\n" % (c, n))
out.write('##INFO=<ID=P,Number=1,Type=Float,Description="mock prediction">\n##INFO=<ID=Q,Number=1,Type=Float,Description="mock confidence">\n')
out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
with open(fa, "rb") as s:
    for pos in range(1, n + 1, 93):
        s.seek(off + (pos - 1) // bases * width + (pos - 1) % bases)
        r = s.read(1).decode().upper()
        if r not in "ACGT":
            continue
        for k, a in enumerate("ACGT"):
            if a != r:
                h = (pos * 48271 + k * 40503) % 1000003
                out.write("chr%s\t%d\t.\t%s\t%s\t.\t.\tP=%.3f;Q=%.2f\n" % (c, pos, r, a, (h % 991) / 991, (h % 97) / 97))
PYE
  both custom_vcf $W/mock.vcf.gz tracka_chr$C --name TrackA && both custom_vcf $W/mock.vcf.gz trackb_chr$C --name TrackB ) & p4=$!
for p in $p1 $p2 $p3 $p4; do wait $p || fail=1; done
for f in gnomad spliceai phylop dbsnp tracka trackb; do
  [ -s $W/16/${f}_chr$C.osa2 ] && [ -s $W/12/${f}_chr$C.osa2 ] || { echo "missing $f"; fail=1; }
done
if [ $C = 1 ]; then
  cp $A/sa/dense/*.oga $A/sa/dense/*.osi $W/16/; cp $A/sa/dense/*.oga $A/sa/dense/*.osi $W/12/
  b --source clinvar -i $S/clinvar.vcf.gz -o $W/16/clinvar || fail=1
  b --source clinvar -i $S/clinvar.vcf.gz -o $W/12/clinvar --chunk-bits 12 || fail=1
fi
[ $fail = 0 ] || { echo "chr$C incomplete"; exit 1; }
save && touch $OUT/done/chr$C || { echo "chr$C not saved"; exit 1; }
log "done: $(ls $W/16 | wc -l) files at 65 kb, $(ls $W/12 | wc -l) narrow"
