#!/bin/bash
# The two mock tracks as ONE genome-wide custom source each. fastVEP treats every
# custom_vcf file as a source of its own (same --name -> fields FV_TRACKA, _2 ... _24),
# so the per-chromosome files of wgs_mission_stack.sh made it look up and print each
# track 24 times. Same generator as there, so the values are identical.
# Writes /results/fastvep05/mission/tracks1/{sa16,sa12}/track{a,b}.osa2
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
ls /results/fastvep05/mission > /dev/null 2>&1
BITS=${BITS:-16 12}   # chunk widths to build
OUT=/results/fastvep05/mission/tracks1; for bits in $BITS; do mkdir -p $OUT/sa$bits; done
W=$(mktemp -d /mnt/data/tmp/tracks1-XXXX); trap 'rm -rf $W' EXIT; cd $W
log() { echo "[$(date -Is)] $*"; }
$PY - $FA <<'PYE' | bgzip -@4 -l1 > mock.vcf.gz
import sys
fa = sys.argv[1]
chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
idx = {}
for line in open(fa + ".fai"):
    f = line.split("\t")
    if f[0] in chroms:
        idx[f[0]] = tuple(map(int, f[1:5]))
out = sys.stdout
out.write("##fileformat=VCFv4.2\n")
for c in chroms:
    out.write("##contig=<ID=chr%s,length=%d>\n" % (c, idx[c][0]))
out.write('##INFO=<ID=P,Number=1,Type=Float,Description="mock prediction">\n##INFO=<ID=Q,Number=1,Type=Float,Description="mock confidence">\n')
out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
with open(fa, "rb") as s:
    for c in chroms:
        n, off, bases, width = idx[c]
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
log "mock track: $(zcat mock.vcf.gz | grep -vc '^#') rows, $(du -h mock.vcf.gz | cut -f1)"
for t in a b; do for bits in $BITS; do
  T=$(echo $t | tr a-z A-Z)
  /usr/bin/time -f "%e s, max RSS %M kB" -o t.txt $FV sa-build --no-progress --assembly GRCh38 --source custom_vcf -i mock.vcf.gz \
    -o $W/track${t}_$bits --name Track$T --chunk-bits $bits > b.log 2>&1 || { tail -5 b.log; cat t.txt; echo "FAILED track$t $bits"; exit 1; }
  log "track$t chunk-bits $bits: $(cat t.txt)"
  for f in $W/track${t}_$bits.*; do cp $f $OUT/sa$bits/track$t.${f#$W/track${t}_$bits.} || exit 1; done
done; done
ls -l $OUT/sa*
echo DONE
