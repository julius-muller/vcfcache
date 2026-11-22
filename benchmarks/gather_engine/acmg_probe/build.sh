#!/bin/bash
# Build the fastVEP 0.5.0 annotation stacks for the chr1 probe, a chr1 cache
# for each recipe tier, then submit the probe units. Runs as a Slurm job.
#
# Tiers (cumulative): A = today's clinical recipe (ClinVar + gnomAD exomes AF);
# B0 = fastVEP's nine-source ACMG stack; B = B0 with --acmg; C = B plus the
# dense sources dbSNP and AlphaMissense.
set -euo pipefail
DEP=$1; COMMIT=$2
D=/results/fastvep05/acmg; S=$D/sources
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
export PATH=/results/input-preserving-20261105-tools/bin:$PATH   # tabix, bgzip
test -e /mnt/data/slurm-results/fastvep05/acmg/sources/DOWNLOADED 2>/dev/null || test -e $S/DOWNLOADED
STACK=$D/sa/stack; DENSE=$D/sa/dense; mkdir -p $STACK $DENSE
W=$(mktemp -d /mnt/data/tmp/acmg-build-XXXX); trap 'rm -rf $W' EXIT
log() { echo "[$(date -Is)] $*"; }
b() { local s=$(date +%s); $FV sa-build --no-progress --assembly GRCh38 "$@"; log "sa-build $* $(( $(date +%s) - s ))s"; }
$FV --version
if [ ! -e $STACK/DONE ]; then
  b --source clinvar -i $S/clinvar.vcf.gz -o $STACK/clinvar
  b --source clinvar_protein -i $S/variant_summary.txt.gz -o $STACK/clinvar_protein
  b --source gnomad_genes -i $S/gnomad.v4.1.constraint_metrics.tsv -o $STACK/gnomad_genes
  python3 $D/scripts/clingen_gdv_to_oga.py $S/clingen_gene_validity.csv $W/clingen_gdv.tsv
  b --source omim -i $W/clingen_gdv.tsv -o $STACK/omim
  # The probe annotates chr1 only, so per-chromosome sources are built for chr1.
  # Workers have no unzip; stream the single member with Python.
  python3 -c 'import shutil,sys,zipfile; shutil.copyfileobj(zipfile.ZipFile(sys.argv[1]).open("revel_with_transcript_ids"), sys.stdout.buffer, 1 << 24)' \
    $S/revel-v1.3_all_chromosomes.zip | awk -F, 'NR==1 || $1=="1"' > $W/revel_chr1.csv
  b --source revel -i $W/revel_chr1.csv -o $STACK/revel_chr1
  python3 $D/scripts/repeatmasker_to_bed.py $S/rmsk.txt.gz > $W/repeatmasker.bed
  b --source custom_bed --name repeatmasker -i $W/repeatmasker.bed -o $STACK/repeatmasker
  b --source gnomad -i $S/gnomad.exomes.v4.1.sites.chr1.vcf.bgz -o $STACK/gnomad_chr1 --chunk-bits 16
  python3 $D/scripts/extract_gnomad_scores.py --chrom 1 --url $S/gnomad.exomes.v4.1.sites.chr1.vcf.bgz --out-dir $W
  b --source spliceai -i $W/spliceai_chr1.vcf -o $STACK/spliceai_chr1 --chunk-bits 16
  b --source phylop -i $W/phylop_chr1.tsv -o $STACK/phylop_chr1 --chunk-bits 16
  rm -f $W/spliceai_chr1.vcf $W/phylop_chr1.tsv $W/revel_chr1.csv
  touch $STACK/DONE
fi
if [ ! -e $DENSE/DONE ]; then
  for f in $STACK/*; do [ "$(basename $f)" = DONE ] || ln -f $f $DENSE/; done
  bcftools view -r NC_000001.11 $S/GCF_000001405.40.gz -Oz -o $W/dbsnp_chr1.vcf.gz --threads 8
  b --source dbsnp -i $W/dbsnp_chr1.vcf.gz -o $DENSE/dbsnp_chr1 --chunk-bits 16
  zcat /results/veprs05/data/AlphaMissense_hg38.tsv.gz | awk -F'\t' '/^#/ || $1=="chr1"' | gzip > $W/am_chr1.tsv.gz
  b --source alphamissense -i $W/am_chr1.tsv.gz -o $DENSE/alphamissense_chr1
  rm -f $W/dbsnp_chr1.vcf.gz $W/am_chr1.tsv.gz
  touch $DENSE/DONE
fi
du -sh $STACK $DENSE; ls -la $STACK $DENSE
python3 $D/scripts/check_acmg_stack.py --fastvep $FV --gff3 $FVD/data/GRCh38/Homo_sapiens.GRCh38.116.gff3 \
  --fasta $FA --sa-dir $STACK || log "stack check reported gaps (expected off chr1)"
# One chr1 cache per tier, from the chr1 part of the same gnomAD blueprint.
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
bcftools view -r chr1 $BP -Ob -W -o $W/bp_chr1.bcf --threads 8
OLD=/results/fastvep05/caches/fastvep050_clinical_grch38_af001/cache/fastvep050_clinical_grch38_af001
for tier in A B0 B C; do
  case $tier in
    A) SA=/results/richness/sa/clinical_GRCh38; EXTRA="";;
    B0) SA=$STACK; EXTRA="";;
    B) SA=$STACK; EXTRA="--acmg";;
    C) SA=$DENSE; EXTRA="--acmg";;
  esac
  NAME=fastvep050_${tier}_chr1; ROOT=$D/caches/$NAME; C=$ROOT/cache/$NAME
  [ -e $ROOT/COMPLETE ] && continue
  rm -rf $ROOT; mkdir -p $C $ROOT/blueprint; cp -r $(dirname $(dirname $OLD))/workflow $ROOT/
  sed -e "s#--sa-dir /results/richness/sa/clinical_GRCh38#--sa-dir $SA${EXTRA:+ $EXTRA}#" $OLD/annotation.yaml > $C/annotation.yaml
  grep -q -- "--sa-dir $SA${EXTRA:+ $EXTRA}" $C/annotation.yaml
  cp $OLD/params.snapshot.yaml $C/
  start=$(date +%s)
  bcftools view $W/bp_chr1.bcf -Ov | env RAYON_NUM_THREADS=16 $FV annotate --input - --output - --output-format vcf \
    --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $SA $EXTRA 2> $W/err_$tier.txt | \
    bcftools view -Ob -o $C/vcfcache_annotated.bcf -W --threads 8
  in=$(bcftools index -n $W/bp_chr1.bcf); out=$(bcftools index -n $C/vcfcache_annotated.bcf)
  log "cache $tier: $in -> $out records, $(( $(date +%s) - start ))s, $(du -h $C/vcfcache_annotated.bcf | cut -f1)"
  tail -5 $W/err_$tier.txt
  test "$in" = "$out"
  bcftools view -h $C/vcfcache_annotated.bcf | grep '^##INFO=<ID=\(CSQ\|FV_\|SpliceAI\)' | cut -c1-100
  echo "{\"tier\": \"$tier\", \"sa_dir\": \"$SA\", \"extra\": \"$EXTRA\", \"records\": $out, \"seconds\": $(( $(date +%s) - start )), \"bytes\": $(stat -c %s $C/vcfcache_annotated.bcf)}" > $C/build_provenance.json
  touch $ROOT/COMPLETE
done
cd /tmp
env PYTHONPATH=$DEP /results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python \
  -m benchmarks.gather_engine.acmg_probe.plan $D/probe-$(date +%Y%m%d) --deployment $DEP --commit $COMMIT
