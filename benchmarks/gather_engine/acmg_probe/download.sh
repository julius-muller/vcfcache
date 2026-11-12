#!/bin/bash
# Fetch the public sources of fastVEP's ACMG stack (plus dbSNP and the chr1
# gnomAD exomes file), then submit the build-and-probe job. Runs on the
# controller as appuser: `nohup bash download.sh <deployment> <commit> &`.
set -euo pipefail
DEP=$1; COMMIT=$2
S=/mnt/data/slurm-results/fastvep05/acmg/sources; mkdir -p $S; cd $S
log() { echo "[$(date -Is)] $*"; }
get() {  # url [name]
  local url=$1 name=${2:-$(basename "$1")}
  [ -s "$name" ] && { log "have $name"; return; }
  for try in 1 2 3 4 5 6; do
    curl -fL --retry 10 --retry-delay 30 -C - -o "$name.part" "$url" && break || log "retry $try $name"
  done
  mv "$name.part" "$name"; log "got $name $(du -h "$name" | cut -f1)"
}
G=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1
N=https://ftp.ncbi.nlm.nih.gov
# The two large files in parallel with the small ones.
( get $G/vcf/exomes/gnomad.exomes.v4.1.sites.chr1.vcf.bgz; get $G/vcf/exomes/gnomad.exomes.v4.1.sites.chr1.vcf.bgz.tbi ) &
( get https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz; get https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi; get https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.md5 ) &
get $N/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz; get $N/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5; get $N/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
get $N/pub/clinvar/tab_delimited/variant_summary.txt.gz
get $G/constraint/gnomad.v4.1.constraint_metrics.tsv
get https://search.clinicalgenome.org/kb/gene-validity/download clingen_gene_validity.csv
get https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
get https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
wait
# Integrity: gzip streams complete, published checksums match.
gzip -t clinvar.vcf.gz variant_summary.txt.gz rmsk.txt.gz
[ "$(md5sum clinvar.vcf.gz | cut -d' ' -f1)" = "$(cut -d' ' -f1 clinvar.vcf.gz.md5)" ]
[ "$(md5sum GCF_000001405.40.gz | cut -d' ' -f1)" = "$(cut -d' ' -f1 GCF_000001405.40.gz.md5)" ]
# Header only: grep -q closing the pipe early would fail under pipefail.
bcftools view -h clinvar.vcf.gz | grep '^##INFO=<ID=AF_EXAC' > /dev/null
unzip -tq revel-v1.3_all_chromosomes.zip
sha256sum * > SHA256SUMS
touch DOWNLOADED; log "downloads complete"
cd /tmp
job=$(sbatch --parsable --partition=benchmark --exclusive --nodes=1 -c16 --mem=28000M --time=08:00:00 \
  --chdir=/tmp --nice=0 --job-name=acmg-build \
  --output=/results/fastvep05/acmg/build-%j.out \
  --wrap "bash $DEP/benchmarks/gather_engine/acmg_probe/build.sh $DEP $COMMIT")
log "build job $job"
