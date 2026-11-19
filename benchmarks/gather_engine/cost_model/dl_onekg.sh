#!/bin/bash
# 1000 Genomes 30x panel of one chromosome, sites only, kept in /results/sources/onekg so that
# builds at other widths do not download it again. Resumes a broken or stalled transfer.
# FTP, not HTTP: on 2026-11-19 the same server gave 3.1 MB/s by FTP and 0.7 MB/s by HTTP.
# Usage: dl_onekg.sh <chromosome without chr>
set -uo pipefail
C=$1; export PATH=/results/input-preserving-20261105-tools/bin:$PATH
ls /results /results/sources > /dev/null 2>&1
D=/results/sources/onekg; mkdir -p $D; T=$D/chr$C.sites.vcf.gz
[ -s $T ] && exit 0
v=""; [ $C = X ] && v=".v2"
K=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr$C.filtered.SNV_INDEL_SV_phased_panel$v.vcf.gz
P=/mnt/data/tmp/onekg_chr$C.vcf.gz
for try in $(seq 1 30); do
  curl -fsSL --speed-limit 20000 --speed-time 120 -C - -o $P $K && break
  echo "try $try ended at $(du -h $P | cut -f1)"; sleep 30
done
bcftools view -G $P -Oz -o $T.part --threads 4 && mv $T.part $T && rm $P || { echo "FAILED 1000 Genomes chr$C"; exit 1; }
echo "chr$C 1000 Genomes sites: $(du -h $T | cut -f1)"
