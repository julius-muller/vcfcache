#!/bin/bash
# Which input formats do fastVEP 0.5.0's onekg, gerp, gnomad (genomes), dbnsfp
# and custom_bed parsers accept, and can two gnomAD sources share one --sa-dir?
# Small chr22 slices only; needs dl_sources.sh's GERP file and tool.
# Usage: probe_sources.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
ls /results /results/sources > /dev/null 2>&1
TOOLS=/results/fastvep05/tools; SRC=/results/sources
W=$(mktemp -d /mnt/data/tmp/probe-XXXX); trap 'rm -rf $W' EXIT; cd $W
REG=chr22:20000000-21000000
try() { echo "--- sa-build $*"; $FV sa-build --no-progress --assembly GRCh38 "$@" 2>&1 | tail -4; ls -la ${@: -1}* 2>/dev/null | head -3; }
mkdir sa
# gnomAD genomes slice, next to the exomes slice
G=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf
bcftools view -r $REG $G/genomes/gnomad.genomes.v4.1.sites.chr22.vcf.bgz -Oz -o gg.vcf.gz && echo "genomes records $(bcftools view -H gg.vcf.gz | wc -l)"
bcftools view -r $REG $G/exomes/gnomad.exomes.v4.1.sites.chr22.vcf.bgz -Oz -o ge.vcf.gz && echo "exomes records $(bcftools view -H ge.vcf.gz | wc -l)"
try --source gnomad -i ge.vcf.gz --chunk-bits 16 -o sa/gnomad_chr22
try --source gnomad -i gg.vcf.gz --chunk-bits 16 -o sa/gnomad_genomes_chr22
# 1000 Genomes 30x panel slice, sites only
K=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
bcftools view -G -r $REG $K -Oz -o kg.vcf.gz; bcftools view -h kg.vcf.gz | grep "^##INFO" | cut -c1-110 | head -40; bcftools view -H kg.vcf.gz | head -2 | cut -c1-400
try --source onekg -i kg.vcf.gz -o sa/onekg_chr22
# the same with superpopulation AFs named as in the phase-3 release (AFR_AF, ...), if the panel lacks them
bcftools view -h kg.vcf.gz | grep -q "ID=AFR_AF" || {
  bcftools view -h kg.vcf.gz | grep -o "ID=AF_[A-Z]*[^,]*" | head; }
# GERP: bigWig -> bedGraph -> a few text layouts
GB=$SRC/gerp/gerp_91_mammals.GRCh38.bw
$TOOLS/bigWigToBedGraph -chrom=22 -start=20000000 -end=21000000 $GB gerp.bg; head -3 gerp.bg; wc -l gerp.bg
awk 'BEGIN{OFS="\t"} {for(p=$2+1;p<=$3;p++) print "chr"$1,p,$4}' gerp.bg > gerp.tsv; gzip -k gerp.tsv
awk 'BEGIN{OFS="\t"} {print "chr"$1,$2,$3,$4}' gerp.bg > gerp.chr.bg
{ echo "fixedStep chrom=chr22 start=20000001 step=1"; cut -f3 gerp.tsv; } > gerp.wig
try --source gerp -i gerp.tsv.gz -o sa/gerp_tsv_chr22
try --source gerp -i gerp.chr.bg -o sa/gerp_bg_chr22
try --source gerp -i gerp.wig -o sa/gerp_wig_chr22
# ENCODE cCREs as a custom BED
head -2 $SRC/encode/GRCh38-cCREs.bed; wc -l < $SRC/encode/GRCh38-cCREs.bed
try --source custom_bed --name ENCODE_cCRE -i $SRC/encode/GRCh38-cCREs.bed -o sa/ccre
# dbNSFP 5.4a slice (URL from the private file, never echoed)
DBN=$(cat $SRC/dbnsfp/URL_GRCH38)
tabix -H $DBN | tr '\t' '\n' | head -12 | tr '\n' ' '; echo; echo "dbNSFP columns $(tabix -H $DBN | tr '\t' '\n' | wc -l)"
{ tabix -H $DBN; tabix $DBN 22:20000000-21000000; } | bgzip > dbn.tsv.gz; echo "dbNSFP rows $(zcat dbn.tsv.gz | wc -l)"
zcat dbn.tsv.gz | sed -n 2p | cut -c1-300
try --source dbnsfp -i dbn.tsv.gz -o sa/dbnsfp_chr22
# annotate 3000 gnomAD-genome sites of the slice from all databases that built
bcftools view -H gg.vcf.gz | awk 'NR%20==0' | cut -f1-5 | head -3000 > v.body
{ echo "##fileformat=VCFv4.2"; echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; awk 'BEGIN{OFS="\t"}{print $1,$2,".",$4,$5,".",".","."}' v.body; } > v.vcf
ls sa
$FV annotate --input v.vcf --output o.vcf --output-format vcf --transcript-cache $TC --fasta $FA --no-progress --sa-dir sa 2> ann.err; tail -3 ann.err
grep "^##INFO" o.vcf | cut -c1-200
grep -v "^#" o.vcf | head -2 | cut -f8 | tr ';' '\n' | grep -v "^CSQ" | cut -c1-300
for k in $(grep "^##INFO" o.vcf | grep -o "ID=[^,]*" | cut -d= -f2); do echo "$k $(grep -v '^#' o.vcf | grep -c "[;	]$k=")"; done
echo DONE
