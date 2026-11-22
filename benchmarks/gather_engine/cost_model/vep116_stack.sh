# Sourced by probe_vep116.sh and wgs_mission_vep116.sh: Ensembl VEP 116.3 with the eleven sources of
# the vep-rs mission stack (wgs_mission_veprs.sh). Tabular sources as --custom VCFs (vep116_custom.sh;
# the SpliceAI-style scores too, since Ensembl's SpliceAI plugin reports scores only where the file's
# gene symbol matches the transcript's and these files carry none), the rest as the same plugins.
# Needs: CU (custom VCF dir), T (vep-rs tables dir, for the SpliceAI VCFs), VDIR (writable --dir: the
# GWAS plugin writes its processed file there on first use), FORKS.
VS=/mnt/data/apps/ensembl-vep/116
VEP="apptainer exec -B /mnt/data:/mnt/data -B /results:/results $VS/vep.sif vep"
FA=$VS/cachedir/homo_sapiens/116_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
VD=/results/veprs05/data
cu() { echo "--custom file=$1,short_name=$2,format=vcf,type=exact,coords=0,fields=$3"; }
G=gnomAD_AF; for f in AC AN nhomalt AF_afr AF_amr AF_asj AF_eas AF_fin AF_mid AF_nfe AF_remaining AF_sas; do G="$G%gnomAD_$f"; done
VFLAGS="--offline --cache --dir_cache $VS/cachedir --dir $VDIR --assembly GRCh38 --fork $FORKS --vcf --no_stats --format vcf --hgvs --fasta $FA"
VFLAGS="$VFLAGS $(cu $CU/clinvar.vcf.gz ClinVar CV_CLNSIG%CV_REVSTAT%CV_CLNDN) $(cu $CU/gnomad.vcf.gz gnomADe $G) $(cu $CU/dbsnp.vcf.gz dbSNP dbSNP_RS)"
VFLAGS="$VFLAGS $(cu $T/spliceai_snv.vcf.gz SpliceAI_snv SpliceAI) $(cu $T/spliceai_indel.vcf.gz SpliceAI_indel SpliceAI) $(cu $CU/phylop.vcf.gz PhyloP PhyloP)"
VFLAGS="$VFLAGS $(cu $CU/trackA.vcf.gz TrackA TrackA_P%TrackA_Q) $(cu $CU/trackB.vcf.gz TrackB TrackB_P%TrackB_Q)"
VFLAGS="$VFLAGS --plugin REVEL,file=$VD/new_tabbed_revel_grch38.tsv.gz --plugin AlphaMissense,file=$VD/AlphaMissense_hg38.tsv.gz"
VFLAGS="$VFLAGS --plugin pLI,$VD/pLI_values.txt --plugin LoFtool,$VD/LoFtool_scores.txt --plugin GWAS,file=$VD/gwas-catalog-download-associations-alt-full.tsv"
