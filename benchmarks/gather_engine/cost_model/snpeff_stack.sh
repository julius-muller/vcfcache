# Sourced by wgs_mission_snpeff.sh: SnpEff 5.4c (database GRCh38.115, Ensembl 115) followed by
# SnpSift annotate with the eight VCF sources that the Ensembl VEP 116.3 mission runs as --custom
# (vep116_stack.sh): ClinVar, gnomAD exomes, dbSNP, SpliceAI-style SNV and indel scores, PhyloP and
# the two mock tracks. One streaming process per stage, as the SnpEff and SnpSift documentation
# shows. The sources are plain VCFs with the index SnpSift builds for them (.sidx, built once by
# snpeff_sources.sh); through tabix indexes SnpSift is 15-20 times slower, same output
# (snpsift_index_probe.sh). Heaps (whole genome of KPGP-00128, 4.6 M records):
# SnpSift keeps the database records it has read, so a stage's memory grows with the input up to its
# heap. The smallest heaps that complete a whole genome are 2 GB per stage, 4 GB for gnomAD and 5 GB
# for SnpEff (snpsift_wg_heap_probe.sh); dbSNP group 1 needs 3 GB for some inputs (part new-0, jobs
# 3710/3711). With 6/5/3/2 GB and the JVM's default collector the 13 processes fill the 30 GB of a
# node and SnpEff is killed (snpeff_pipe_probe.sh, job 3718, 29.2 GB); with the serial collector,
# which returns memory the heap no longer needs (free ratio 10-30%), the SnpSift stages leave room:
# 27.2 GB peak, 3.5 GB of the node free (job 3719). Collector and heap do not change the output.
# Time on chr1-4 of KPGP-00128, where the default collector fits too (snpeff_gc_check.sh, sl-w4,
# jobs 3721/3722): default 1864-2158 s (4 runs), serial 1909-1969, G1 returning memory
# (-XX:G1PeriodicGCInterval=10000 and the same free ratios; whole genome 25.2 GB, job 3723)
# 1935-1978: the collectors that let the whole genome fit cost about 2-3% against the default one,
# within its spread; G1 returning memory is no faster than serial. SnpSift
# needs the SnpEff config to start. -noId keeps the sample's ID column (the sources carry no IDs).
# VEP's plugin sources (REVEL, AlphaMissense, pLI, LoFtool, GWAS) have no SnpSift counterpart for
# these files and are left out.
# Needs: CU (custom VCF dir), T (SpliceAI-style VCF dir), ERR (stderr file).
SE=/results/snpeff/dist
# No JVM perf-data file, and JVM warnings to stderr: the JVM writes its warnings to stdout by default,
# and when the perf-data file of a PID was locked ("hsperfdata ... locked by another process", seen on
# the cluster) the warning landed in the VCF stream between two stages (snpeff_dbsnp_split.sh, job 3663)
JAVA="$SE/jre/bin/java -XX:-UsePerfData -Xlog:disable -Xlog:all=warning:stderr"
SNPEFF="$JAVA -Xmx${SNPEFF_HEAP:-6g} -jar $SE/snpEff/snpEff.jar -noLog -noStats -nodownload -c $SE/snpEff/snpEff.config -dataDir $SE/snpEff/data GRCh38.115 - 2>> $ERR"
sift() { echo "| $JAVA ${SIFT_JOPT--XX:+UseSerialGC -XX:MinHeapFreeRatio=10 -XX:MaxHeapFreeRatio=30} -Xmx${3:-2g} -jar $SE/SnpSift.jar annotate -noLog -c $SE/snpEff/snpEff.config -noId ${2:+-name $2} $1 - 2>> $ERR"; }
# dbSNP: five chromosome groups with their own .sidx (snpeff_dbsnp_split.sh: the index of the whole
# file does not fit a 26 GB heap), else the whole file through tabix
if [ -e $CU/dbsnp_g1.vcf.sidx ]; then DBSNP=$(for g in 1 2 3 4 5; do echo -n "$(sift $CU/dbsnp_g$g.vcf "" ${DBSNP_HEAP:-3g}) "; done)
else DBSNP=$(sift $CU/dbsnp.vcf.gz); fi
SPIPE="$SNPEFF $(sift $CU/clinvar.vcf) $(sift $CU/gnomad.vcf "" ${GNOMAD_HEAP:-5g}) $DBSNP"
SPIPE="$SPIPE $(sift $T/spliceai_snv.vcf SNV_) $(sift $T/spliceai_indel.vcf INDEL_) $(sift $CU/phylop.vcf)"
SPIPE="$SPIPE $(sift $CU/trackA.vcf) $(sift $CU/trackB.vcf)"
SNPEFF_VERSION="$JAVA -jar $SE/snpEff/snpEff.jar -version"
