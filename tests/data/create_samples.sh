# Extract header lines (which start with "#")
grep '^#' input.vcf > header.vcf

# Randomly sample 1 million variant lines (non-header)
grep -v '^#' input.vcf | shuf -n 1000000 > s1_1m.vcf

# Combine header and sampled variant lines into one VCF
cat header.vcf s1_1m.vcf | bcftools sort -m 5G -Ob --write-index -o s1_1m.bcf
