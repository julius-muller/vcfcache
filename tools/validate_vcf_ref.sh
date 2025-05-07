#!/bin/bash
# validate_bcf_v10.sh - BCF/reference compatibility validator with configurable sampling
# Usage: ./validate_bcf_v10.sh <bcf_file> <reference_fasta> <chr_add_file> [num_variants_to_sample]

set -e  # Exit on error

# Colors and formatting
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Symbols
CHECK_MARK="✓"
CROSS_MARK="✗"
ARROW="→"

# Function to print a step
function print_step() {
    echo -e "${BLUE}${ARROW} $1${NC}"
}

# Function to print success
function print_success() {
    echo -e "${GREEN}${CHECK_MARK} $1${NC}"
}

# Function to print error and exit
function print_error() {
    echo -e "${RED}${CROSS_MARK} $1${NC}" >&2
    exit 1
}

# Function to print warning
function print_warning() {
    echo -e "${YELLOW}! $1${NC}" >&2
}

# Check arguments
if [[ $# -lt 3 ]]; then
    print_error "Usage: $0 <bcf_file> <reference_fasta> <chr_add_file> [num_variants_to_sample]"
fi

BCF_FILE="$1"
REFERENCE="$2"
CHR_ADD="$3"
NUM_VARIANTS=${4:-3}  # Default to 3 variants if not specified

print_step "Checking input files..."

# Check if files exist and are readable
if [[ ! -f "$BCF_FILE" ]]; then
    print_error "BCF file not found: $BCF_FILE"
fi

if [[ ! -f "$REFERENCE" ]]; then
    print_error "Reference file not found: $REFERENCE"
fi

if [[ ! -f "$CHR_ADD" ]]; then
    print_error "Chromosome mapping file not found: $CHR_ADD"
fi

# Check for indexes
if [[ -f "${BCF_FILE}.csi" ]]; then
    BCF_INDEX="${BCF_FILE}.csi"
elif [[ -f "${BCF_FILE}.tbi" ]]; then
    BCF_INDEX="${BCF_FILE}.tbi"
else
    print_error "No index found for BCF file. Create it with: bcftools index $BCF_FILE"
fi

if [[ ! -f "${REFERENCE}.fai" ]]; then
    print_error "Reference index not found. Create it with: samtools faidx $REFERENCE"
fi

print_success "All input files exist"
print_success "Found both VCF and reference indexes"

# Extract reference genome information
print_step "Extracting reference genome information..."
REF_CHROMS=($(cut -f1 "${REFERENCE}.fai"))
REF_CHROM_COUNT=${#REF_CHROMS[@]}
echo "  Found $REF_CHROM_COUNT chromosomes in reference genome"

# Extract VCF chromosome information from index
print_step "Extracting VCF chromosome information..."
echo "  Retrieved chromosomes from VCF index"
VCF_CHROMS=($(bcftools index -s "$BCF_FILE" 2>/dev/null | cut -f1 | sort -u))
VCF_CHROM_COUNT=${#VCF_CHROMS[@]}
echo "  Found $VCF_CHROM_COUNT unique chromosomes in VCF"

# Load chromosome mapping
print_step "Processing chromosome mappings..."
declare -A CHR_MAP
while read -r FROM TO; do
    [[ -z "$FROM" || -z "$TO" ]] && continue
    CHR_MAP["$FROM"]="$TO"
done < "$CHR_ADD"
MAPPING_COUNT=${#CHR_MAP[@]}
echo "  Loaded $MAPPING_COUNT chromosome mappings from $CHR_ADD"

# Determine chromosome naming convention
VCF_HAS_CHR_PREFIX=0
REF_HAS_CHR_PREFIX=0

for chrom in "${VCF_CHROMS[@]:0:1}"; do
    if [[ "$chrom" == chr* ]]; then
        VCF_HAS_CHR_PREFIX=1
        break
    fi
done

for chrom in "${REF_CHROMS[@]:0:1}"; do
    if [[ "$chrom" == chr* ]]; then
        REF_HAS_CHR_PREFIX=1
        break
    fi
done

if [[ $VCF_HAS_CHR_PREFIX -eq $REF_HAS_CHR_PREFIX ]]; then
    echo "  Both VCF and reference use the same chromosome naming convention"
else
    if [[ $VCF_HAS_CHR_PREFIX -eq 1 ]]; then
        echo "  VCF uses 'chr' prefix but reference doesn't"
    else
        echo "  Reference uses 'chr' prefix but VCF doesn't"
    fi
fi

# Select variants for validation
print_step "Selecting variants for reference sequence validation..."

# First try with common chromosomes
COMMON_CHROMS=("chr1" "1" "chrX" "X")
for chr in "${COMMON_CHROMS[@]}"; do
    # Skip if already validated
    [[ $VALIDATED -eq 1 ]] && break

    # Try to find variants on this chromosome
    VARIANT_LINES=$(bcftools view -H -r "$chr" "$BCF_FILE" 2>/dev/null | head -n $NUM_VARIANTS)

    if [[ -n "$VARIANT_LINES" ]]; then
        break
    fi
done

# If no variants found in common chromosomes, try any chromosome from VCF
if [[ -z "$VARIANT_LINES" ]]; then
    for chr in "${VCF_CHROMS[@]}"; do
        VARIANT_LINES=$(bcftools view -H -r "$chr" "$BCF_FILE" 2>/dev/null | head -n $NUM_VARIANTS)
        if [[ -n "$VARIANT_LINES" ]]; then
            break
        fi
    done
fi

# If still no variants found, error out
if [[ -z "$VARIANT_LINES" ]]; then
    print_error "Could not extract any variants from the BCF file for validation."
fi

# Function to map chromosome name
function map_chromosome() {
    local chr="$1"
    # Direct mapping from chr_add
    if [[ -n "${CHR_MAP[$chr]}" ]]; then
        echo "${CHR_MAP[$chr]}"
        return 0
    fi

    # Try without 'chr' prefix if it exists
    if [[ "$chr" == chr* && -n "${CHR_MAP[${chr#chr}]}" ]]; then
        echo "${CHR_MAP[${chr#chr}]}"
        return 0
    fi

    # Try with 'chr' prefix if it doesn't exist
    if [[ "$chr" != chr* && -n "${CHR_MAP[chr$chr]}" ]]; then
        echo "${CHR_MAP[chr$chr]}"
        return 0
    fi

    # Handle chromosome naming convention differences
    if [[ $VCF_HAS_CHR_PREFIX -ne $REF_HAS_CHR_PREFIX ]]; then
        if [[ $VCF_HAS_CHR_PREFIX -eq 1 && "$chr" == chr* ]]; then
            echo "${chr#chr}"
            return 0
        elif [[ $VCF_HAS_CHR_PREFIX -eq 0 && "$chr" != chr* ]]; then
            echo "chr$chr"
            return 0
        fi
    fi

    # Return original if no mapping found
    echo "$chr"
}

# Process each variant for validation
VALIDATED=0
TOTAL_CHECKED=0

echo "$VARIANT_LINES" | while read -r line; do
    CHROM=$(echo "$line" | cut -f1)
    POS=$(echo "$line" | cut -f2)
    REF_ALLELE=$(echo "$line" | cut -f4)

    TOTAL_CHECKED=$((TOTAL_CHECKED + 1))

    # Map chromosome name for reference lookup
    REF_CHROM=$(map_chromosome "$CHROM")

    echo "  Checking variant $TOTAL_CHECKED: $CHROM:$POS (REF=$REF_ALLELE)"
    if [[ "$REF_CHROM" != "$CHROM" ]]; then
        echo "    Mapped to $REF_CHROM in reference"
    fi

    # Get reference sequence
    REF_SEQ=$(samtools faidx "$REFERENCE" "$REF_CHROM:$POS-$((POS + ${#REF_ALLELE} - 1))" 2>/dev/null | grep -v ">" | tr -d '\n')

    # If failed, try alternative chromosome names
    if [[ -z "$REF_SEQ" ]]; then
        for alt_chrom in "${CHROM#chr}" "chr${CHROM#chr}"; do
            if [[ "$alt_chrom" != "$REF_CHROM" ]]; then
                echo "    Trying alternative chromosome name: $alt_chrom"
                REF_SEQ=$(samtools faidx "$REFERENCE" "$alt_chrom:$POS-$((POS + ${#REF_ALLELE} - 1))" 2>/dev/null | grep -v ">" | tr -d '\n')
                if [[ -n "$REF_SEQ" ]]; then
                    REF_CHROM="$alt_chrom"
                    break
                fi
            fi
        done
    fi

    if [[ -z "$REF_SEQ" ]]; then
        print_warning "    Could not retrieve reference sequence for $CHROM at position $POS"
        continue
    fi

    # Compare sequences (case-insensitive)
    if [[ "${REF_ALLELE,,}" == "${REF_SEQ,,}" ]]; then
        echo "    Reference sequence match: $REF_SEQ"
        VALIDATED=1
    else
        print_warning "    Reference sequence mismatch. VCF: $REF_ALLELE, Reference: $REF_SEQ"
    fi
done

if [[ $VALIDATED -eq 1 ]]; then
    print_success "Validation complete. BCF is compatible with the reference genome."
    exit 0
else
    print_error "Could not validate any variants against the reference genome. Check chromosome naming and coordinates."
fi