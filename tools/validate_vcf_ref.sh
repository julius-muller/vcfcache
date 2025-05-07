
#!/bin/bash
# validate_vcf_ref.sh - BCF/reference compatibility validator with configurable sampling
# Usage: ./validate_vcf_ref.sh <bcf_file> <reference_fasta> <chr_add_file> [num_variants_to_sample]

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
    STEP_START_TIME=$(date +%s.%N)
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

# Function to end step timing
function end_step() {
    local STEP_END_TIME=$(date +%s.%N)
    local STEP_DURATION=$(echo "$STEP_END_TIME - $STEP_START_TIME" | bc)
    echo "  Time: ${STEP_DURATION}s" >> timing.log
}

# Start overall timing
SCRIPT_START_TIME=$(date +%s.%N)
echo "Starting validation script at $(date)" > timing.log

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
end_step

# Extract reference genome information
print_step "Extracting reference genome information..."
REF_STEP_START=$(date +%s.%N)
REF_CHROMS=($(cut -f1 "${REFERENCE}.fai"))
REF_CHROM_COUNT=${#REF_CHROMS[@]}
echo "  Found $REF_CHROM_COUNT chromosomes in reference genome"
REF_STEP_END=$(date +%s.%N)
REF_STEP_DURATION=$(echo "$REF_STEP_END - $REF_STEP_START" | bc)
echo "  Reference extraction time: ${REF_STEP_DURATION}s" >> timing.log
end_step

# Extract VCF chromosome information from index
print_step "Extracting VCF chromosome information..."
VCF_INDEX_START=$(date +%s.%N)
echo "  Retrieved chromosomes from VCF index"
VCF_CHROMS=($(bcftools index -s "$BCF_FILE" 2>/dev/null | cut -f1 | sort -u))
VCF_CHROM_COUNT=${#VCF_CHROMS[@]}
echo "  Found $VCF_CHROM_COUNT unique chromosomes in VCF"
VCF_INDEX_END=$(date +%s.%N)
VCF_INDEX_DURATION=$(echo "$VCF_INDEX_END - $VCF_INDEX_START" | bc)
echo "  VCF index extraction time: ${VCF_INDEX_DURATION}s" >> timing.log
end_step

# Load chromosome mapping
print_step "Processing chromosome mappings..."
MAPPING_START=$(date +%s.%N)
declare -A CHR_MAP
while read -r FROM TO; do
    [[ -z "$FROM" || "$FROM" == \#* || -z "$TO" ]] && continue
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
MAPPING_END=$(date +%s.%N)
MAPPING_DURATION=$(echo "$MAPPING_END - $MAPPING_START" | bc)
echo "  Chromosome mapping time: ${MAPPING_DURATION}s" >> timing.log
end_step

# Select variants for validation
print_step "Selecting variants for reference sequence validation..."
VARIANT_SELECTION_START=$(date +%s.%N)

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

# First try with common chromosomes for finding variants
VARIANT_LINES=""
COMMON_CHROMS=("chr1" "1" "chrX" "X")
for chr in "${COMMON_CHROMS[@]}"; do
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

VARIANT_SELECTION_END=$(date +%s.%N)
VARIANT_SELECTION_DURATION=$(echo "$VARIANT_SELECTION_END - $VARIANT_SELECTION_START" | bc)
echo "  Variant selection time: ${VARIANT_SELECTION_DURATION}s" >> timing.log

# Process each variant for validation
VALIDATED=0

# Create a temporary file to store the parsed variants
TEMP_VARIANTS=$(mktemp)
echo "$VARIANT_LINES" > "$TEMP_VARIANTS"

# Process the variants line by line (avoiding the pipeline subshell issue)
TOTAL_CHECKED=0
VALIDATION_START=$(date +%s.%N)
while IFS= read -r line; do
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
done < "$TEMP_VARIANTS"
VALIDATION_END=$(date +%s.%N)
VALIDATION_DURATION=$(echo "$VALIDATION_END - $VALIDATION_START" | bc)
echo "  Variant validation time: ${VALIDATION_DURATION}s" >> timing.log

# Clean up temp file
rm -f "$TEMP_VARIANTS"

# Calculate total runtime
SCRIPT_END_TIME=$(date +%s.%N)
TOTAL_RUNTIME=$(echo "$SCRIPT_END_TIME - $SCRIPT_START_TIME" | bc)
echo "Total script runtime: ${TOTAL_RUNTIME}s" >> timing.log
echo "Script completed at $(date)" >> timing.log

if [[ $VALIDATED -eq 1 ]]; then
    print_success "Validation complete. BCF is compatible with the reference genome."
    # Print summary of timing information at the end
    echo -e "\n== PERFORMANCE SUMMARY =="
    cat timing.log
    exit 0
else
    print_error "Could not validate any variants against the reference genome. Check chromosome naming and coordinates."
fi