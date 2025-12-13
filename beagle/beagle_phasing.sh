#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J beagle_phasing
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 4
# Usage: ./beagle_phasing.sh [--panel 1000gp|ncbn1000gp|ncbn] <input.vcf.gz> <output_dir>

# Parse arguments
PANEL="1000gp"

while [[ $# -gt 0 ]]; do
    case $1 in
        --panel)
            if [ "$2" = "1000gp" ] || [ "$2" = "ncbn1000gp" ] || [ "$2" = "ncbn" ]; then
                PANEL="$2"
            else
                echo "Error: --panel must be '1000gp', 'ncbn1000gp', or 'ncbn'"
                exit 1
            fi
            shift 2
            ;;
        *)
            break
            ;;
    esac
done

INPUT_VCF="$1"
OUTPUT_DIR="$2"

if [ -z "$INPUT_VCF" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 [--panel 1000gp|ncbn1000gp|ncbn] <input.vcf.gz> <output_dir>"
    exit 1
fi

# Convert to absolute paths
INPUT_VCF=$(realpath "$INPUT_VCF")
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

# Check if input VCF exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF not found: $INPUT_VCF"
    exit 1
fi

# Get input VCF basename and construct output VCF name
INPUT_BASENAME=$(basename "$INPUT_VCF" .vcf.gz)
OUTPUT_VCF="${OUTPUT_DIR}/${INPUT_BASENAME}_phased.vcf.gz"

# Tool paths
BEAGLE_JAR="/home/itoyu8/bin/beagle/beagle.17Dec24.224.jar"
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"

# Set reference panel paths based on --panel option
if [ "$PANEL" = "1000gp" ]; then
    REF_PREFIX="/home/itoyu8/database/1000genomes/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_"
elif [ "$PANEL" = "ncbn1000gp" ]; then
    REF_DIR="/home/itoyu8/database/1000genomes/ncbn_window40_Map"
elif [ "$PANEL" = "ncbn" ]; then
    REF_DIR="/home/itoyu8/database/1000genomes/ncbn_pure"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p log

# Create temporary directory for intermediate files
TEMP_DIR="${OUTPUT_DIR}/beagle_temp"
mkdir -p "$TEMP_DIR"

# Define chromosomes to process (chr1-22, chrX)
CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)

# Array to store phased VCF files
PHASED_VCFS=()

# Process each chromosome
for CHR in "${CHROMOSOMES[@]}"; do
    # Determine reference VCF filename based on panel
    if [ "$PANEL" = "1000gp" ]; then
        if [ "$CHR" = "chrX" ]; then
            REF_VCF="${REF_PREFIX}${CHR}.filtered.eagle2-phased.v2.vcf.gz"
        else
            REF_VCF="${REF_PREFIX}${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz"
        fi
    else
        REF_VCF="${REF_DIR}/NCBN_genotype.beagle_${CHR}.vcf.gz"
    fi

    # Output file for this chromosome
    PHASED_VCF_PREFIX="${TEMP_DIR}/phased.${CHR}"
    PHASED_VCF="${PHASED_VCF_PREFIX}.vcf.gz"

    # Run Beagle with error handling
    if java -Xmx24g -jar "$BEAGLE_JAR" \
        gt="$INPUT_VCF" \
        ref="$REF_VCF" \
        out="$PHASED_VCF_PREFIX" \
        chrom="$CHR" \
        nthreads=4; then

        # Index the phased VCF
        if tabix -p vcf "$PHASED_VCF"; then
            # Add to array only if both Beagle and tabix succeeded
            PHASED_VCFS+=("$PHASED_VCF")
        else
            echo "Warning: Failed to index $CHR, skipping"
        fi
    else
        echo "Warning: Beagle failed for $CHR, skipping"
    fi
done

# Check if any chromosomes were successfully phased
if [ ${#PHASED_VCFS[@]} -eq 0 ]; then
    echo "Error: No chromosomes were successfully phased"
    exit 1
fi

# Concatenate all phased VCFs and filter to keep only heterozygous SNPs
# (homozygous sites are not useful for downstream haplotagging)
"$BCFTOOLS" concat "${PHASED_VCFS[@]}" | \
    "$BCFTOOLS" view -v snps -g het -O z -o "$OUTPUT_VCF"

# Index final VCF
tabix -p vcf "$OUTPUT_VCF"
