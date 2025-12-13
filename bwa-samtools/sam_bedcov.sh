#!/bin/bash
#SBATCH -J sam_bedcov
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 1
# Usage: sbatch sam_bedcov.sh [--reference hg38|chm13] [--window-size SIZE] [-o OUTPUT_DIR] <BAM1> [BAM2] [BAM3] ...

# Parse arguments
REFERENCE_TYPE="hg38"
WINDOW_SIZE=10000
OUTPUT_DIR=""
BAM_FILES=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --reference)
            if [ "$2" = "chm13" ]; then
                REFERENCE_TYPE="chm13"
            elif [ "$2" = "hg38" ]; then
                REFERENCE_TYPE="hg38"
            else
                echo "Error: --reference must be 'hg38' or 'chm13'"
                exit 1
            fi
            shift 2
            ;;
        --window-size)
            WINDOW_SIZE="$2"
            shift 2
            ;;
        -o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            BAM_FILES+=("$1")
            shift
            ;;
    esac
done

if [ ${#BAM_FILES[@]} -lt 1 ]; then
    echo "Usage: sbatch $0 [--reference hg38|chm13] [--window-size SIZE] [-o OUTPUT_DIR] <BAM1> [BAM2] [BAM3] ..."
    exit 1
fi

# Set output directory: use specified dir or same directory as the first BAM file
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR=$(dirname "${BAM_FILES[0]}")
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Set reference-specific genome file paths
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    GENOME_FILE="/home/itoyu8/database/reference/chm13/v2.0/human.chm13.genome"
else
    GENOME_FILE="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/human.hg38.genome"
fi

# Create log directory
mkdir -p log

# Generate genome file if needed (commented out since files are manually created)
# REF_FASTA_FAI path would go here; use: awk -v OFS='\t' '{print $1, $2}' "${REF_FASTA_FAI}" > "${GENOME_FILE}"

# Define the list of chromosomes to process (same as haplotag_pipeline_hifi_allchrom.sh)
CHROMS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

# Create BED file with windows
WINDOWS_BED="${OUTPUT_DIR}/windows.bed"

# Generate windows for all chromosomes, then filter to only include target chromosomes
/home/itoyu8/bin/bedtools/bedtools2/bin/bedtools makewindows -g "${GENOME_FILE}" -w "${WINDOW_SIZE}" > "${OUTPUT_DIR}/all_genome_windows.bed"

# Filter to only include target chromosomes in the specified order
> "${WINDOWS_BED}" # Create empty file
for CHR_NAME in "${CHROMS[@]}"; do
    grep "^${CHR_NAME}\s" "${OUTPUT_DIR}/all_genome_windows.bed" >> "${WINDOWS_BED}" || true
done

# Clean up temporary file
rm -f "${OUTPUT_DIR}/all_genome_windows.bed"

if [ ! -s "${WINDOWS_BED}" ]; then
    echo "Error: Windows file is empty after filtering for target chromosomes."
    exit 1
fi

# Run samtools bedcov
OUTPUT_FILE="${OUTPUT_DIR}/bedcov.tsv"
/home/itoyu8/bin/samtools/samtools-1.19/samtools bedcov -c "${WINDOWS_BED}" "${BAM_FILES[@]}" > "${OUTPUT_FILE}"