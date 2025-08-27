#!/bin/bash
#SBATCH -J sam_bedcov
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 1

# Usage: sbatch sam_bedcov.sh <BAM1> [BAM2] [BAM3] ...

if [ "$#" -lt 1 ]; then
    echo "Usage: sbatch $0 <BAM1> [BAM2] [BAM3] ..."
    exit 1
fi

BAM_FILES=("$@")
# Set output directory to the same directory as the first BAM file
OUTPUT_DIR=$(dirname "${BAM_FILES[0]}")

# Fixed paths
GENOME_FILE="/home/itoyu8/database/reference/hg38/v0/human.hg38.genome"

# Create log directory
mkdir -p log

# Generate genome file if needed
REF_FASTA_FAI="/home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
if [ ! -f "${GENOME_FILE}" ]; then
    awk -v OFS='\t' '{print $1, $2}' "${REF_FASTA_FAI}" > "${GENOME_FILE}"
fi

# Define the list of chromosomes to process (same as haplotag_pipeline_hifi_allchrom.sh)
CHROMS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

# Create BED file with 10kb windows
WINDOWS_BED="${OUTPUT_DIR}/windows.bed"

# Generate 10kb windows for all chromosomes, then filter to only include target chromosomes
/home/itoyu8/bin/bedtools/bedtools2/bin/bedtools makewindows -g "${GENOME_FILE}" -w 10000 > "${OUTPUT_DIR}/all_genome_windows.bed"

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
OUTPUT_FILE="${OUTPUT_DIR}/bedcov_results.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools bedcov -c "${WINDOWS_BED}" "${BAM_FILES[@]}" > "${OUTPUT_FILE}"