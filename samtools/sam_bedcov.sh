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

# Create BED file with 10kb windows
WINDOWS_BED="${OUTPUT_DIR}/windows.bed"
/home/itoyu8/bin/bedtools/bedtools2/bin/bedtools makewindows -g "${GENOME_FILE}" -w 10000 > "${WINDOWS_BED}"

# Run samtools bedcov
OUTPUT_FILE="${OUTPUT_DIR}/bedcov_results.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools bedcov -c "${WINDOWS_BED}" "${BAM_FILES[@]}" > "${OUTPUT_FILE}"