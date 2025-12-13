#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J count_hpdepth
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 4

# Usage: sbatch count_haplotype_depth.sh --vcf <phased.vcf.gz> --tumor-h1 <tumor.H1.bam> --tumor-h2 <tumor.H2.bam> --normal-h1 <normal.H1.bam> --normal-h2 <normal.H2.bam> --output <output.tsv>

set -e

# Parse arguments
VCF=""
TUMOR_H1=""
TUMOR_H2=""
NORMAL_H1=""
NORMAL_H2=""
OUTPUT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf)
            VCF="$2"
            shift 2
            ;;
        --tumor-h1)
            TUMOR_H1="$2"
            shift 2
            ;;
        --tumor-h2)
            TUMOR_H2="$2"
            shift 2
            ;;
        --normal-h1)
            NORMAL_H1="$2"
            shift 2
            ;;
        --normal-h2)
            NORMAL_H2="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown option $1"
            echo "Usage: $0 --vcf <phased.vcf.gz> --tumor-h1 <tumor.H1.bam> --tumor-h2 <tumor.H2.bam> --normal-h1 <normal.H1.bam> --normal-h2 <normal.H2.bam> --output <output.tsv>"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$VCF" ] || [ -z "$TUMOR_H1" ] || [ -z "$TUMOR_H2" ] || [ -z "$NORMAL_H1" ] || [ -z "$NORMAL_H2" ] || [ -z "$OUTPUT" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --vcf <phased.vcf.gz> --tumor-h1 <tumor.H1.bam> --tumor-h2 <tumor.H2.bam> --normal-h1 <normal.H1.bam> --normal-h2 <normal.H2.bam> --output <output.tsv>"
    exit 1
fi

# Tool paths
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"
SAMTOOLS="/home/itoyu8/bin/samtools/samtools-1.19/samtools"

# Use output path directly as specified by user
OUTPUT_FILE="${OUTPUT}"

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Create temporary directory
TEMP_DIR=$(mktemp -d)
trap "rm -rf ${TEMP_DIR}" EXIT

# Step 1: Extract hetero SNP positions from VCF (filter for chr1-22, chrX, chrY)
# Output format: CHROM POS0 (0-based position, 2-column BED format)
# Filter: SNPs only, phased heterozygous genotypes (including multiallelic sites like 1|2)
echo "Extracting hetero SNP positions from VCF..."
"${BCFTOOLS}" view -v snps -g het -p "${VCF}" | \
    "${BCFTOOLS}" query -f '%CHROM\t%POS0\n' | \
    awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/' \
    > "${TEMP_DIR}/targets.bed"

# Step 2: Run samtools depth for each BAM file
echo "Calculating depth for tumor H1..."
"${SAMTOOLS}" depth -a -b "${TEMP_DIR}/targets.bed" "${TUMOR_H1}" > "${TEMP_DIR}/depth_tumor_h1.txt"

echo "Calculating depth for tumor H2..."
"${SAMTOOLS}" depth -a -b "${TEMP_DIR}/targets.bed" "${TUMOR_H2}" > "${TEMP_DIR}/depth_tumor_h2.txt"

echo "Calculating depth for normal H1..."
"${SAMTOOLS}" depth -a -b "${TEMP_DIR}/targets.bed" "${NORMAL_H1}" > "${TEMP_DIR}/depth_normal_h1.txt"

echo "Calculating depth for normal H2..."
"${SAMTOOLS}" depth -a -b "${TEMP_DIR}/targets.bed" "${NORMAL_H2}" > "${TEMP_DIR}/depth_normal_h2.txt"

# Step 3: Combine all data with index
echo "Combining data..."
paste \
    <(cut -f 1,2 "${TEMP_DIR}/depth_tumor_h1.txt") \
    <(cut -f 3 "${TEMP_DIR}/depth_tumor_h1.txt") \
    <(cut -f 3 "${TEMP_DIR}/depth_tumor_h2.txt") \
    <(cut -f 3 "${TEMP_DIR}/depth_normal_h1.txt") \
    <(cut -f 3 "${TEMP_DIR}/depth_normal_h2.txt") | \
awk 'BEGIN {
    OFS="\t"
    idx=0
    print "index", "chrom", "position", "tumor_hp1", "tumor_hp2", "normal_hp1", "normal_hp2"
}
{
    print idx, $1, $2, $3, $4, $5, $6
    idx++
}' > "${OUTPUT_FILE}"

echo "Output written to: ${OUTPUT_FILE}"
echo "Total SNVs: $(tail -n +2 "${OUTPUT_FILE}" | wc -l)"
