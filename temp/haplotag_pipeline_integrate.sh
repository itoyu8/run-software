#!/bin/bash
#SBATCH -J haplotag_pipeline_integrate
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 4

set -e

# haplotag_pipeline_integrate.sh
# Integrates existing chromosome-specific outputs into concatenated files
# This script performs Step 5 from haplotag_pipeline_allchrom.sh
# 
# Usage: ./haplotag_pipeline_integrate.sh <output_dir>
# 
# Input files expected in output_dir:
#   - chr1/chr1_readcount.txt through chr22/chr22_readcount.txt
#   - chr1/phased.chr1.vcf.gz through chr22/phased.chr22.vcf.gz
#
# Output files created in output_dir:
#   - all_chroms_readcount.txt (concatenated readcount files)
#   - phased.snp.vcf.gz (concatenated and indexed VCF files)

if [ $# -ne 1 ]; then
    echo "Usage: $0 <output_dir>"
    echo "Example: $0 output"
    exit 1
fi

OUTPUT_DIR="$1"

# Check if output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Error: Output directory $OUTPUT_DIR does not exist"
    exit 1
fi

# Tool paths
BCFTOOLS_PATH="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"
APPTAINER_SIF_PATH="/home/itoyu8/singularity/compat_parabricks-0.2.2.sif"

# Check if tools exist
if [ ! -x "$BCFTOOLS_PATH" ]; then
    echo "Error: bcftools not found or not executable at $BCFTOOLS_PATH"
    exit 1
fi

if [ ! -f "$APPTAINER_SIF_PATH" ]; then
    echo "Error: Apptainer SIF not found at $APPTAINER_SIF_PATH"
    exit 1
fi

echo "--- haplotag_pipeline_integrate.sh Configuration ---"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Bcftools Path: ${BCFTOOLS_PATH}"
echo "Apptainer SIF: ${APPTAINER_SIF_PATH}"
echo "Processing chromosomes: chr1-chr22"
echo "-----------------------------------------------"

# Initialize arrays to store file paths
PHASED_VCF_FILES=()
READCOUNT_FILES=()

# Check and collect chromosome files
echo "Checking for chromosome-specific files..."
for CHR_NUM in {1..22}; do
    CHR_DIR="${OUTPUT_DIR}/chr${CHR_NUM}"
    READCOUNT_FILE="${CHR_DIR}/chr${CHR_NUM}_readcount.txt"
    PHASED_VCF_FILE="${CHR_DIR}/phased.chr${CHR_NUM}.vcf.gz"
    
    echo "Checking chr${CHR_NUM}..."
    
    # Check readcount file
    if [ ! -f "$READCOUNT_FILE" ]; then
        echo "Error: Readcount file not found: $READCOUNT_FILE"
        exit 1
    fi
    echo "  Found readcount file: $READCOUNT_FILE"
    READCOUNT_FILES+=("$READCOUNT_FILE")
    
    # Check phased VCF file
    if [ ! -f "$PHASED_VCF_FILE" ]; then
        echo "Error: Phased VCF file not found: $PHASED_VCF_FILE"
        exit 1
    fi
    echo "  Found phased VCF file: $PHASED_VCF_FILE"
    PHASED_VCF_FILES+=("$PHASED_VCF_FILE")
done

echo "All chromosome files found. Proceeding with integration..."

# --- Step 5: Concatenate all chromosome files ---
echo "Step 5: Concatenating all chromosome files..."

# Concatenate all readcount files
echo "Concatenating readcount files..."
ALL_CHROMS_READCOUNT="${OUTPUT_DIR}/all_chroms_readcount.txt"
cat "${READCOUNT_FILES[@]}" > "${ALL_CHROMS_READCOUNT}" \
  || { echo "Failed to concatenate readcount files"; exit 1; }
echo "Combined readcount file created: ${ALL_CHROMS_READCOUNT}"

# Concatenate all phased VCF files
echo "Concatenating phased VCF files..."
FINAL_PHASED_VCF="${OUTPUT_DIR}/phased.snp.vcf.gz"

# Use bcftools concat to properly merge VCF files
"${BCFTOOLS_PATH}" concat -O z -o "${FINAL_PHASED_VCF}" "${PHASED_VCF_FILES[@]}" \
  || { echo "Failed to concatenate phased VCF files"; exit 1; }
echo "Combined phased VCF file created: ${FINAL_PHASED_VCF}"

# Index the final concatenated VCF
echo "Indexing final phased VCF..."
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" tabix -p vcf "${FINAL_PHASED_VCF}" \
  || { echo "Failed to index final phased VCF"; exit 1; }
echo "Final phased VCF indexed."

echo "Step 5 finished."

# Generate summary
echo "--- Integration Summary ---"
echo "Total chromosomes processed: 22"
echo "Total readcount files concatenated: ${#READCOUNT_FILES[@]}"
echo "Total VCF files concatenated: ${#PHASED_VCF_FILES[@]}"
echo ""
echo "Final output files created in: ${OUTPUT_DIR}"
echo "  - Combined readcount file: ${ALL_CHROMS_READCOUNT}"
echo "  - Combined phased VCF: ${FINAL_PHASED_VCF}"
echo "  - Combined phased VCF index: ${FINAL_PHASED_VCF}.tbi"
echo ""
echo "Integration completed successfully!"