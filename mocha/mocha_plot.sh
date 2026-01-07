#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J mocha_plot
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 4
# Usage: sbatch mocha_plot.sh --chr <1-22> <sample_name>
# Output: /home/itoyu8/project/genome_methy/sampleselect/output_mocha/<sample_name>/figure/chr<N>.png

set -euxo pipefail

# GRCh38 chromosome lengths (chr1-22)
# Index 0 is dummy for 1-indexed access
CHR_LEN=(
    0
    248956422 242193529 198295559 190214555 181538259
    170805979 159345973 145138636 138394717 133797422
    135086622 133275309 114364328 107043718 101991189
    90338345 83257441 80373285 58617616 64444167
    46709983 50818468
)

# Tool paths
MOCHA_CONTAINER="/home/itoyu8/singularity/mocha_0.1.0.sif"
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"

# Base directory
BASE_DIR="/home/itoyu8/project/genome_methy/sampleselect/output_mocha"

# Default values
CHR=""
SAMPLE_NAME=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --chr)
            CHR="$2"
            shift 2
            ;;
        *)
            SAMPLE_NAME="$1"
            shift
            ;;
    esac
done

# Validation
if [ -z "$CHR" ]; then
    echo "Error: --chr is required"
    exit 1
fi

if [ -z "$SAMPLE_NAME" ]; then
    echo "Error: sample_name is required"
    exit 1
fi

if [ "$CHR" -lt 1 ] || [ "$CHR" -gt 22 ]; then
    echo "Error: --chr must be between 1 and 22"
    exit 1
fi

# Setup paths
SAMPLE_DIR="${BASE_DIR}/${SAMPLE_NAME}"
INPUT_VCF="${SAMPLE_DIR}/mocha_output.vcf.gz"
OUTPUT_DIR="${SAMPLE_DIR}/figure"
mkdir -p "${OUTPUT_DIR}"

# Validate input file exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF not found: ${INPUT_VCF}"
    exit 1
fi

# Get sample name from VCF header (first sample column)
VCF_SAMPLE=$("${BCFTOOLS}" query -l "${INPUT_VCF}" | head -n1)

# Build region string
REGION="chr${CHR}:1-${CHR_LEN[$CHR]}"
OUTPUT_FILE="${OUTPUT_DIR}/chr${CHR}.png"

# Create filtered VCF (exclude variants with missing AD)
FILTERED_VCF="${SAMPLE_DIR}/mocha_output_filtered_chr${CHR}.vcf.gz"
"${BCFTOOLS}" view -i 'FORMAT/AD[0:0]!="." & FORMAT/AD[0:1]!="."' \
    -r "${REGION}" \
    -Oz -o "${FILTERED_VCF}" \
    "${INPUT_VCF}"
"${BCFTOOLS}" index -t "${FILTERED_VCF}"

# Run mocha_plot.R
time singularity exec --bind /lustre1/:/lustre1/,/home/itoyu8/:/home/itoyu8/ "${MOCHA_CONTAINER}" \
    /usr/local/share/mocha/mocha_plot.R \
        --genome GRCh38 \
        --wgs \
        --mocha \
        --vcf "${FILTERED_VCF}" \
        --samples "${VCF_SAMPLE}" \
        --regions "${REGION}" \
        --png "${OUTPUT_FILE}"

# Cleanup temporary filtered VCF
rm -f "${FILTERED_VCF}" "${FILTERED_VCF}.tbi"

echo "Exit status: $?"
