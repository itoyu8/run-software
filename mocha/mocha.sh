#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J mocha
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 8
# Usage: sbatch mocha.sh <GATK_VCF> <PHASED_VCF> <OUTPUT_PREFIX>

set -euxo pipefail

# Parse arguments
GATK_VCF=$1
PHASED_VCF=$2
OUTPUT_PREFIX=$3

if [ -z "$GATK_VCF" ] || [ -z "$PHASED_VCF" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Usage: sbatch $0 <GATK_VCF> <PHASED_VCF> <OUTPUT_PREFIX>"
    exit 1
fi

# Tool paths
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"
MOCHA_CONTAINER="/home/itoyu8/singularity/mocha_0.1.0.sif"

# Resource files
SEGDUPS_EXCLUDE="/home/itoyu8/database/tools/mocha/segdups_exclude.bed.gz"
CNPS="/home/itoyu8/database/tools/mocha/cnps.bed"
REFERENCE="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"

# Output directory setup
OUTPUT_DIR="${OUTPUT_PREFIX}"
mkdir -p "${OUTPUT_DIR}"

# Output files
AD_ANNOTATED_VCF="${OUTPUT_DIR}/mocha_input_ad.vcf.gz"
GC_ANNOTATED_VCF="${OUTPUT_DIR}/mocha_input.vcf.gz"
MOCHA_OUTPUT_VCF="${OUTPUT_DIR}/mocha_output.vcf.gz"
MOCHA_TSV="${OUTPUT_DIR}/mocha_calls.tsv"
MOCHA_STATS="${OUTPUT_DIR}/mocha_stats.tsv"

# Create log directory
mkdir -p log

# Step 1: Annotate phased VCF with AD from GATK VCF
STEP1_START=$(date +%s)

"${BCFTOOLS}" annotate \
    -a "${GATK_VCF}" \
    -c FORMAT/AD \
    -Oz -o "${AD_ANNOTATED_VCF}" \
    "${PHASED_VCF}"

"${BCFTOOLS}" index -t "${AD_ANNOTATED_VCF}"

STEP1_END=$(date +%s)
STEP1_ELAPSED=$((STEP1_END - STEP1_START))

# Step 2: Add GC content field
STEP2_START=$(date +%s)

singularity exec --bind /lustre1/:/lustre1/,/home/itoyu8/:/home/itoyu8/ "${MOCHA_CONTAINER}" \
    bcftools +mochatools --no-version \
        -Oz -o "${GC_ANNOTATED_VCF}" \
        "${AD_ANNOTATED_VCF}" \
        -- -t GC -f "${REFERENCE}"

"${BCFTOOLS}" index -t "${GC_ANNOTATED_VCF}"

STEP2_END=$(date +%s)
STEP2_ELAPSED=$((STEP2_END - STEP2_START))

# Step 3: Run MoChA analysis
STEP3_START=$(date +%s)

singularity exec --bind /lustre1/:/lustre1/,/home/itoyu8/:/home/itoyu8/ "${MOCHA_CONTAINER}" \
    bcftools +mocha \
        -g GRCh38 \
        -p "${CNPS}" \
        -T "^${SEGDUPS_EXCLUDE}" \
        -z "${MOCHA_STATS}" \
        -c "${MOCHA_TSV}" \
        -Oz -o "${MOCHA_OUTPUT_VCF}" \
        "${GC_ANNOTATED_VCF}"

"${BCFTOOLS}" index -t "${MOCHA_OUTPUT_VCF}"

STEP3_END=$(date +%s)
STEP3_ELAPSED=$((STEP3_END - STEP3_START))

# End time and summary
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - STEP1_START))

echo "Step 1 (AD annotation): ${STEP1_ELAPSED} seconds"
echo "Step 2 (GC annotation): ${STEP2_ELAPSED} seconds"
echo "Step 3 (MoChA analysis): ${STEP3_ELAPSED} seconds"
echo "Total time: ${TOTAL_ELAPSED} seconds"
echo "Output: ${OUTPUT_DIR}"
