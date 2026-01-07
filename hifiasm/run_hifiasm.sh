#!/bin/bash
#SBATCH -p rjobs
#SBATCH -J hifiasm_pipeline
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=7G
#SBATCH -c 56
# Usage: sbatch run_hifiasm_pipeline.sh

set -euxo pipefail

#############################
# Edit these variables
#############################
ONT_FASTQ=""
H1_BAM=""
H2_BAM=""
OUTPUT_DIR=""
OUTPUT_NAME="output"
REFERENCE_TYPE="chm13"  # hg38 or chm13
PURGE_LEVEL=0           # 0, 1, 2, or 3
DUAL_SCAF=""            # "" or "--dual-scaf"
#############################

# Validate
if [ -z "$ONT_FASTQ" ] || [ -z "$H1_BAM" ] || [ -z "$H2_BAM" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Please set ONT_FASTQ, H1_BAM, H2_BAM, OUTPUT_DIR"
    exit 1
fi

# Get script directory
SCRIPT_DIR=$(dirname "$(realpath "$0")")

OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
mkdir -p "${OUTPUT_DIR}"

# Step 1: Run hifiasm trio assembly
bash "${SCRIPT_DIR}/hifiasm_trio_yak.sh" \
    --ont "${ONT_FASTQ}" \
    --h1 "${H1_BAM}" \
    --h2 "${H2_BAM}" \
    -d "${OUTPUT_DIR}" \
    -o "${OUTPUT_NAME}" \
    --l "${PURGE_LEVEL}" \
    ${DUAL_SCAF}

# Step 2: Process hap1 GFA
HAP1_GFA="${OUTPUT_DIR}/${OUTPUT_NAME}.dip.hap1.p_ctg.gfa"
bash "${SCRIPT_DIR}/gfaprocess.sh" \
    -d "${OUTPUT_DIR}" \
    --reference "${REFERENCE_TYPE}" \
    "${HAP1_GFA}"

# Step 3: Process hap2 GFA
HAP2_GFA="${OUTPUT_DIR}/${OUTPUT_NAME}.dip.hap2.p_ctg.gfa"
bash "${SCRIPT_DIR}/gfaprocess.sh" \
    -d "${OUTPUT_DIR}" \
    --reference "${REFERENCE_TYPE}" \
    "${HAP2_GFA}"

echo "Exit status: $?"
