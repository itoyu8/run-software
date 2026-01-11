#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J gfaprocess
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=16G
#SBATCH -c 4
# Usage: sbatch gfaprocess.sh -d <output_dir> --reference <hg38|chm13> <input.gfa|input.fa|input.fasta>
# Output: <basename>.fa (if gfa input), <basename>.stats.tsv, <basename>.paf

set -euxo pipefail

# Container
LH3_TOOLS_SIF="/home/itoyu8/singularity/lh3-tools_0.1.0.sif"
THREADS=${SLURM_CPUS_PER_TASK:-16}

# Parse arguments
OUTPUT_DIR=""
REFERENCE_TYPE=""
INPUT_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -d)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --reference)
            if [ "$2" = "chm13" ] || [ "$2" = "hg38" ]; then
                REFERENCE_TYPE="$2"
            else
                echo "Error: --reference must be 'hg38' or 'chm13'"
                exit 1
            fi
            shift 2
            ;;
        *)
            INPUT_FILE="$1"
            shift
            ;;
    esac
done

# Validate required arguments
if [ -z "$OUTPUT_DIR" ] || [ -z "$REFERENCE_TYPE" ] || [ -z "$INPUT_FILE" ]; then
    echo "Usage: $0 -d <output_dir> --reference <hg38|chm13> <input.gfa|input.fa|input.fasta>"
    exit 1
fi

# Set reference path
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

# Set up output directory
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")

# Determine input type and set basename/FA_FILE accordingly
case "${INPUT_FILE}" in
    *.gfa)
        BASENAME=$(basename "${INPUT_FILE}" .gfa)
        FA_FILE="${OUTPUT_DIR}/${BASENAME}.fa"
        IS_GFA=true
        ;;
    *.fasta)
        BASENAME=$(basename "${INPUT_FILE}" .fasta)
        FA_FILE=$(realpath "${INPUT_FILE}")
        IS_GFA=false
        ;;
    *.fa)
        BASENAME=$(basename "${INPUT_FILE}" .fa)
        FA_FILE=$(realpath "${INPUT_FILE}")
        IS_GFA=false
        ;;
    *)
        echo "Error: Input file must be .gfa, .fa, or .fasta"
        exit 1
        ;;
esac

# Output file paths
STATS_FILE="${OUTPUT_DIR}/${BASENAME}.stats.tsv"
PAF_FILE="${OUTPUT_DIR}/${BASENAME}.paf"

# Step 1: Convert GFA to FASTA (skip if input is already FASTA)
if [ "$IS_GFA" = true ]; then
    time singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ \
        "${LH3_TOOLS_SIF}" gfatools gfa2fa "${INPUT_FILE}" > "${FA_FILE}"
fi

# Step 2: Calculate assembly statistics
time singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ \
    "${LH3_TOOLS_SIF}" seqkit stats -a "${FA_FILE}" > "${STATS_FILE}"

# Step 3: Align to reference
time singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ \
    "${LH3_TOOLS_SIF}" minimap2 -x asm5 -t "${THREADS}" "${REFERENCE_PATH}" "${FA_FILE}" > "${PAF_FILE}"

echo "Exit status: $?"
