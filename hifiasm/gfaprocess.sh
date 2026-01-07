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

# Tool paths
GFATOOLS="/home/itoyu8/bin/gfatools/gfatools/gfatools"
SEQKIT="/home/itoyu8/bin/seqkit/seqkit"
MINIMAP2="/home/itoyu8/bin/minimap2/minimap2-2.28/minimap2"
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
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
mkdir -p "${OUTPUT_DIR}"

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
    time "${GFATOOLS}" gfa2fa "${INPUT_FILE}" > "${FA_FILE}"
fi

# Step 2: Calculate assembly statistics
time "${SEQKIT}" stats -a "${FA_FILE}" > "${STATS_FILE}"

# Step 3: Align to reference
time "${MINIMAP2}" -x asm5 -t "${THREADS}" "${REFERENCE_PATH}" "${FA_FILE}" > "${PAF_FILE}"

echo "Exit status: $?"
