#!/bin/bash
#SBATCH -p rjobs
#SBATCH -J hifiasm_ont
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=7G
#SBATCH -c 56
# Usage: sbatch hifiasm_ontonly.sh --ont <ont.fastq.gz> -d <output_dir> -o <basename> [--l 0|1|2|3] [--dual-scaf]
# Output: <basename>.bp.p_ctg.gfa, <basename>.bp.p_utg.gfa

set -euxo pipefail

# Parse arguments
ONT_FASTQ=""
OUTPUT_DIR=""
OUTPUT_NAME="output"
PURGE_LEVEL=3
DUAL_SCAF=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --ont)
            ONT_FASTQ="$2"
            shift 2
            ;;
        -d)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -o)
            OUTPUT_NAME="$2"
            shift 2
            ;;
        --l)
            if [[ "$2" =~ ^[0-3]$ ]]; then
                PURGE_LEVEL="$2"
            else
                echo "Error: --l must be 0, 1, 2, or 3"
                exit 1
            fi
            shift 2
            ;;
        --dual-scaf)
            DUAL_SCAF="--dual-scaf"
            shift
            ;;
        *)
            echo "Error: Unknown option $1"
            echo "Usage: $0 --ont <ont.fastq.gz> -d <output_dir> -o <basename> [--l 0|1|2|3] [--dual-scaf]"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$ONT_FASTQ" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 --ont <ont.fastq.gz> -d <output_dir> -o <basename> [--l 0|1|2|3] [--dual-scaf]"
    exit 1
fi

# Tool paths
HIFIASM_PATH="/home/itoyu8/bin/hifiasm/hifiasm-0.25.0/hifiasm"
THREADS=${SLURM_CPUS_PER_TASK:-56}

# Set up output directory
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
LOG_FILE="${OUTPUT_DIR}/assembly.log"

# Run hifiasm
time "${HIFIASM_PATH}" --ont -i -l${PURGE_LEVEL} ${DUAL_SCAF} -t"${THREADS}" \
    -o "${OUTPUT_DIR}/${OUTPUT_NAME}" \
    "${ONT_FASTQ}" > "${LOG_FILE}" 2>&1

echo "Exit status: $?"
