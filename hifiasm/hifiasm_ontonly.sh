#!/bin/bash
#SBATCH -p rjobs
#SBATCH -J hifiasm_ont
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=7G
#SBATCH -c 56

# Usage: sbatch hifiasm_ontonly.sh --ont <ont.fastq.gz> -o <output_prefix> [--l 0|1|2|3] [--dual-scaf]
# !!CAUTION!! use -o ~/absolute_path/output_prefix

set -e

# Start time
START_TIME=$(date +%s)

# Parse arguments
ONT_FASTQ=""
OUTPUT_PREFIX=""
PURGE_LEVEL=3  # Default purge level
DUAL_SCAF=""   # Empty by default (off)

while [[ $# -gt 0 ]]; do
    case $1 in
        --ont)
            ONT_FASTQ="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_PREFIX="$2"
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
            echo "Usage: $0 --ont <ont.fastq.gz> -o <output_prefix> [--l 0|1|2|3] [--dual-scaf]"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$ONT_FASTQ" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --ont <ont.fastq.gz> -o <output_prefix> [--l 0|1|2|3] [--dual-scaf]"
    exit 1
fi

# Hifiasm path and settings
HIFIASM_PATH="/home/itoyu8/bin/hifiasm/hifiasm-0.25.0/hifiasm"
THREADS=${SLURM_CPUS_PER_TASK:-56}

# Set up output directory and log file
# Create output directory and place all files inside it
OUTPUT_DIR="${OUTPUT_PREFIX}"
mkdir -p "${OUTPUT_DIR}"
OUTPUT_BASE=$(basename "${OUTPUT_PREFIX}")
LOG_FILE="${OUTPUT_DIR}/assembly.log"

# Log parameters
echo "Purge level: -l${PURGE_LEVEL}"
if [ -n "$DUAL_SCAF" ]; then
    echo "Dual scaffolding: enabled"
else
    echo "Dual scaffolding: disabled"
fi

# Run hifiasm (output files will be in OUTPUT_DIR/)
"${HIFIASM_PATH}" --ont -i -l${PURGE_LEVEL} ${DUAL_SCAF} -t"${THREADS}" \
    -o "${OUTPUT_DIR}/${OUTPUT_BASE}" \
    "${ONT_FASTQ}" > "${LOG_FILE}" 2>&1

# End time and calculate duration
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

printf "Hifiasm assembly completed in %02d:%02d:%02d\n" $HOURS $MINUTES $SECONDS
