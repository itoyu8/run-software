#!/bin/bash
#SBATCH -p rjobs
#SBATCH -J hifiasm_trio_yak
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=7G
#SBATCH -c 56
# Usage: sbatch hifiasm_trio_yak.sh --ont <ont.fastq.gz> --h1 <h1.bam> --h2 <h2.bam> -o <output_prefix> [--l 0|1|2|3] [--dual-scaf]

set -e

# Start time
START_TIME=$(date +%s)

# Parse arguments
ONT_FASTQ=""
H1_BAM=""
H2_BAM=""
OUTPUT_PREFIX=""
PURGE_LEVEL=0  # Default purge level for trio binning
DUAL_SCAF=""   # Empty by default (off)

while [[ $# -gt 0 ]]; do
    case $1 in
        --ont)
            ONT_FASTQ="$2"
            shift 2
            ;;
        --h1)
            H1_BAM="$2"
            shift 2
            ;;
        --h2)
            H2_BAM="$2"
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
            echo "Usage: $0 --ont <ont.fastq.gz> --h1 <h1.bam> --h2 <h2.bam> -o <output_prefix> [--l 0|1|2|3] [--dual-scaf]"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$ONT_FASTQ" ] || [ -z "$H1_BAM" ] || [ -z "$H2_BAM" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --ont <ont.fastq.gz> --h1 <h1.bam> --h2 <h2.bam> -o <output_prefix> [--l 0|1|2|3] [--dual-scaf]"
    exit 1
fi

# Tool paths
SAMTOOLS="/home/itoyu8/bin/samtools/samtools-1.19/samtools"
YAK="/home/itoyu8/bin/yak/yak-0.1/yak"
HIFIASM="/home/itoyu8/bin/hifiasm/hifiasm-0.25.0/hifiasm"
THREADS=${SLURM_CPUS_PER_TASK:-56}

# Set up output directory and place all files inside it
OUTPUT_DIR="${OUTPUT_PREFIX}"
mkdir -p "${OUTPUT_DIR}"
OUTPUT_BASE=$(basename "${OUTPUT_PREFIX}")

# Define output files
H1_FASTQ="${OUTPUT_DIR}/h1.fq.gz"
H2_FASTQ="${OUTPUT_DIR}/h2.fq.gz"
H1_YAK="${OUTPUT_DIR}/h1.yak"
H2_YAK="${OUTPUT_DIR}/h2.yak"
LOG_FILE="${OUTPUT_DIR}/assembly.log"

echo "=== Converting BAM to FASTQ ==="
echo "H1 BAM: ${H1_BAM}"
echo "H2 BAM: ${H2_BAM}"
echo "Output directory: ${OUTPUT_DIR}"

# Convert H1 BAM to FASTQ
echo "Converting H1 BAM to FASTQ..."
"${SAMTOOLS}" fastq -@ "${THREADS}" "${H1_BAM}" | gzip > "${H1_FASTQ}"
H1_READS=$(zcat "${H1_FASTQ}" | wc -l | awk '{print int($1/4)}')
echo "H1 reads: ${H1_READS}"

# Convert H2 BAM to FASTQ
echo "Converting H2 BAM to FASTQ..."
"${SAMTOOLS}" fastq -@ "${THREADS}" "${H2_BAM}" | gzip > "${H2_FASTQ}"
H2_READS=$(zcat "${H2_FASTQ}" | wc -l | awk '{print int($1/4)}')
echo "H2 reads: ${H2_READS}"

echo ""
echo "=== Running yak count (k=31) ==="

# Yak count for H1
echo "Running yak count for H1: ${H1_FASTQ} -> ${H1_YAK}"
"${YAK}" count -k31 -b37 -t"${THREADS}" -o "${H1_YAK}" "${H1_FASTQ}"
if [ $? -ne 0 ]; then
    echo "Error: yak count failed for H1"
    exit 1
fi
echo "Successfully created ${H1_YAK}"

# Yak count for H2
echo "Running yak count for H2: ${H2_FASTQ} -> ${H2_YAK}"
"${YAK}" count -k31 -b37 -t"${THREADS}" -o "${H2_YAK}" "${H2_FASTQ}"
if [ $? -ne 0 ]; then
    echo "Error: yak count failed for H2"
    exit 1
fi
echo "Successfully created ${H2_YAK}"

echo ""
echo "=== Running hifiasm with yak-based trio binning ==="
echo "ONT FASTQ: ${ONT_FASTQ}"
echo "H1 yak: ${H1_YAK}"
echo "H2 yak: ${H2_YAK}"
echo "Purge level: -l${PURGE_LEVEL}"
if [ -n "$DUAL_SCAF" ]; then
    echo "Dual scaffolding: enabled"
else
    echo "Dual scaffolding: disabled"
fi
echo "Threads: ${THREADS}"

# Run hifiasm with -1/-2 options for yak-based trio binning
"${HIFIASM}" --ont -i -l${PURGE_LEVEL} ${DUAL_SCAF} -t"${THREADS}" \
    -1 "${H1_YAK}" \
    -2 "${H2_YAK}" \
    -o "${OUTPUT_DIR}/${OUTPUT_BASE}" \
    "${ONT_FASTQ}" > "${LOG_FILE}" 2>&1

if [ $? -ne 0 ]; then
    echo "Error: hifiasm failed. Check log ${LOG_FILE}"
    exit 1
fi

# End time and calculate duration
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

echo ""
echo "=== Assembly completed ==="
printf "Total time: %02d:%02d:%02d\n" $HOURS $MINUTES $SECONDS
echo "Output files are in: ${OUTPUT_DIR}"
echo "  - H1 FASTQ: ${H1_FASTQ}"
echo "  - H2 FASTQ: ${H2_FASTQ}"
echo "  - H1 yak: ${H1_YAK}"
echo "  - H2 yak: ${H2_YAK}"
echo "  - Assembly log: ${LOG_FILE}"
