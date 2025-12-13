#!/bin/bash
#SBATCH -p rjobs
#SBATCH -J hifiasm_trio
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=7G
#SBATCH -c 56
# Usage: sbatch hifiasm_trio_bamlist.sh --ont <ont.fastq.gz> --bam <haplotagged.bam> -o <output_prefix> [--l 0|1|2|3] [--dual-scaf]

set -e

# Start time
START_TIME=$(date +%s)

# Parse arguments
ONT_FASTQ=""
PHASED_BAM=""
OUTPUT_PREFIX=""
PURGE_LEVEL=0  # Default purge level for trio binning
DUAL_SCAF=""   # Empty by default (off)

while [[ $# -gt 0 ]]; do
    case $1 in
        --ont)
            ONT_FASTQ="$2"
            shift 2
            ;;
        --bam)
            PHASED_BAM="$2"
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
            echo "Usage: $0 --ont <ont.fastq.gz> --bam <haplotagged.bam> -o <output_prefix> [--l 0|1|2|3] [--dual-scaf]"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$ONT_FASTQ" ] || [ -z "$PHASED_BAM" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --ont <ont.fastq.gz> --bam <haplotagged.bam> -o <output_prefix> [--l 0|1|2|3] [--dual-scaf]"
    exit 1
fi

# Tool paths
SAMTOOLS="/home/itoyu8/bin/samtools/samtools-1.19/samtools"
HIFIASM="/home/itoyu8/bin/hifiasm/hifiasm-0.25.0/hifiasm"
THREADS=${SLURM_CPUS_PER_TASK:-56}

# Set up output directory and place all files inside it
OUTPUT_DIR="${OUTPUT_PREFIX}"
mkdir -p "${OUTPUT_DIR}"
OUTPUT_BASE=$(basename "${OUTPUT_PREFIX}")

# Define output files
HAP1_READS="${OUTPUT_DIR}/hap1_reads.txt"
HAP2_READS="${OUTPUT_DIR}/hap2_reads.txt"
LOG_FILE="${OUTPUT_DIR}/assembly.log"

echo "=== Extracting haplotype-tagged reads from BAM ==="
echo "Haplotagged BAM: ${PHASED_BAM}"
echo "Output directory: ${OUTPUT_DIR}"

# Extract HP:1 (haplotype 1) read names
echo "Extracting HP:1 read names..."
"${SAMTOOLS}" view -d HP:1 "${PHASED_BAM}" | cut -f1 > "${HAP1_READS}"
HAP1_COUNT=$(wc -l < "${HAP1_READS}")
echo "HP:1 reads: ${HAP1_COUNT}"

# Extract HP:2 (haplotype 2) read names
echo "Extracting HP:2 read names..."
"${SAMTOOLS}" view -d HP:2 "${PHASED_BAM}" | cut -f1 > "${HAP2_READS}"
HAP2_COUNT=$(wc -l < "${HAP2_READS}")
echo "HP:2 reads: ${HAP2_COUNT}"

echo ""
echo "=== Running hifiasm with pseudo trio binning ==="
echo "ONT FASTQ: ${ONT_FASTQ}"
echo "Purge level: -l${PURGE_LEVEL}"
if [ -n "$DUAL_SCAF" ]; then
    echo "Dual scaffolding: enabled"
else
    echo "Dual scaffolding: disabled"
fi
echo "Threads: ${THREADS}"

# Run hifiasm with -3/-4 options for pseudo trio binning
"${HIFIASM}" --ont -i -l${PURGE_LEVEL} ${DUAL_SCAF} -t"${THREADS}" \
    -3 "${HAP1_READS}" \
    -4 "${HAP2_READS}" \
    -o "${OUTPUT_DIR}/${OUTPUT_BASE}" \
    "${ONT_FASTQ}" > "${LOG_FILE}" 2>&1

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
