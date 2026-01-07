#!/bin/bash
#SBATCH -p rjobs
#SBATCH -J hifiasm_trio
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=7G
#SBATCH -c 56
# Usage: sbatch hifiasm_trio_bamlist.sh --ont <ont.fastq.gz> --bam <haplotagged.bam> -d <output_dir> -o <basename> [--l 0|1|2|3] [--dual-scaf]
# Output: <basename>.dip.hap1.p_ctg.gfa, <basename>.dip.hap2.p_ctg.gfa, <basename>.dip.p_utg.gfa

set -euxo pipefail

# Parse arguments
ONT_FASTQ=""
PHASED_BAM=""
OUTPUT_DIR=""
OUTPUT_NAME="output"
PURGE_LEVEL=0
DUAL_SCAF=""

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
            echo "Usage: $0 --ont <ont.fastq.gz> --bam <haplotagged.bam> -d <output_dir> -o <basename> [--l 0|1|2|3] [--dual-scaf]"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$ONT_FASTQ" ] || [ -z "$PHASED_BAM" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 --ont <ont.fastq.gz> --bam <haplotagged.bam> -d <output_dir> -o <basename> [--l 0|1|2|3] [--dual-scaf]"
    exit 1
fi

# Tool paths
CONTAINER="/home/itoyu8/singularity/compat_parabricks-0.2.2.sif"
HIFIASM="/home/itoyu8/bin/hifiasm/hifiasm-0.25.0/hifiasm"
THREADS=${SLURM_CPUS_PER_TASK:-56}

# Set up output directory
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
mkdir -p "${OUTPUT_DIR}"

# Define output files
HAP1_READS="${OUTPUT_DIR}/hap1_reads.txt"
HAP2_READS="${OUTPUT_DIR}/hap2_reads.txt"
LOG_FILE="${OUTPUT_DIR}/assembly.log"

# Extract HP:1 (haplotype 1) read names
singularity exec "${CONTAINER}" samtools view -d HP:1 "${PHASED_BAM}" | cut -f1 > "${HAP1_READS}"

# Extract HP:2 (haplotype 2) read names
singularity exec "${CONTAINER}" samtools view -d HP:2 "${PHASED_BAM}" | cut -f1 > "${HAP2_READS}"

# Run hifiasm with -3/-4 options for pseudo trio binning
time "${HIFIASM}" --ont -i -l${PURGE_LEVEL} ${DUAL_SCAF} -t"${THREADS}" \
    -3 "${HAP1_READS}" \
    -4 "${HAP2_READS}" \
    -o "${OUTPUT_DIR}/${OUTPUT_NAME}" \
    "${ONT_FASTQ}" > "${LOG_FILE}" 2>&1

echo "Exit status: $?"
