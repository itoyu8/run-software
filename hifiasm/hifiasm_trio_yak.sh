#!/bin/bash
#SBATCH -p rjobs
#SBATCH -J hifiasm_trio_yak
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=7G
#SBATCH -c 56
# Usage: sbatch hifiasm_trio_yak.sh --ont <ont.fastq.gz> --h1 <h1.bam> --h2 <h2.bam> -d <output_dir> -o <basename> [--l 0|1|2|3] [--dual-scaf]
# Output: <basename>.dip.hap1.p_ctg.gfa, <basename>.dip.hap2.p_ctg.gfa, <basename>.dip.p_utg.gfa

set -euxo pipefail

# Parse arguments
ONT_FASTQ=""
H1_BAM=""
H2_BAM=""
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
        --h1)
            H1_BAM="$2"
            shift 2
            ;;
        --h2)
            H2_BAM="$2"
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
            echo "Usage: $0 --ont <ont.fastq.gz> --h1 <h1.bam> --h2 <h2.bam> -d <output_dir> -o <basename> [--l 0|1|2|3] [--dual-scaf]"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$ONT_FASTQ" ] || [ -z "$H1_BAM" ] || [ -z "$H2_BAM" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 --ont <ont.fastq.gz> --h1 <h1.bam> --h2 <h2.bam> -d <output_dir> -o <basename> [--l 0|1|2|3] [--dual-scaf]"
    exit 1
fi

# Tool paths
CONTAINER="/home/itoyu8/singularity/compat_parabricks-0.2.2.sif"
YAK="/home/itoyu8/bin/yak/yak-0.1/yak"
HIFIASM="/home/itoyu8/bin/hifiasm/hifiasm-0.25.0/hifiasm"
THREADS=${SLURM_CPUS_PER_TASK:-56}

# Set up output directory
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
mkdir -p "${OUTPUT_DIR}"

# Define output files
H1_FASTQ="${OUTPUT_DIR}/h1.fq.gz"
H2_FASTQ="${OUTPUT_DIR}/h2.fq.gz"
H1_YAK="${OUTPUT_DIR}/h1.yak"
H2_YAK="${OUTPUT_DIR}/h2.yak"
LOG_FILE="${OUTPUT_DIR}/assembly.log"

# Convert H1 BAM to FASTQ
time singularity exec "${CONTAINER}" samtools fastq -@ "${THREADS}" "${H1_BAM}" | gzip > "${H1_FASTQ}"

# Convert H2 BAM to FASTQ
time singularity exec "${CONTAINER}" samtools fastq -@ "${THREADS}" "${H2_BAM}" | gzip > "${H2_FASTQ}"

# Yak count for H1
time "${YAK}" count -k31 -b37 -t"${THREADS}" -o "${H1_YAK}" "${H1_FASTQ}"

# Yak count for H2
time "${YAK}" count -k31 -b37 -t"${THREADS}" -o "${H2_YAK}" "${H2_FASTQ}"

# Run hifiasm with -1/-2 options for yak-based trio binning
time "${HIFIASM}" --ont -i -l${PURGE_LEVEL} ${DUAL_SCAF} -t"${THREADS}" \
    -1 "${H1_YAK}" \
    -2 "${H2_YAK}" \
    -o "${OUTPUT_DIR}/${OUTPUT_NAME}" \
    "${ONT_FASTQ}" > "${LOG_FILE}" 2>&1

echo "Exit status: $?"
