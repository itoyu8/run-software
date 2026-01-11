#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J minimap2_samsort
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 16
# Usage: ./minimap2_samsort.sh --type <ont|hifi> [--reference hg38|chm13] [-d output_dir] [-o output_name] <input_fastq>
# Output: <output_dir>/<output_name>.bam, <output_dir>/<output_name>.bam.bai

set -euxo pipefail

# Parse arguments
SEQ_TYPE=""
INPUT_FASTQ=""
OUTPUT_DIR="."
OUTPUT_NAME="output"
REFERENCE_TYPE="hg38"

while [[ $# -gt 0 ]]; do
    case $1 in
        --type)
            if [ "$2" = "ont" ] || [ "$2" = "hifi" ]; then
                SEQ_TYPE="$2"
            else
                echo "Error: --type must be either 'ont' or 'hifi'"
                exit 1
            fi
            shift 2
            ;;
        --reference)
            if [ "$2" = "chm13" ]; then
                REFERENCE_TYPE="chm13"
            elif [ "$2" = "hg38" ]; then
                REFERENCE_TYPE="hg38"
            else
                echo "Error: --reference must be either 'hg38' or 'chm13'"
                exit 1
            fi
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
        -*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            if [ -z "$INPUT_FASTQ" ]; then
                INPUT_FASTQ="$1"
            else
                echo "Too many arguments"
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$SEQ_TYPE" ] || [ -z "$INPUT_FASTQ" ]; then
    echo "Usage: $0 --type <ont|hifi> [--reference hg38|chm13] [-d output_dir] [-o output_name] <input_fastq>"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")

THREADS=${SLURM_CPUS_PER_TASK:-16}

if [ "$SEQ_TYPE" = "ont" ]; then
    MINIMAP2_PRESET="-ax map-ont"
elif [ "$SEQ_TYPE" = "hifi" ]; then
    MINIMAP2_PRESET="-ax map-hifi"
fi

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

REFERENCE_DIR=$(dirname "${REFERENCE_GENOME_PATH}")
REFERENCE_BASENAME=$(basename "${REFERENCE_GENOME_PATH}" .fa)
REFERENCE_MMI_PATH="${REFERENCE_DIR}/${REFERENCE_BASENAME}.mmi"

TEMP_BAM_FILE="${OUTPUT_DIR}/${OUTPUT_NAME}.temp.bam"
FINAL_BAM_FILE="${OUTPUT_DIR}/${OUTPUT_NAME}.bam"
FINAL_BAM_INDEX="${OUTPUT_DIR}/${OUTPUT_NAME}.bam.bai"

mkdir -p ./log

if [ ! -f "${REFERENCE_MMI_PATH}" ]; then
    /home/itoyu8/bin/minimap2/minimap2-2.28/minimap2 -t "${THREADS}" -d "${REFERENCE_MMI_PATH}" "${REFERENCE_GENOME_PATH}"
fi

time /home/itoyu8/bin/minimap2/minimap2-2.28/minimap2 ${MINIMAP2_PRESET} -t ${THREADS} "${REFERENCE_MMI_PATH}" "${INPUT_FASTQ}" | \
    /home/itoyu8/bin/samtools/samtools-1.19/samtools view -bh -@ ${THREADS} - > "${TEMP_BAM_FILE}"

time /home/itoyu8/bin/samtools/samtools-1.19/samtools sort -@ "${THREADS}" -o "${FINAL_BAM_FILE}" "${TEMP_BAM_FILE}"

/home/itoyu8/bin/samtools/samtools-1.19/samtools index "${FINAL_BAM_FILE}" "${FINAL_BAM_INDEX}"

rm -f "${TEMP_BAM_FILE}"

echo "Exit status: $?"
