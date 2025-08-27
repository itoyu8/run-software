#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J minimap2_samsort
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 32

# Usage: ./minimap2_samsort.sh --type <ont|hifi> [--reference hg38|chm13] <input_fastq> <output_base_name>

# Parse arguments
SEQ_TYPE=""
INPUT_FASTQ=""
OUTPUT_BASE_NAME=""
REFERENCE_TYPE="hg38"

while [[ $# -gt 0 ]]; do
    case $1 in
        --type)
            if [ "$2" = "ont" ] || [ "$2" = "hifi" ]; then
                SEQ_TYPE="$2"
            else
                echo "Error: --type must be either 'ont' or 'hifi'"
                echo "Usage: $0 --type <ont|hifi> [--reference hg38|chm13] <input_fastq> <output_base_name>"
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
                echo "Usage: $0 --type <ont|hifi> [--reference hg38|chm13] <input_fastq> <output_base_name>"
                exit 1
            fi
            shift 2
            ;;
        -*)
            echo "Unknown option $1"
            echo "Usage: $0 --type <ont|hifi> [--reference hg38|chm13] <input_fastq> <output_base_name>"
            exit 1
            ;;
        *)
            if [ -z "$INPUT_FASTQ" ]; then
                INPUT_FASTQ="$1"
            elif [ -z "$OUTPUT_BASE_NAME" ]; then
                OUTPUT_BASE_NAME="$1"
            else
                echo "Too many arguments"
                echo "Usage: $0 --type <ont|hifi> [--reference hg38|chm13] <input_fastq> <output_base_name>"
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$SEQ_TYPE" ] || [ -z "$INPUT_FASTQ" ] || [ -z "$OUTPUT_BASE_NAME" ]; then
    echo "Usage: $0 --type <ont|hifi> [--reference hg38|chm13] <input_fastq> <output_base_name>"
    exit 1
fi
THREADS=${SLURM_CPUS_PER_TASK:-32}

INPUT_DIR=$(dirname "$INPUT_FASTQ")

if [ "$SEQ_TYPE" = "ont" ]; then
    MINIMAP2_PRESET="-ax map-ont"
elif [ "$SEQ_TYPE" = "hifi" ]; then
    MINIMAP2_PRESET="-ax map-pb"
fi

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

REFERENCE_DIR=$(dirname "${REFERENCE_GENOME_PATH}")
REFERENCE_BASENAME=$(basename "${REFERENCE_GENOME_PATH}" .fasta)
REFERENCE_MMI_PATH="${REFERENCE_DIR}/${REFERENCE_BASENAME}.mmi"

TEMP_BAM_FILE="${INPUT_DIR}/${OUTPUT_BASE_NAME}.temp.bam"
FINAL_BAM_FILE="${INPUT_DIR}/${OUTPUT_BASE_NAME}.bam"
FINAL_BAM_INDEX="${INPUT_DIR}/${OUTPUT_BASE_NAME}.bam.bai"

mkdir -p ./log

if [ ! -f "${REFERENCE_MMI_PATH}" ]; then
    /home/itoyu8/bin/minimap2/minimap2-2.28/minimap2 -t "${THREADS}" -d "${REFERENCE_MMI_PATH}" "${REFERENCE_GENOME_PATH}" \
        || { echo "ERROR: Reference genome indexing failed."; exit 1; }
fi

/home/itoyu8/bin/minimap2/minimap2-2.28/minimap2 ${MINIMAP2_PRESET} -t ${THREADS} "${REFERENCE_MMI_PATH}" "${INPUT_FASTQ}" | \
/home/itoyu8/bin/samtools/samtools-1.19/samtools view -bh -@ ${THREADS} - \
    > "${TEMP_BAM_FILE}" \
    || { echo "ERROR: Alignment and BAM conversion failed."; exit 1; }

/home/itoyu8/bin/samtools/samtools-1.19/samtools sort -@ "${THREADS}" -o "${FINAL_BAM_FILE}" "${TEMP_BAM_FILE}" \
    || { echo "ERROR: BAM sorting failed."; exit 1; }

/home/itoyu8/bin/samtools/samtools-1.19/samtools index "${FINAL_BAM_FILE}" "${FINAL_BAM_INDEX}" \
    || { echo "ERROR: BAM indexing failed."; exit 1; }

rm -f "${TEMP_BAM_FILE}"