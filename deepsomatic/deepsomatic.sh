#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J deepsomatic
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 32
# Usage: bash deepsomatic.sh -d /output/dir [--reference hg38|chm13] [--platform ont|hifi] tumor.bam normal.bam
# Output: ${OUTPUT_DIR}/output.vcf.gz

set -euxo pipefail

OUTPUT_DIR="."
REFERENCE_TYPE="chm13"
PLATFORM="ont"
INPUT_FILES=()

while [[ $# -gt 0 ]]; do
    case $1 in
        -d)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --reference)
            if [ "$2" = "chm13" ]; then
                REFERENCE_TYPE="chm13"
            elif [ "$2" = "hg38" ]; then
                REFERENCE_TYPE="hg38"
            else
                echo "Error: --reference must be 'hg38' or 'chm13'"
                exit 1
            fi
            shift 2
            ;;
        --platform)
            if [ "$2" = "ont" ] || [ "$2" = "hifi" ]; then
                PLATFORM="$2"
            else
                echo "Error: --platform must be 'ont' or 'hifi'"
                exit 1
            fi
            shift 2
            ;;
        *)
            INPUT_FILES+=("$1")
            shift
            ;;
    esac
done

if [ ${#INPUT_FILES[@]} -ne 2 ]; then
    echo "Error: Requires exactly 2 BAM files (tumor.bam normal.bam)"
    exit 1
fi

TUMOR_BAM="${INPUT_FILES[0]}"
NORMAL_BAM="${INPUT_FILES[1]}"

if [ ! -f "$TUMOR_BAM" ]; then
    echo "Error: Tumor BAM file not found: $TUMOR_BAM"
    exit 1
fi

if [ ! -f "$NORMAL_BAM" ]; then
    echo "Error: Normal BAM file not found: $NORMAL_BAM"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
TUMOR_BAM=$(realpath "${TUMOR_BAM}")
NORMAL_BAM=$(realpath "${NORMAL_BAM}")

THREADS=${SLURM_CPUS_PER_TASK:-32}

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

if [ "$PLATFORM" = "ont" ]; then
    MODEL_TYPE="ONT"
else
    MODEL_TYPE="PACBIO"
fi

CONTAINER_PATH="/home/itoyu8/singularity/deepsomatic_0.1.0.sif"

BIND_PATHS="${OUTPUT_DIR}:${OUTPUT_DIR}"
BIND_PATHS="${BIND_PATHS},$(dirname "${TUMOR_BAM}"):$(dirname "${TUMOR_BAM}")"
BIND_PATHS="${BIND_PATHS},$(dirname "${NORMAL_BAM}"):$(dirname "${NORMAL_BAM}")"
BIND_PATHS="${BIND_PATHS},$(dirname "${REFERENCE_GENOME_PATH}"):$(dirname "${REFERENCE_GENOME_PATH}")"

time singularity exec \
    --bind "${BIND_PATHS}" \
    "${CONTAINER_PATH}" \
    run_deepsomatic \
    --model_type="${MODEL_TYPE}" \
    --ref="${REFERENCE_GENOME_PATH}" \
    --reads_tumor="${TUMOR_BAM}" \
    --reads_normal="${NORMAL_BAM}" \
    --sample_name_tumor="tumor" \
    --sample_name_normal="normal" \
    --output_vcf="${OUTPUT_DIR}/output.vcf.gz" \
    --output_gvcf="${OUTPUT_DIR}/output.g.vcf.gz" \
    --num_shards="${THREADS}"

echo "Exit status: $?"
