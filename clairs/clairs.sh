#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J clairs
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 32
# Usage1: bash clairs.sh -d /output/dir tumor.bam normal.bam
# Usage2: bash clairs.sh -d /output/dir --normal-vcf germline.vcf.gz tumor.bam normal.bam
# Usage3: bash clairs.sh -d /output/dir --normal-vcf germline.vcf.gz --haplotagged tumor_haplotagged.bam normal.bam
# Output: ${OUTPUT_DIR}/output.vcf.gz

set -euxo pipefail

OUTPUT_DIR="."
REFERENCE_TYPE="chm13"
PLATFORM="ont"
HAPLOTAGGED=""
NORMAL_VCF=""
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
        --haplotagged)
            HAPLOTAGGED="true"
            shift
            ;;
        --normal-vcf)
            NORMAL_VCF="$2"
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

if [ -n "$HAPLOTAGGED" ] && [ -z "$NORMAL_VCF" ]; then
    echo "Error: --haplotagged requires --normal-vcf"
    exit 1
fi

# Note: --normal-vcf alone is valid (skip Clair3 germline calling, but do phasing)

if [ -n "$NORMAL_VCF" ] && [ ! -f "$NORMAL_VCF" ]; then
    echo "Error: Normal VCF file not found: $NORMAL_VCF"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
TUMOR_BAM=$(realpath "${TUMOR_BAM}")
NORMAL_BAM=$(realpath "${NORMAL_BAM}")
[ -n "$NORMAL_VCF" ] && NORMAL_VCF=$(realpath "${NORMAL_VCF}")

THREADS=${SLURM_CPUS_PER_TASK:-32}

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

if [ "$PLATFORM" = "ont" ]; then
    PLATFORM_MODEL="ont_r10_dorado_sup_5khz"
else
    PLATFORM_MODEL="hifi_revio"
fi

CONTAINER_PATH="/home/itoyu8/singularity/clairs_0.1.0.sif"

BIND_PATHS="${OUTPUT_DIR}:${OUTPUT_DIR}"
BIND_PATHS="${BIND_PATHS},$(dirname "${TUMOR_BAM}"):$(dirname "${TUMOR_BAM}")"
BIND_PATHS="${BIND_PATHS},$(dirname "${NORMAL_BAM}"):$(dirname "${NORMAL_BAM}")"
BIND_PATHS="${BIND_PATHS},$(dirname "${REFERENCE_GENOME_PATH}"):$(dirname "${REFERENCE_GENOME_PATH}")"
[ -n "$NORMAL_VCF" ] && BIND_PATHS="${BIND_PATHS},$(dirname "${NORMAL_VCF}"):$(dirname "${NORMAL_VCF}")"

CLAIRS_OPTS=(
    --tumor_bam_fn="${TUMOR_BAM}"
    --normal_bam_fn="${NORMAL_BAM}"
    --ref_fn="${REFERENCE_GENOME_PATH}"
    --threads="${THREADS}"
    --platform="${PLATFORM_MODEL}"
    --output_dir="${OUTPUT_DIR}"
    --enable_indel_calling
)

if [ -n "$NORMAL_VCF" ]; then
    CLAIRS_OPTS+=(--normal_vcf_fn="${NORMAL_VCF}")
fi

if [ -n "$HAPLOTAGGED" ]; then
    CLAIRS_OPTS+=(--haplotagged_tumor_bam_provided_so_skip_intermediate_phasing_and_haplotagging)
fi

time singularity exec \
    --bind "${BIND_PATHS}" \
    "${CONTAINER_PATH}" \
    /opt/bin/run_clairs \
    "${CLAIRS_OPTS[@]}"

echo "Exit status: $?"
