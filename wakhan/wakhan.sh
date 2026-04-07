#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J wakhan
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 32
# Usage: bash wakhan/wakhan.sh [--reference hg38|chm13] [--tumor-only] [--breakpoints severus.vcf] [--cpd] [-d output_dir] [-o genome_name] --phased-vcf phased.vcf.gz tumor.bam
# Output: ${OUTPUT_DIR}/${OUTPUT_NAME}/

set -euxo pipefail

OUTPUT_DIR="."
OUTPUT_NAME="output"
REFERENCE_TYPE="hg38"
PHASED_VCF=""
BREAKPOINTS=""
CPD=""
TUMOR_ONLY=""
TUMOR_BAM=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -d)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -o)
            OUTPUT_NAME="$2"
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
        --phased-vcf)
            PHASED_VCF="$2"
            shift 2
            ;;
        --breakpoints)
            BREAKPOINTS="$2"
            shift 2
            ;;
        --cpd)
            CPD="true"
            shift
            ;;
        --tumor-only)
            TUMOR_ONLY="true"
            shift
            ;;
        *)
            TUMOR_BAM="$1"
            shift
            ;;
    esac
done

if [ -z "$PHASED_VCF" ]; then
    echo "Error: --phased-vcf is required"
    exit 1
fi

if [ -z "$TUMOR_BAM" ]; then
    echo "Error: Tumor BAM file is required"
    exit 1
fi

if [ -z "$BREAKPOINTS" ] && [ -z "$CPD" ]; then
    echo "Error: Either --breakpoints or --cpd is required"
    exit 1
fi

if [ ! -f "$PHASED_VCF" ]; then
    echo "Error: Phased VCF file not found: $PHASED_VCF"
    exit 1
fi

if [ ! -f "$TUMOR_BAM" ]; then
    echo "Error: Tumor BAM file not found: $TUMOR_BAM"
    exit 1
fi

if [ -n "$BREAKPOINTS" ] && [ ! -f "$BREAKPOINTS" ]; then
    echo "Error: Breakpoints VCF file not found: $BREAKPOINTS"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
TUMOR_BAM=$(realpath "${TUMOR_BAM}")
PHASED_VCF=$(realpath "${PHASED_VCF}")
[ -n "$BREAKPOINTS" ] && BREAKPOINTS=$(realpath "${BREAKPOINTS}")

THREADS=${SLURM_CPUS_PER_TASK:-32}

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

CONTAINER_PATH="/home/itoyu8/singularity/wakhan_0.4.2.sif"

WAKHAN_OPTS=(
    --target-bam "${TUMOR_BAM}"
    --reference "${REFERENCE_GENOME_PATH}"
    --genome-name "${OUTPUT_NAME}"
    --out-dir "${OUTPUT_DIR}"
    --threads "${THREADS}"
)

if [ -n "$TUMOR_ONLY" ]; then
    WAKHAN_OPTS+=(--tumor-phased-vcf "${PHASED_VCF}")
else
    WAKHAN_OPTS+=(--normal-phased-vcf "${PHASED_VCF}")
fi

if [ -n "$BREAKPOINTS" ]; then
    WAKHAN_OPTS+=(--breakpoints "${BREAKPOINTS}")
fi

if [ -n "$CPD" ]; then
    WAKHAN_OPTS+=(--cpd)
fi

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    WAKHAN_SITE_PACKAGES=$(singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1:/lustre1 \
        "${CONTAINER_PATH}" python -c "import wakhan; import os; print(os.path.dirname(wakhan.__file__))")
    WAKHAN_OPTS+=(
        --centromere-bed "${WAKHAN_SITE_PACKAGES}/data/chm13_centromere.bed"
        --cancer-genes "${WAKHAN_SITE_PACKAGES}/data/cancer_genes.tsv"
        --reference-name chm13
    )
fi

time singularity exec \
    --bind /home/itoyu8/:/home/itoyu8/,/lustre1:/lustre1 \
    "${CONTAINER_PATH}" \
    wakhan all \
    "${WAKHAN_OPTS[@]}"

echo "Exit status: $?"
