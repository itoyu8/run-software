#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J dvr9_whphase
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
# Usage: ./dvr9_whphase.sh [--reference hg38|chm13] [--strict-filter] [-d output_dir] [-o output_name] <input.bam>
# Output: <output_dir>/<output_name>.dv.vcf.gz, <output_dir>/<output_name>.phased.vcf.gz

set -euxo pipefail

# Parse arguments
INPUT_BAM=""
OUTPUT_DIR="."
OUTPUT_NAME="output"
REFERENCE_TYPE="hg38"
STRICT_FILTER=false

while [[ $# -gt 0 ]]; do
    case $1 in
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
        --strict-filter)
            STRICT_FILTER=true
            shift
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
            if [ -z "$INPUT_BAM" ]; then
                INPUT_BAM="$1"
            else
                echo "Too many arguments"
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$INPUT_BAM" ]; then
    echo "Usage: $0 [--reference hg38|chm13] [--strict-filter] [-d output_dir] [-o output_name] <input.bam>"
    exit 1
fi

INPUT_BAM=$(realpath "$INPUT_BAM")
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")

THREADS=${SLURM_CPUS_PER_TASK:-16}

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

DVR9_WHATSHAP_SIF="/home/itoyu8/singularity/dvr9-whatshap_0.1.0.sif"
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"

DV_TEMP_DIR="${OUTPUT_DIR}/${OUTPUT_NAME}_dv_intermediate"
DV_OUTPUT="${OUTPUT_DIR}/${OUTPUT_NAME}.dv.vcf.gz"
FILTERED_OUTPUT="${OUTPUT_DIR}/${OUTPUT_NAME}.dv.filtered.vcf.gz"
PHASED_OUTPUT="${OUTPUT_DIR}/${OUTPUT_NAME}.phased.vcf.gz"

mkdir -p "${DV_TEMP_DIR}"
mkdir -p ./log

# Step 1: Run PEPPER-Margin-DeepVariant
time singularity exec \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1/ \
    "${DVR9_WHATSHAP_SIF}" run_pepper_margin_deepvariant call_variant \
    -b "${INPUT_BAM}" \
    -f "${REFERENCE_GENOME_PATH}" \
    -o "${DV_TEMP_DIR}" \
    -t "${THREADS}" \
    --ont_r9_guppy5_sup

# Move output VCF to expected location
mv "${DV_TEMP_DIR}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz" "${DV_OUTPUT}"
mv "${DV_TEMP_DIR}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz.tbi" "${DV_OUTPUT}.tbi"

# Step 2: Optional strict filter
WHATSHAP_INPUT="${DV_OUTPUT}"

if [ "$STRICT_FILTER" = true ]; then
    MIN_GQ=20
    MIN_VAF=0.3
    MAX_VAF=0.7

    time "${BCFTOOLS}" view \
        -f PASS \
        -m2 -M2 \
        --genotype het \
        "${DV_OUTPUT}" | \
    "${BCFTOOLS}" filter \
        -i "FORMAT/GQ >= ${MIN_GQ} && FORMAT/VAF >= ${MIN_VAF} && FORMAT/VAF <= ${MAX_VAF}" \
        -O z -o "${FILTERED_OUTPUT}"

    tabix -p vcf "${FILTERED_OUTPUT}"
    WHATSHAP_INPUT="${FILTERED_OUTPUT}"
fi

# Step 3: Run WhatsHap for additional phasing
time singularity exec \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1/ \
    "${DVR9_WHATSHAP_SIF}" whatshap phase \
    --reference "${REFERENCE_GENOME_PATH}" \
    --ignore-read-groups \
    --distrust-genotypes \
    -o "${PHASED_OUTPUT}" \
    "${WHATSHAP_INPUT}" \
    "${INPUT_BAM}"

# Step 4: Index phased VCF
tabix -p vcf "${PHASED_OUTPUT}"

echo "Exit status: $?"
