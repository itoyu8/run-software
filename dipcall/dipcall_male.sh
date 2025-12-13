#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J dipcall
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 8

# Usage: sbatch dipcall_male.sh --hap1 <hap1.fa> --hap2 <hap2.fa> -o <output_prefix> [--reference hg38|chm13]

set -e

# Start time
START_TIME=$(date +%s)

# Parse arguments
HAP1=""
HAP2=""
OUTPUT_PREFIX=""
REFERENCE_TYPE="hg38"

while [[ $# -gt 0 ]]; do
    case $1 in
        --hap1)
            HAP1="$2"
            shift 2
            ;;
        --hap2)
            HAP2="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_PREFIX="$2"
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
        *)
            echo "Error: Unknown option $1"
            echo "Usage: $0 --hap1 <hap1.fa> --hap2 <hap2.fa> -o <output_prefix> [--reference hg38|chm13]"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$HAP1" ] || [ -z "$HAP2" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --hap1 <hap1.fa> --hap2 <hap2.fa> -o <output_prefix> [--reference hg38|chm13]"
    exit 1
fi

# Set reference genome and PAR file
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
    PAR_BED="/home/itoyu8/bin/dipcall/dipcall.kit/chm13v2.0_PAR.bed"
else
    REFERENCE_GENOME="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
    PAR_BED="/home/itoyu8/bin/dipcall/dipcall.kit/hs38.PAR.bed"
fi

# Dipcall paths
DIPCALL_KIT="/home/itoyu8/bin/dipcall/dipcall.kit"
RUN_DIPCALL="${DIPCALL_KIT}/run-dipcall"

# Threads
THREADS=${SLURM_CPUS_PER_TASK:-32}

# Set up output directory and place all files inside it
OUTPUT_DIR="${OUTPUT_PREFIX}"
mkdir -p "${OUTPUT_DIR}"
OUTPUT_BASE=$(basename "${OUTPUT_PREFIX}")

# Generate Makefile
echo "Generating Makefile for dipcall..."
"${RUN_DIPCALL}" -x "${PAR_BED}" "${OUTPUT_DIR}/${OUTPUT_BASE}" "${REFERENCE_GENOME}" "${HAP1}" "${HAP2}" > "${OUTPUT_DIR}/${OUTPUT_BASE}.mak"

# Run dipcall (sequential to avoid OOM from multiple samtools sort)
echo "Running dipcall sequentially to avoid memory issues..."
make -j 1 -f "${OUTPUT_DIR}/${OUTPUT_BASE}.mak"

# End time and calculate duration
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

printf "Dipcall completed in %02d:%02d:%02d\n" $HOURS $MINUTES $SECONDS
echo "Final outputs:"
echo "  VCF: ${OUTPUT_DIR}/${OUTPUT_BASE}.dip.vcf.gz"
echo "  BED: ${OUTPUT_DIR}/${OUTPUT_BASE}.dip.bed"
