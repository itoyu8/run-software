#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J whtag_split
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=64G
#SBATCH -c 16
# Usage: sbatch whtag_split.sh [--reference hg38|chm13] <PHASED_VCF> <BAM/CRAM_FILE>

# whatshap does not allow multi-threading system (except polyphase)
# CRAM files are converted to BAM first using samtools for faster processing

# Start time
START_TIME=$(date +%s)

# Parse arguments
REFERENCE_TYPE="hg38"

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
        *)
            break
            ;;
    esac
done

PHASED_VCF=$1
INPUT_FILE=$2

if [ -z "$PHASED_VCF" ] || [ -z "$INPUT_FILE" ]; then
    echo "Usage: sbatch $0 [--reference hg38|chm13] <PHASED_VCF> <BAM/CRAM_FILE>"
    exit 1
fi

# Detect file format and extract basename without extension
if [[ "$INPUT_FILE" == *.cram ]]; then
    FILE_BASENAME=$(basename "$INPUT_FILE" .cram)
    FILE_EXT="cram"
elif [[ "$INPUT_FILE" == *.bam ]]; then
    FILE_BASENAME=$(basename "$INPUT_FILE" .bam)
    FILE_EXT="bam"
else
    echo "Error: Input file must be .bam or .cram"
    exit 1
fi

# Set up output directory in the same directory as input file
INPUT_DIR=$(dirname "$INPUT_FILE")
OUTPUT_DIR="${INPUT_DIR}/${FILE_BASENAME}_hptag"
FILE_PREFIX="${FILE_BASENAME}"

# Set reference-specific paths
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REF_FASTA_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
    GENOME_FILE="/home/itoyu8/database/reference/chm13/v2.0/human.chm13.genome"
else
    REF_FASTA_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
    GENOME_FILE="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/human.hg38.genome"
fi

WHATSHAP_CONTAINER_SIF_PATH="/home/itoyu8/singularity/scarpia-python_0.2.0.sif"
SAMTOOLS="/home/itoyu8/bin/samtools/samtools-1.19/samtools"
THREADS=${SLURM_CPUS_PER_TASK:-16}

# Create output directory and log directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p log

# Define output files
HAPLOTAG_BAM="${OUTPUT_DIR}/${FILE_PREFIX}.hptag.bam"
HAPLOTAG_TSV="${OUTPUT_DIR}/${FILE_PREFIX}.haplotype.tsv.gz"
H1_BAM="${OUTPUT_DIR}/${FILE_PREFIX}.h1.bam"
H2_BAM="${OUTPUT_DIR}/${FILE_PREFIX}.h2.bam"

# Index phased VCF if needed
if [ ! -f "${PHASED_VCF}.tbi" ]; then
    tabix -p vcf "${PHASED_VCF}"
fi

# Step 0: Convert CRAM to BAM if input is CRAM
if [ "$FILE_EXT" = "cram" ]; then
    echo "Step 0: Converting CRAM to BAM..."
    STEP0_START=$(date +%s)

    CONVERTED_BAM="${INPUT_DIR}/${FILE_BASENAME}.bam"

    "${SAMTOOLS}" view -@ "${THREADS}" -b -T "${REF_FASTA_PATH}" -o "${CONVERTED_BAM}" "${INPUT_FILE}"
    "${SAMTOOLS}" index -@ "${THREADS}" "${CONVERTED_BAM}"

    STEP0_END=$(date +%s)
    STEP0_ELAPSED=$((STEP0_END - STEP0_START))
    echo "Step 0 completed: ${STEP0_ELAPSED} seconds"

    PROCESSING_FILE="${CONVERTED_BAM}"
else
    PROCESSING_FILE="${INPUT_FILE}"
    STEP0_ELAPSED=0
fi

# Step 1: Haplotagging
echo "Step 1: Running whatshap haplotag..."
STEP1_START=$(date +%s)

singularity exec --bind /lustre1/:/lustre1/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap haplotag \
    "${PHASED_VCF}" \
    "${PROCESSING_FILE}" \
    -r "${REF_FASTA_PATH}" \
    -o "${HAPLOTAG_BAM}" \
    --output-threads ${THREADS} \
    --output-haplotag-list "${HAPLOTAG_TSV}" \
    --ignore-read-groups \
    --tag-supplementary \
    --skip-missing-contigs

# Index haplotag BAM
"${SAMTOOLS}" index "${HAPLOTAG_BAM}"

STEP1_END=$(date +%s)
STEP1_ELAPSED=$((STEP1_END - STEP1_START))
echo "Step 1 completed: ${STEP1_ELAPSED} seconds"

# Step 2: Decompress TSV for splitting
echo "Step 2: Decompressing TSV..."
STEP2_START=$(date +%s)

gunzip -f "${HAPLOTAG_TSV}"
HAPLOTAG_TSV_UNZIPPED="${HAPLOTAG_TSV%.gz}"

STEP2_END=$(date +%s)
STEP2_ELAPSED=$((STEP2_END - STEP2_START))
echo "Step 2 completed: ${STEP2_ELAPSED} seconds"

# Step 3: Split BAM by haplotype
echo "Step 3: Splitting BAM by haplotype..."
STEP3_START=$(date +%s)

singularity exec --bind /lustre1/:/lustre1/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap split \
    --output-h1 "${H1_BAM}" \
    --output-h2 "${H2_BAM}" \
    "${HAPLOTAG_BAM}" \
    "${HAPLOTAG_TSV_UNZIPPED}"

# Index split BAM files
"${SAMTOOLS}" index "${H1_BAM}"
"${SAMTOOLS}" index "${H2_BAM}"

STEP3_END=$(date +%s)
STEP3_ELAPSED=$((STEP3_END - STEP3_START))
echo "Step 3 completed: ${STEP3_ELAPSED} seconds"

# Clean up decompressed TSV
rm -f "${HAPLOTAG_TSV_UNZIPPED}"

# End time and summary
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=== Processing Summary ==="
if [ "$FILE_EXT" = "cram" ]; then
    echo "Step 0 (CRAM to BAM conversion): ${STEP0_ELAPSED} seconds"
fi
echo "Step 1 (Haplotagging): ${STEP1_ELAPSED} seconds"
echo "Step 2 (Decompressing TSV): ${STEP2_ELAPSED} seconds"
echo "Step 3 (Splitting BAM): ${STEP3_ELAPSED} seconds"
echo "Total time: ${TOTAL_ELAPSED} seconds"
echo ""
echo "Pipeline completed successfully"
echo "Output directory: ${OUTPUT_DIR}"