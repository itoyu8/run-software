#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J whtag_split
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1

# Usage: sbatch whtag_split.sh <PHASED_VCF> <BAM_FILE>

# whatshap does not allow multi-threading system (except polyphase)

# Start time
START_TIME=$(date +%s)

PHASED_VCF=$1
BAM_FILE=$2

# Extract BAM file basename without extension
BAM_BASENAME=$(basename "$BAM_FILE" .bam)

# Set up output directory in the same directory as BAM file
BAM_DIR=$(dirname "$BAM_FILE")
OUTPUT_DIR="${BAM_DIR}/${BAM_BASENAME}_hptag"
FILE_PREFIX="${BAM_BASENAME}"

# Fixed paths
REF_FASTA_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
WHATSHAP_CONTAINER_SIF_PATH="/home/itoyu8/singularity/scarpia-python_0.2.0.sif"
GENOME_FILE="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/human.hg38.genome"

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

# Step 1: Haplotagging
echo "Step 1: Running whatshap haplotag..."
STEP1_START=$(date +%s)

singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap haplotag \
    "${PHASED_VCF}" \
    "${BAM_FILE}" \
    -r "${REF_FASTA_PATH}" \
    -o "${HAPLOTAG_BAM}" \
    --output-threads ${NSLOTS:-8} \
    --output-haplotag-list "${HAPLOTAG_TSV}" \
    --ignore-read-groups \
    --tag-supplementary \
    --skip-missing-contigs

# Index haplotag BAM
/home/itoyu8/bin/samtools/samtools-1.19/samtools index "${HAPLOTAG_BAM}"

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

singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap split \
    --output-h1 "${H1_BAM}" \
    --output-h2 "${H2_BAM}" \
    "${HAPLOTAG_BAM}" \
    "${HAPLOTAG_TSV_UNZIPPED}"

# Index split BAM files
/home/itoyu8/bin/samtools/samtools-1.19/samtools index "${H1_BAM}"
/home/itoyu8/bin/samtools/samtools-1.19/samtools index "${H2_BAM}"

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
echo "Step 1 (Haplotagging): ${STEP1_ELAPSED} seconds"
echo "Step 2 (Decompressing TSV): ${STEP2_ELAPSED} seconds"
echo "Step 3 (Splitting BAM): ${STEP3_ELAPSED} seconds"
echo "Total time: ${TOTAL_ELAPSED} seconds"
echo ""
echo "Pipeline completed successfully"
echo "Output directory: ${OUTPUT_DIR}"