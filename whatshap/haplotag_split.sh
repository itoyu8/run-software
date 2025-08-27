#!/bin/bash
#SBATCH -J haplotag_split
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=16G
#SBATCH -c 8

# Usage: sbatch haplotag_split.sh <PHASED_VCF> <BAM_FILE> [output_base_name] [file_prefix]

PHASED_VCF=$1
BAM_FILE=$2
OUTPUT_BASE_NAME=${3:-"whatshap_output"}
FILE_PREFIX=${4:-"sample"}

# Set up output directory in the same directory as BAM file
BAM_DIR=$(dirname "$BAM_FILE")
OUTPUT_DIR="${BAM_DIR}/${OUTPUT_BASE_NAME}"

# Fixed paths
REF_FASTA_PATH="/home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta"
WHATSHAP_CONTAINER_SIF_PATH="/home/itoyu8/singularity/scarpia-python_0.2.0.sif"
GENOME_FILE="/home/itoyu8/database/reference/hg38/v0/human.hg38.genome"

# Create output directory and log directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p log

# Define output files
HAPLOTAG_BAM="${OUTPUT_DIR}/${FILE_PREFIX}.haplotag.bam"
HAPLOTAG_TSV="${OUTPUT_DIR}/${FILE_PREFIX}.haplotype.tsv.gz"
H1_BAM="${OUTPUT_DIR}/${FILE_PREFIX}.h1.bam"
H2_BAM="${OUTPUT_DIR}/${FILE_PREFIX}.h2.bam"

# Index phased VCF if needed
if [ ! -f "${PHASED_VCF}.tbi" ]; then
    tabix -p vcf "${PHASED_VCF}"
fi

# Step 1: Haplotagging
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

# Step 2: Decompress TSV for splitting
gunzip -f "${HAPLOTAG_TSV}"
HAPLOTAG_TSV_UNZIPPED="${HAPLOTAG_TSV%.gz}"

# Step 3: Split BAM by haplotype
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap split \
    --output-h1 "${H1_BAM}" \
    --output-h2 "${H2_BAM}" \
    "${HAPLOTAG_BAM}" \
    "${HAPLOTAG_TSV_UNZIPPED}"

# Index split BAM files
/home/itoyu8/bin/samtools/samtools-1.19/samtools index "${H1_BAM}"
/home/itoyu8/bin/samtools/samtools-1.19/samtools index "${H2_BAM}"

# Clean up decompressed TSV
rm -f "${HAPLOTAG_TSV_UNZIPPED}"