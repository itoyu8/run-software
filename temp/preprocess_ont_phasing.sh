#!/bin/bash
#SBATCH -J ont_phasing
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G # Request 4GB RAM per core
#SBATCH -c 16             # Request 56 cores (matches THREADS)
# Exit immediately if a command exits with a non-zero status.

set -e

# Check for NORMAL_BAM_PATH argument
if [ -z "$1" ]; then
    echo "ERROR: NORMAL_BAM_PATH (input BAM file) not provided."
    echo "Usage: $0 <path_to_normal_bam>"
    exit 1
fi
NORMAL_BAM_PATH="$1"
echo "INFO: NORMAL_BAM_PATH set to: ${NORMAL_BAM_PATH}"

# --- Configuration ---
echo "INFO: Configuring script..."
echo "INFO: Please ensure all paths below are correctly set for your environment."

DEEPVARIANT_SIF_PATH="/home/itoyu8/singularity/deepvariant_1.9.0.sif"
WHATSHAP_CONTAINER_SIF_PATH="/home/itoyu8/singularity/scarpia-python_0.2.0.sif" # For whatshap

REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta"
# NORMAL_BAM_PATH is now set from the first command-line argument
BCFTOOLS_PATH="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"
GNOMAD_VCF_PATH="/home/itoyu8/database/reference/gnomAD_4.1/gnomad.genomes.v4.1.sites.merged.light.vcf.bgz"

HOST_OUTPUT_PARENT_DIR="$(pwd)" # Output files will be relative to the current directory
MAIN_OUTPUT_DIR_NAME="ont_phasing_output" # Main subdirectory for outputs
HOST_MAIN_OUTPUT_DIR="${HOST_OUTPUT_PARENT_DIR}/${MAIN_OUTPUT_DIR_NAME}"

BIND_PATH_HOST_1="/home/itoyu8" # For Singularity
BIND_PATH_HOST_2="/lustre1"   # For Singularity

THREADS=16 # For DeepVariant sharding

OUTPUT_VCF_FILENAME="normal_dv.vcf.gz" # Output from DeepVariant
OUTPUT_VCF_PATH="${HOST_MAIN_OUTPUT_DIR}/${OUTPUT_VCF_FILENAME}"
FILTERED_VCF_FILENAME="normal_dv.filtered.vcf.gz" # Output from bcftools
FILTERED_VCF_PATH="${HOST_MAIN_OUTPUT_DIR}/${FILTERED_VCF_FILENAME}"

WHATSHAP_OUTPUT_SUBDIR_NAME="wh_phase_results" # Subdirectory for Whatshap outputs
WHATSHAP_OUTPUT_DIR_ON_HOST="${HOST_MAIN_OUTPUT_DIR}/${WHATSHAP_OUTPUT_SUBDIR_NAME}"
PHASED_VCF_FILENAME="phased.snp.vcf.gz" # Output from Whatshap
PHASED_VCF_OUTPUT_PATH="${WHATSHAP_OUTPUT_DIR_ON_HOST}/${PHASED_VCF_FILENAME}"

# --- Script Start ---
echo "----------------------------------------"
echo "Starting ONT Phasing Script"
echo "----------------------------------------"
echo "Configuration Summary:"
echo "  Input Normal BAM: ${NORMAL_BAM_PATH}"
echo "  DeepVariant SIF: ${DEEPVARIANT_SIF_PATH}"
echo "  Whatshap SIF: ${WHATSHAP_CONTAINER_SIF_PATH}"
echo "  Reference Genome: ${REFERENCE_GENOME_PATH}"
echo "  BCFtools Path: ${BCFTOOLS_PATH}"
echo "  GnomAD VCF Path: ${GNOMAD_VCF_PATH}"
echo "  Main Output Directory (Host): ${HOST_MAIN_OUTPUT_DIR}"
echo "  DeepVariant Output VCF (Host): ${OUTPUT_VCF_PATH}"
echo "  Filtered VCF (Host): ${FILTERED_VCF_PATH}"
echo "  Whatshap Output Dir (Host): ${WHATSHAP_OUTPUT_DIR_ON_HOST}"
echo "  Phased VCF (Host): ${PHASED_VCF_OUTPUT_PATH}"
echo "  Threads/Shards (for DeepVariant): ${THREADS}"
echo "  Host Bind Path 1: ${BIND_PATH_HOST_1}"
echo "  Host Bind Path 2: ${BIND_PATH_HOST_2}"
# Add more bind paths to summary if configured
echo "----------------------------------------"

# Create output directories on the host
echo "INFO: Creating output directories..."
mkdir -p "${HOST_MAIN_OUTPUT_DIR}"
mkdir -p "${WHATSHAP_OUTPUT_DIR_ON_HOST}"
echo "INFO: Output directories created:"
echo "  Main: ${HOST_MAIN_OUTPUT_DIR}"
echo "  Whatshap Output: ${WHATSHAP_OUTPUT_DIR_ON_HOST}"
echo "----------------------------------------"

# --- 1. Run DeepVariant ---
echo "Step 1: Running DeepVariant for variant calling..."

echo "INFO: Using fixed bind paths for DeepVariant: /home/itoyu8/ and /lustre1/"
echo "INFO: Executing DeepVariant Singularity command..."
singularity run --nv \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1/ \
    "${DEEPVARIANT_SIF_PATH}" run_deepvariant \
    --model_type ONT_R104 \
    --ref "${REFERENCE_GENOME_PATH}" \
    --reads "${NORMAL_BAM_PATH}" \
    --output_vcf "${OUTPUT_VCF_PATH}" \
    --num_shards "${THREADS}" \
    || { echo "ERROR: DeepVariant failed. Check logs."; exit 1; }

echo "DeepVariant finished successfully."
echo "Raw Output VCF from DeepVariant: ${OUTPUT_VCF_PATH}"
echo "----------------------------------------"

# --- 2. Filter VCF with bcftools ---
echo "Step 2: Filtering VCF for known SNPs using bcftools..."
echo "INFO: Input VCF for filtering: ${OUTPUT_VCF_PATH}"
echo "INFO: GnomAD VCF for filtering: ${GNOMAD_VCF_PATH}"
echo "INFO: Output filtered VCF: ${FILTERED_VCF_PATH}"

# Command to filter OUTPUT_VCF_PATH to SNPs also present in GNOMAD_VCF_PATH
# bcftools isec options:
# -n=2 : output records present in exactly two files (i.e., the intersection)
# -w2  : write records from the second input file (normal_dv.vcf.gz)
# -v snps : only include SNP variants
# -c all : If a record is present in multiple files, collapse them. Here, with -n=2, it ensures we get the record from the second file.
# -Oz  : output compressed VCF
"${BCFTOOLS_PATH}" isec -n=2 -w2 -c snps -O z -o "${FILTERED_VCF_PATH}" "${GNOMAD_VCF_PATH}" "${OUTPUT_VCF_PATH}" \
    || { echo "ERROR: bcftools filtering failed. Check logs."; exit 1; }

echo "bcftools filtering finished successfully."
echo "Filtered VCF: ${FILTERED_VCF_PATH}"
echo "----------------------------------------"

# --- 3. Run Whatshap ---
echo "Step 3: Running Whatshap for phasing..."
echo "INFO: Input VCF for Whatshap: ${FILTERED_VCF_PATH}"
echo "INFO: Input BAM for Whatshap: ${NORMAL_BAM_PATH}" # This is from the command line arg
echo "INFO: Reference Genome: ${REFERENCE_GENOME_PATH}"
echo "INFO: Output Phased VCF: ${PHASED_VCF_OUTPUT_PATH}"
echo "INFO: Whatshap Output Directory: ${WHATSHAP_OUTPUT_DIR_ON_HOST}"

echo "INFO: Executing Whatshap Singularity command..."
singularity exec --nv \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1/ \
    "${WHATSHAP_CONTAINER_SIF_PATH}" whatshap phase \
    --reference "${REFERENCE_GENOME_PATH}" \
    --ignore-read-groups \
    --distrust-genotypes \
    -o "${PHASED_VCF_OUTPUT_PATH}" \
    "${FILTERED_VCF_PATH}" \
    "${NORMAL_BAM_PATH}" \
    || { echo "ERROR: Whatshap failed. Check logs."; exit 1; }

echo "Whatshap finished successfully."
echo "Phased VCF Output: ${PHASED_VCF_OUTPUT_PATH}"
echo "Output Directory: ${WHATSHAP_OUTPUT_DIR_ON_HOST}"
echo "----------------------------------------"

echo "ONT Phasing Script finished successfully."
