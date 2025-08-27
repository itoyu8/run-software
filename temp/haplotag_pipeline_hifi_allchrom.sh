#!/bin/bash
#SBATCH -J haplotag_pipeline_hifi_allchrom
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G # Request 8GB RAM per core
#SBATCH -c 8             # Request 8 cores

# --- Script Configuration & Input Arguments ---
set -e # Exit immediately if a command exits with a non-zero status.

if [ "$#" -ne 4 ]; then # Require exactly 4 arguments now (CHR_NUM removed)
    echo "Usage: qsub $0 <PHASED_VCF> <NORMAL_MARKDUP_PATH> <TUMOR_MARKDUP_PATH> <OUTPUT_DIR_NAME>"
    exit 1
fi

# --- Fixed Paths (Edit these directly if defaults change) ---
# Reference files are now fixed within the script
#REF_FASTA_PATH="/home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta"
REF_FASTA_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
REF_1000G_PATH_PREFIX="/home/itoyu8/database/1000genomes/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_"

# Tool Paths
BEAGLE_JAR_PATH="/home/itoyu8/bin/beagle/beagle.17Dec24.224.jar"
APPTAINER_SIF_PATH="/home/itoyu8/singularity/compat_parabricks-0.2.2.sif" # For GATK, samtools, etc.
WHATSHAP_CONTAINER_SIF_PATH="/home/itoyu8/singularity/scarpia-python_0.2.0.sif" # For whatshap
GATK_JAR_PATH="/tools/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar" # Assuming this path inside container is fixed
# WHATSHAP_PATH="/home/itoyu8/bin/whatshap/whatshap" # No longer needed for direct execution
BEDTOOLS_PATH="/home/itoyu8/bin/bedtools/bedtools2/bin/bedtools" # Note: Bedtools module might provide this too
BCFTOOLS_PATH="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools" # Point directly to the executable

# --- Input Argument Processing ---
PHASED_VCF=$1
NORMAL_MARKDUP_PATH=$2
TUMOR_MARKDUP_PATH=$3
OUTPUT_DIR_NAME=$4
# CHR_NUM is no longer an argument

# Check if any arguments are empty (basic check)
if [ -z "$PHASED_VCF" ] || [ -z "$NORMAL_MARKDUP_PATH" ] || [ -z "$TUMOR_MARKDUP_PATH" ] || [ -z "$OUTPUT_DIR_NAME" ]; then
    echo "Error: One or more required arguments are empty."
    echo "Usage: qsub $0 <PHASED_VCF> <NORMAL_MARKDUP_PATH> <TUMOR_MARKDUP_PATH> <OUTPUT_DIR_NAME>"
    exit 1
fi

# Define the list of chromosomes to process
CHROMS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
# CHROMS=("chr22") # For testing

# Define output directory path relative to the current working directory
OUTPUT_DIR="./${OUTPUT_DIR_NAME}"

echo "--- Configuration ---"
echo "Input phased VCF: ${PHASED_VCF}"
echo "Normal Markdup: ${NORMAL_MARKDUP_PATH}"
echo "Tumor Markdup: ${TUMOR_MARKDUP_PATH}"
echo "Chromosomes to process: ${CHROMS[*]}"
echo "Output Directory (Final files): ${OUTPUT_DIR}" # Clarified output dir usage
echo "Reference FASTA (fixed): ${REF_FASTA_PATH}"
echo "1000G Ref Prefix (fixed): ${REF_1000G_PATH_PREFIX}"
echo "GATK JAR (in container): ${GATK_JAR_PATH}"
echo "Beagle JAR: ${BEAGLE_JAR_PATH}"
echo "Apptainer SIF (GATK, etc.): ${APPTAINER_SIF_PATH}"
echo "Whatshap Container SIF: ${WHATSHAP_CONTAINER_SIF_PATH}" # Updated to show container path
echo "Bedtools Path: ${BEDTOOLS_PATH}" # Module might override this
echo "Bcftools Path: ${BCFTOOLS_PATH}"
echo "---------------------"

# --- Setup ---
echo "Ensuring output directory exists: ${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"
echo "Ensuring log directory exists: ./log/" # Log dir relative to script CWD
mkdir -p ./log/


echo "Loading required modules..."
# module use /usr/local/package/modulefiles/
# Load modules that provide dependencies like java, potentially samtools, bedtools if not using explicit path
# module load apptainer samtools bedtools # Keep loading modules as they might set up environment/other deps

echo "Verifying tool paths..."
# Optional: Add checks to ensure the executables exist at the specified paths
[ -f "$BEAGLE_JAR_PATH" ] || { echo "Error: Beagle JAR not found at $BEAGLE_JAR_PATH"; exit 1; }
# [ -x "$WHATSHAP_PATH" ] || { echo "Error: Whatshap executable not found or not executable at $WHATSHAP_PATH"; exit 1; } # Commented out as whatshap runs in container
[ -x "$BCFTOOLS_PATH" ] || { echo "Error: Bcftools executable not found or not executable at $BCFTOOLS_PATH"; exit 1; }
[ -f "$APPTAINER_SIF_PATH" ] || { echo "Error: Apptainer SIF not found at $APPTAINER_SIF_PATH"; exit 1; }
[ -f "$WHATSHAP_CONTAINER_SIF_PATH" ] || { echo "Error: Whatshap Container SIF not found at $WHATSHAP_CONTAINER_SIF_PATH"; exit 1; } # Added check for the whatshap container
# Add more checks as needed (e.g., for reference files)

# --- Generate Genome File if needed ---
GENOME_FILE="./human.hg38.genome" # Genome file in the script's CWD
REF_FASTA_FAI="${REF_FASTA_PATH}.fai"

if [ ! -f "${GENOME_FILE}" ]; then
    echo "Generating genome file: ${GENOME_FILE}..."
    if [ ! -f "${REF_FASTA_FAI}" ]; then
        echo "Error: Reference FASTA index not found: ${REF_FASTA_FAI}"
        echo "Please index the reference FASTA using 'samtools faidx ${REF_FASTA_PATH}'"
        exit 1
    fi
    # Run awk on the host to create the genome file in CWD
    echo "Running awk on host to generate genome file..."
    awk -v OFS='\t' '{print $1, $2}' "${REF_FASTA_FAI}" > "${GENOME_FILE}" \
      || { echo "Failed to generate genome file"; exit 1; }
    echo "Genome file generated: ${GENOME_FILE}"
else
    echo "Genome file already exists: ${GENOME_FILE}"
fi

echo "Indexing phased VCF..."
# Index PHASED_VCF if .tbi file does not exist.
# PHASED_VCF should now be the .vcf.gz file if compression occurred.
# The index file (.tbi) will be created in the same directory as PHASED_VCF.
if [ ! -f "${PHASED_VCF}.tbi" ]; then
    echo "Index file ${PHASED_VCF}.tbi not found. Creating index..."
    tabix -p vcf "${PHASED_VCF}"
    
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create index for ${PHASED_VCF}. Exiting."
        exit 1
    else
        echo "Successfully created index for ${PHASED_VCF}."
    fi
else
    echo "Index file ${PHASED_VCF}.tbi already exists. Skipping indexing."
fi

echo "Step 2 finished."

# --- Step 3: Haplotagging (all chromosomes at once) ---
echo "Step 3: Haplotagging BAM files (all chromosomes)..."

NORMAL_HAPLOTAG_BAM="${OUTPUT_DIR}/normal.haplotag.all_chroms.bam"
NORMAL_HAPLOTAG_TSV="${OUTPUT_DIR}/normal.haplotag.all_chroms.haplotype.tsv.gz"
TUMOR_HAPLOTAG_BAM="${OUTPUT_DIR}/tumor.haplotag.all_chroms.bam"
TUMOR_HAPLOTAG_TSV="${OUTPUT_DIR}/tumor.haplotag.all_chroms.haplotype.tsv.gz"

echo "Haplotagging Normal BAM (all chromosomes)..."
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap haplotag \
    "${PHASED_VCF}" \
    "${NORMAL_MARKDUP_PATH}" \
    -r "${REF_FASTA_PATH}" \
    -o "${NORMAL_HAPLOTAG_BAM}" \
    --output-threads ${NSLOTS:-8} \
    --output-haplotag-list "${NORMAL_HAPLOTAG_TSV}" \
    --ignore-read-groups \
    --tag-supplementary \
    --skip-missing-contigs \

echo "Haplotagging Tumor BAM (all chromosomes)..."
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap haplotag \
    "${PHASED_VCF}" \
    "${TUMOR_MARKDUP_PATH}" \
    -r "${REF_FASTA_PATH}" \
    -o "${TUMOR_HAPLOTAG_BAM}" \
    --output-threads ${NSLOTS:-8} \
    --output-haplotag-list "${TUMOR_HAPLOTAG_TSV}" \
    --ignore-read-groups \
    --tag-supplementary \
    --skip-missing-contigs \

echo "Indexing haplotagged BAM files..."
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${NORMAL_HAPLOTAG_BAM}"
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${TUMOR_HAPLOTAG_BAM}"
echo "Haplotagged BAM indexing complete."

echo "Decompressing Normal haplotype TSV for splitting..."
gunzip -f "${NORMAL_HAPLOTAG_TSV}"
NORMAL_HAPLOTAG_TSV_UNZIPPED="${NORMAL_HAPLOTAG_TSV%.gz}"

echo "Splitting Normal BAM by haplotype (all chromosomes)..."
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap split \
    --output-h1 "${OUTPUT_DIR}/normal.haplotag.all_chroms.h1.bam" \
    --output-h2 "${OUTPUT_DIR}/normal.haplotag.all_chroms.h2.bam" \
    "${NORMAL_HAPLOTAG_BAM}" \
    "${NORMAL_HAPLOTAG_TSV_UNZIPPED}"

echo "Decompressing Tumor haplotype TSV for splitting..."
gunzip -f "${TUMOR_HAPLOTAG_TSV}"
TUMOR_HAPLOTAG_TSV_UNZIPPED="${TUMOR_HAPLOTAG_TSV%.gz}"

echo "Splitting Tumor BAM by haplotype (all chromosomes)..."
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap split \
    --output-h1 "${OUTPUT_DIR}/tumor.haplotag.all_chroms.h1.bam" \
    --output-h2 "${OUTPUT_DIR}/tumor.haplotag.all_chroms.h2.bam" \
    "${TUMOR_HAPLOTAG_BAM}" \
    "${TUMOR_HAPLOTAG_TSV_UNZIPPED}"

echo "Indexing split BAM files..."
NORMAL_H1_BAM="${OUTPUT_DIR}/normal.haplotag.all_chroms.h1.bam"
NORMAL_H2_BAM="${OUTPUT_DIR}/normal.haplotag.all_chroms.h2.bam"
TUMOR_H1_BAM="${OUTPUT_DIR}/tumor.haplotag.all_chroms.h1.bam"
TUMOR_H2_BAM="${OUTPUT_DIR}/tumor.haplotag.all_chroms.h2.bam"

singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${NORMAL_H1_BAM}"
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${NORMAL_H2_BAM}"
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${TUMOR_H1_BAM}"
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${TUMOR_H2_BAM}"
echo "Split BAM indexing complete."

# --- Step 4: Create windows file for all chromosomes and calculate read counts ---
echo "Step 4: Creating windows file for all target chromosomes and calculating read counts..."
ALL_CHROMS_WINDOWS_FILE="${OUTPUT_DIR}/all_chroms.windows.bed"

if [ -f "$GENOME_FILE" ]; then
    echo "Creating windows file for all target chromosomes..."
    # Generate 10kb windows for all chromosomes, then filter to only include target chromosomes
    "${BEDTOOLS_PATH}" makewindows -g "${GENOME_FILE}" -w 10000 > "${OUTPUT_DIR}/all_genome_windows.bed" \
      || { echo "Failed to create genome-wide windows file"; exit 1; }
    
    # Filter to only include target chromosomes in the specified order
    > "${ALL_CHROMS_WINDOWS_FILE}" # Create empty file
    for CHR_NAME in "${CHROMS[@]}"; do
        grep "^${CHR_NAME}\s" "${OUTPUT_DIR}/all_genome_windows.bed" >> "${ALL_CHROMS_WINDOWS_FILE}" || true
    done
    
    # Clean up temporary file
    rm -f "${OUTPUT_DIR}/all_genome_windows.bed"
    
    if [ -s "${ALL_CHROMS_WINDOWS_FILE}" ]; then
        echo "Windows file created: ${ALL_CHROMS_WINDOWS_FILE}"
        TOTAL_WINDOWS=$(wc -l < "${ALL_CHROMS_WINDOWS_FILE}")
        echo "Total windows: ${TOTAL_WINDOWS}"
    else
        echo "Error: Windows file is empty after filtering for target chromosomes."
        exit 1
    fi
else
    echo "Error: Genome file ${GENOME_FILE} not found."
    exit 1
fi

echo "Calculating read counts per window using samtools bedcov (all chromosomes)..."
FINAL_READCOUNT_FILE="${OUTPUT_DIR}/all_chroms_readcount.txt"

# Ensure all required BAM files exist
if [ ! -f "$TUMOR_H1_BAM" ] || [ ! -f "$NORMAL_H1_BAM" ] || [ ! -f "$TUMOR_H2_BAM" ] || [ ! -f "$NORMAL_H2_BAM" ]; then
    echo "Error: One or more split BAM files not found for bedcov step."
    exit 1
fi

echo "Running samtools bedcov (all chromosomes: Tumor H1, Normal H1, Tumor H2, Normal H2)..."
/home/itoyu8/bin/samtools/samtools-1.19/samtools bedcov -c "${ALL_CHROMS_WINDOWS_FILE}" "${TUMOR_H1_BAM}" "${NORMAL_H1_BAM}" "${TUMOR_H2_BAM}" "${NORMAL_H2_BAM}" \
  > "${FINAL_READCOUNT_FILE}" \
  || { echo "Failed to run samtools bedcov"; exit 1; }

echo "Read count file created: ${FINAL_READCOUNT_FILE}"
TOTAL_LINES=$(wc -l < "${FINAL_READCOUNT_FILE}")
echo "Total lines in readcount file: ${TOTAL_LINES}"

echo "Step 4 finished."

# --- Step 5: Clean up intermediate files (optional) ---
echo "Step 5: Cleaning up intermediate files..."

echo "Removing intermediate haplotagged BAM files (keeping only split BAMs)..."
rm -f "${NORMAL_HAPLOTAG_BAM}" "${NORMAL_HAPLOTAG_BAM}.bai" \
      "${TUMOR_HAPLOTAG_BAM}" "${TUMOR_HAPLOTAG_BAM}.bai" \
      "${NORMAL_HAPLOTAG_TSV_UNZIPPED}" "${TUMOR_HAPLOTAG_TSV_UNZIPPED}"

echo "Intermediate file cleanup complete."
echo "Step 5 finished."

# Final completion message
echo "--- Pipeline Completed Successfully ---"
echo "All chromosomes processed in batch mode."
echo "Final readcount file: ${FINAL_READCOUNT_FILE}"
echo "Split BAM files:"
echo "  Normal H1: ${NORMAL_H1_BAM}"
echo "  Normal H2: ${NORMAL_H2_BAM}"
echo "  Tumor H1: ${TUMOR_H1_BAM}"
echo "  Tumor H2: ${TUMOR_H2_BAM}"
echo "Windows file: ${ALL_CHROMS_WINDOWS_FILE}"
echo "All output files are in: ${OUTPUT_DIR}"