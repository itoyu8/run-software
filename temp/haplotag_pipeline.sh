#!/bin/bash
#SBATCH -J script_job
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=12G # Request 8GB RAM per core
#SBATCH -c 8             # Request 16 cores

# --- Script Configuration & Input Arguments ---
set -e # Exit immediately if a command exits with a non-zero status.

if [ "$#" -ne 5 ]; then # Require exactly 5 arguments
    echo "Usage: qsub $0 <GVCF_PATH> <NORMAL_MARKDUP_PATH> <TUMOR_MARKDUP_PATH> <CHR_NUM> <OUTPUT_DIR_NAME>"
    exit 1
fi

# --- Fixed Paths (Edit these directly if defaults change) ---
# Reference files are now fixed within the script
REF_FASTA_PATH="/home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta"
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
GVCF_PATH=$1
NORMAL_MARKDUP_PATH=$2
TUMOR_MARKDUP_PATH=$3
CHR_NUM=$4
OUTPUT_DIR_NAME=$5

# Check if any arguments are empty (basic check)
if [ -z "$GVCF_PATH" ] || [ -z "$NORMAL_MARKDUP_PATH" ] || [ -z "$TUMOR_MARKDUP_PATH" ] || [ -z "$CHR_NUM" ] || [ -z "$OUTPUT_DIR_NAME" ]; then
    echo "Error: One or more required arguments are empty."
    echo "Usage: qsub $0 <GVCF_PATH> <NORMAL_MARKDUP_PATH> <TUMOR_MARKDUP_PATH> <CHR_NUM> <OUTPUT_DIR_NAME>"
    exit 1
fi

# Define output directory path relative to the current working directory
OUTPUT_DIR="./${OUTPUT_DIR_NAME}"

echo "--- Configuration ---"
echo "Input GVCF: ${GVCF_PATH}"
echo "Normal Markdup: ${NORMAL_MARKDUP_PATH}"
echo "Tumor Markdup: ${TUMOR_MARKDUP_PATH}"
echo "Chromosome: ${CHR_NUM}"
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

# Define intermediate file paths in the Current Working Directory (CWD, where qsub is run)
RAW_VCF="./genotyped.vcf.gz"
SNP_VCF="./genotyped.snp.vcf.gz"
FILT_SNP_VCF="./genotyped.snp.filt.vcf.gz"
PASS_FILT_SNP_VCF="./genotyped.snp.filt.pass.vcf.gz" # This is used by Step 2

# --- Step 1: VCF Filtering (Conditional Execution) ---
# Check if the final output of Step 1 already exists
if [ -f "$PASS_FILT_SNP_VCF" ] && [ -f "${PASS_FILT_SNP_VCF}.tbi" ]; then
    echo "Step 1: Found existing ${PASS_FILT_SNP_VCF} and index. Skipping VCF filtering."
else
    echo "Step 1: Filtering VCF (Intermediate files in CWD)..."

    # --- Check Input GVCF and Index ---
    echo "Checking input GVCF format and index..."
    if [[ ! "${GVCF_PATH}" == *.g.vcf.gz ]]; then
        echo "Error: Input GVCF file '${GVCF_PATH}' does not end with .g.vcf.gz"
        exit 1
    fi

    GVCF_INDEX_PATH="${GVCF_PATH}.tbi"
    if [ ! -f "${GVCF_INDEX_PATH}" ]; then
        echo "Warning: Index file '${GVCF_INDEX_PATH}' not found for input GVCF."
        echo "Generating index using tabix (via Apptainer)..."
        # Only need /home/itoyu8 bind as GVCF path is assumed to be under it
        singularity exec \
          --bind /home/itoyu8/:/home/itoyu8/ \
          "${APPTAINER_SIF_PATH}" tabix -p vcf "${GVCF_PATH}"
        echo "Index generated."
    else
        echo "Index file found: ${GVCF_INDEX_PATH}"
    fi

    echo "Running GenotypeGVCFs..."
    # Output RAW_VCF to CWD
    singularity exec \
      --bind /home/itoyu8/:/home/itoyu8/ \
      "${APPTAINER_SIF_PATH}" /usr/bin/java \
      -Xmx4G -jar "${GATK_JAR_PATH}" GenotypeGVCFs \
      -R "${REF_FASTA_PATH}" \
      -V "${GVCF_PATH}" \
      -O "${RAW_VCF}"

    # echo "Indexing raw VCF (using Apptainer)..." # Commented out as GATK often indexes

    echo "Running SelectVariants (SNP)..."
    # Output SNP_VCF to CWD
    singularity exec \
      --bind /home/itoyu8/:/home/itoyu8/ \
      "${APPTAINER_SIF_PATH}" /usr/bin/java \
      -Xmx4G -jar "${GATK_JAR_PATH}" SelectVariants \
        -V "${RAW_VCF}" \
        -select-type SNP \
        -O "${SNP_VCF}"

    echo "Running VariantFiltration..."
    # Output FILT_SNP_VCF to CWD
    singularity exec \
      --bind /home/itoyu8/:/home/itoyu8/ \
      "${APPTAINER_SIF_PATH}" /usr/bin/java \
        -Xmx4G -jar "${GATK_JAR_PATH}" VariantFiltration \
        -V "${SNP_VCF}" \
        -O "${FILT_SNP_VCF}" \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

    echo "Filtering for PASS variants with bcftools..."
    # Output PASS_FILT_SNP_VCF to CWD
    "${BCFTOOLS_PATH}" view -f PASS -O z "${FILT_SNP_VCF}" > "${PASS_FILT_SNP_VCF}"

    echo "Indexing PASS filtered VCF (using Apptainer)..."
    # Index PASS_FILT_SNP_VCF in CWD
    singularity exec \
      --bind /home/itoyu8/:/home/itoyu8/ \
      "${APPTAINER_SIF_PATH}" tabix -p vcf "${PASS_FILT_SNP_VCF}" # Assumes tabix is in container PATH

    echo "Step 1 finished."
fi # End of Step 1 conditional execution

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

# --- Step 2: Phasing (from beagle.sh) ---
echo "Step 2: Phasing VCF with Beagle..."

# Input PASS_FILT_SNP_VCF is now expected in CWD (./genotyped.snp.filt.pass.vcf.gz)
# Output PHASED_VCF goes to OUTPUT_DIR
PHASED_VCF_PREFIX="${OUTPUT_DIR}/phased.chr${CHR_NUM}"
PHASED_VCF="${PHASED_VCF_PREFIX}.vcf.gz"
REF_1000G_VCF="${REF_1000G_PATH_PREFIX}chr${CHR_NUM}.filtered.shapeit2-duohmm-phased.vcf.gz" # Construct full path

# --- Check if reference VCF exists and Input VCF from Step 1 exists ---
[ -f "$REF_1000G_VCF" ] || { echo "Error: 1000G Reference VCF not found at $REF_1000G_VCF"; exit 1; }
[ -f "$PASS_FILT_SNP_VCF" ] || { echo "Error: PASS filtered SNP VCF ${PASS_FILT_SNP_VCF} not found in CWD"; exit 1; }

echo "Running Beagle..."
# Beagle runs outside container, needs access to CWD for input VCF and OUTPUT_DIR for output
java -Xmx24g -jar "${BEAGLE_JAR_PATH}" \
    gt="${PASS_FILT_SNP_VCF}" \
    ref="${REF_1000G_VCF}" \
    out="${PHASED_VCF_PREFIX}" \
    chrom="chr${CHR_NUM}" \
    nthreads=4 # Fixed number of threads for Beagle

echo "Indexing phased VCF (using Apptainer)..."
# Index PHASED_VCF in OUTPUT_DIR
singularity exec \
  --bind /home/itoyu8/:/home/itoyu8/ \
  "${APPTAINER_SIF_PATH}" tabix -p vcf "${PHASED_VCF}" # Assumes tabix is in container PATH

echo "Step 2 finished."

# --- Step 3: Haplotagging (from split.sh) ---
echo "Step 3: Haplotagging BAM files..."

# Input PHASED_VCF is in OUTPUT_DIR
# Outputs go to OUTPUT_DIR
NORMAL_HAPLOTAG_BAM="${OUTPUT_DIR}/normal.haplotag.chr${CHR_NUM}.bam"
NORMAL_HAPLOTAG_TSV="${OUTPUT_DIR}/normal.haplotag.chr${CHR_NUM}.haplotype.tsv.gz"
TUMOR_HAPLOTAG_BAM="${OUTPUT_DIR}/tumor.haplotag.chr${CHR_NUM}.bam"
TUMOR_HAPLOTAG_TSV="${OUTPUT_DIR}/tumor.haplotag.chr${CHR_NUM}.haplotype.tsv.gz"

echo "Haplotagging Normal BAM..."
# Run whatshap inside its dedicated container
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap haplotag \
    "${PHASED_VCF}" \
    "${NORMAL_MARKDUP_PATH}" \
    -r "${REF_FASTA_PATH}" \
    --regions "chr${CHR_NUM}" \
    -o "${NORMAL_HAPLOTAG_BAM}" \
    --output-threads ${NSLOTS:-8} \
    --output-haplotag-list "${NORMAL_HAPLOTAG_TSV}" \
    --ignore-read-groups

echo "Haplotagging Tumor BAM..."
# Run whatshap inside its dedicated container
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap haplotag \
    "${PHASED_VCF}" \
    "${TUMOR_MARKDUP_PATH}" \
    -r "${REF_FASTA_PATH}" \
    --regions "chr${CHR_NUM}" \
    -o "${TUMOR_HAPLOTAG_BAM}" \
    --output-threads ${NSLOTS:-8} \
    --output-haplotag-list "${TUMOR_HAPLOTAG_TSV}" \
    --ignore-read-groups

echo "Indexing haplotagged BAM files (using Apptainer and samtools)..."
# Index BAMs in OUTPUT_DIR
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${NORMAL_HAPLOTAG_BAM}"
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${TUMOR_HAPLOTAG_BAM}"
echo "Haplotagged BAM indexing complete."

echo "Decompressing Normal haplotype TSV for splitting..."
gunzip -f "${NORMAL_HAPLOTAG_TSV}" # TSV is in OUTPUT_DIR
NORMAL_HAPLOTAG_TSV_UNZIPPED="${NORMAL_HAPLOTAG_TSV%.gz}" # Remove .gz suffix

echo "Splitting Normal BAM by haplotype..."
# Run whatshap inside its dedicated container
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap split \
    --output-h1 "${OUTPUT_DIR}/normal.haplotag.chr${CHR_NUM}.h1.bam" \
    --output-h2 "${OUTPUT_DIR}/normal.haplotag.chr${CHR_NUM}.h2.bam" \
    "${NORMAL_HAPLOTAG_BAM}" \
    "${NORMAL_HAPLOTAG_TSV_UNZIPPED}" # Use the unzipped TSV file

echo "Decompressing Tumor haplotype TSV for splitting..."
gunzip -f "${TUMOR_HAPLOTAG_TSV}" # TSV is in OUTPUT_DIR
TUMOR_HAPLOTAG_TSV_UNZIPPED="${TUMOR_HAPLOTAG_TSV%.gz}" # Remove .gz suffix

echo "Splitting Tumor BAM by haplotype..."
# Run whatshap inside its dedicated container
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap split \
    --output-h1 "${OUTPUT_DIR}/tumor.haplotag.chr${CHR_NUM}.h1.bam" \
    --output-h2 "${OUTPUT_DIR}/tumor.haplotag.chr${CHR_NUM}.h2.bam" \
    "${TUMOR_HAPLOTAG_BAM}" \
    "${TUMOR_HAPLOTAG_TSV_UNZIPPED}" # Use the unzipped TSV file

echo "Indexing split BAM files (using Apptainer and samtools)..."
# Define split BAM paths in OUTPUT_DIR
NORMAL_H1_BAM="${OUTPUT_DIR}/normal.haplotag.chr${CHR_NUM}.h1.bam"
NORMAL_H2_BAM="${OUTPUT_DIR}/normal.haplotag.chr${CHR_NUM}.h2.bam"
TUMOR_H1_BAM="${OUTPUT_DIR}/tumor.haplotag.chr${CHR_NUM}.h1.bam"
TUMOR_H2_BAM="${OUTPUT_DIR}/tumor.haplotag.chr${CHR_NUM}.h2.bam"

# Index split BAMs in OUTPUT_DIR
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${NORMAL_H1_BAM}"
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${NORMAL_H2_BAM}"
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${TUMOR_H1_BAM}"
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" samtools index "${TUMOR_H2_BAM}"
echo "Split BAM indexing complete."

# Create windows file using bedtools and awk
WINDOWS_FILE="${OUTPUT_DIR}/chr${CHR_NUM}.windows.bed"
# GENOME_FILE is defined in CWD after Step 1
if [ -f "$GENOME_FILE" ]; then
  echo "Creating windows file for chr${CHR_NUM}..."
  # bedtools runs outside container
  # Needs GENOME_FILE (CWD) and output to WINDOWS_FILE (OUTPUT_DIR)
  "${BEDTOOLS_PATH}" makewindows -g "${GENOME_FILE}" -w 10000 | awk -v chr="chr${CHR_NUM}" '$1 == chr' > "${WINDOWS_FILE}" \
    || { echo "Failed to create windows file"; exit 1; }
  echo "Windows file created: ${WINDOWS_FILE}"
else
  echo "Warning: Genome file ${GENOME_FILE} not found, skipping window creation."
fi

# --- Step 4: Calculate Read Counts per Window ---
echo "Step 4: Calculating read counts per window using samtools bedcov..."
READCOUNT_OUTPUT_FILE="${OUTPUT_DIR}/chr${CHR_NUM}_readcount.txt"

# Ensure input files exist before running bedcov
# WINDOWS_FILE and split BAMs are in OUTPUT_DIR
if [ ! -f "$WINDOWS_FILE" ]; then
  echo "Error: Windows file ${WINDOWS_FILE} not found for bedcov step."
  exit 1
fi
if [ ! -f "$TUMOR_H1_BAM" ] || [ ! -f "$NORMAL_H1_BAM" ] || [ ! -f "$TUMOR_H2_BAM" ] || [ ! -f "$NORMAL_H2_BAM" ]; then
  echo "Error: One or more split BAM files not found for bedcov step."
  exit 1
fi

# Execute bedcov inside the container, redirecting output
# Order: Tumor H1, Normal H1, Tumor H2, Normal H2
echo "Running samtools bedcov (Tumor H1, Normal H1, Tumor H2, Normal H2)..."
# Inputs are in OUTPUT_DIR, Output goes to OUTPUT_DIR
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${APPTAINER_SIF_PATH}" \
  samtools bedcov "${WINDOWS_FILE}" "${TUMOR_H1_BAM}" "${NORMAL_H1_BAM}" "${TUMOR_H2_BAM}" "${NORMAL_H2_BAM}" \
  > "${READCOUNT_OUTPUT_FILE}" \
  || { echo "Failed to run samtools bedcov"; exit 1; }

echo "Read count file created: ${READCOUNT_OUTPUT_FILE}"
echo "Step 4 finished."

# Final completion message
echo "--- Pipeline Completed Successfully ---"
echo "All steps (1-4) completed."

echo "Cleaning up intermediate files from CWD..."
# Remove intermediate VCF files from Step 1 (now in CWD)
# Keep PASS_FILT_SNP_VCF and its index as they are reusable
rm -f "${RAW_VCF}" "${RAW_VCF}.tbi" \
      "${SNP_VCF}" "${SNP_VCF}.tbi" \
      "${FILT_SNP_VCF}" "${FILT_SNP_VCF}.tbi"
echo "Intermediate file cleanup finished."

echo "Final output files are in: ${OUTPUT_DIR}"
