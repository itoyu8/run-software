#!/bin/bash
#SBATCH -J hifiasm_asssembly
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=7G # Request 6GB RAM per core
#SBATCH -c 56            # Request 56 cores

# --- Script Configuration & Input Arguments ---
set -e # Exit immediately if a command exits with a non-zero status.

# --- Script Information ---
# Usage: ./hifiasm_assembly.sh <phased_vcf_gz> <stats_file> <normal_bam> <ont_fastq_gz>
#
# This script performs:
# 1. Haplotagging on a normal BAM file using a phased VCF.
# 2. Creates a BED file from a statistics file, filtering for regions where 'diff' is True.
# 3. Extracts reads from the haplotagged BAM based on HP tag and the BED file.
# 4. Converts the extracted reads to FASTQ format.
#
# Arguments:
#   <phased_vcf_gz>: Path to the phased and gzipped VCF file (e.g., input.vcf.gz).
#                    The script will attempt to index it if a .tbi file is not found.
#   <stats_file>:    Path to the statistics file.
#                    Example format:
#                    chr_num	start_block_pos	end_block_pos	h1_ratio_mean	h1_ratio_var	h2_ratio_mean	h2_ratio_var	p-value	diff
#                    chr1	10000	45850000	0.145346	0.578499	6.474829	4.309993	0.0	True
#   <normal_bam>:    Path to the sorted normal BAM file.
#   <ont_fastq_gz>:  Path to the ONT FASTQ.gz file for hifiasm assembly.
#
# Outputs (in OUTPUT_DIR):
#   - normal.re_haplotag.bam: Haplotagged BAM file (as named in current script).
#   - normal.re_haplotag.haplotype.tsv.gz: Haplotag list (as named in current script).
#   - regions_to_extract.bed: BED file generated from stats_file.
#   - HP1_filtered.bam: Reads with HP:1 from regions_to_extract.bed.
#   - HP2_filtered.bam: Reads with HP:2 from regions_to_extract.bed.
#   - HP1_filtered.fq.gz: FASTQ from HP1_filtered.bam.
#   - HP2_filtered.fq.gz: FASTQ from HP2_filtered.bam.
#   - hp1.yak: Yak count file for HP1 reads.
#   - hp2.yak: Yak count file for HP2 reads.
#   - hifiasm_assembly.*: Hifiasm assembly output files (e.g., .gfa, .log).
#
# Note: You may need to make this script executable using: chmod +x hifiasm_assembly.sh
#       Ensure WHATSHAP_CONTAINER_SIF_PATH and REF_FASTA_PATH are correctly set below.

# --- Configuration ---
# !!! IMPORTANT: Adjust these paths to your environment !!!
WHATSHAP_CONTAINER_SIF_PATH="/home/itoyu8/singularity/scarpia-python_0.2.0.sif" # For whatshap
REF_FASTA_PATH="/home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta"
YAK_PATH="/home/itoyu8/bin/yak/yak-0.1/yak"
HIFIASM_PATH="/home/itoyu8/bin/hifiasm/hifiasm-0.25.0/hifiasm"
OUTPUT_DIR="./hifiasm_output"                                      # Output directory
NSLOTS=${NSLOTS:-56}                                               # Number of threads for CPU-bound tasks (whatshap, samtools, yak, hifiasm)
SAMTOOLS_PATH="/home/itoyu8/bin/samtools/samtools-1.19/samtools"   # Path to samtools
SAMTOOLS_THREADS=${NSLOTS}                                         # Number of threads for samtools, aligned with NSLOTS
YAK_THREADS=${NSLOTS}                                              # Number of threads for yak, aligned with NSLOTS
HIFIASM_THREADS=${NSLOTS}                                          # Number of threads for hifiasm, aligned with NSLOTS

# --- Argument Parsing ---
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <phased_vcf_gz> <stats_file> <normal_bam> <ont_fastq_gz>"
    exit 1
fi

PHASED_VCF_GZ="$1"
STATS_FILE="$2"
NORMAL_BAM="$3"
ONT_FASTQ_GZ="$4"

echo "--- Input Parameters ---"
echo "Phased VCF (gzipped): ${PHASED_VCF_GZ}"
echo "Stats File: ${STATS_FILE}"
echo "Normal BAM: ${NORMAL_BAM}"
echo "ONT FASTQ for Assembly: ${ONT_FASTQ_GZ}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Whatshap Container: ${WHATSHAP_CONTAINER_SIF_PATH}"
echo "Reference FASTA: ${REF_FASTA_PATH}"
echo "Yak Path: ${YAK_PATH}"
echo "Hifiasm Path: ${HIFIASM_PATH}"
echo "Threads (NSLOTS for CPU-bound): ${NSLOTS}"
echo "Samtools Path: ${SAMTOOLS_PATH}"
echo "Samtools Threads: ${SAMTOOLS_THREADS}"
echo "Yak Threads: ${YAK_THREADS}"
echo "Hifiasm Threads: ${HIFIASM_THREADS}"
echo "------------------------"

# --- Preliminary Checks ---
if [ ! -x "${SAMTOOLS_PATH}" ]; then
    echo "Error: Samtools not found or not executable at ${SAMTOOLS_PATH}. Exiting."
fi
if [ ! -x "${YAK_PATH}" ]; then
    echo "Error: yak not found or not executable at ${YAK_PATH}. Exiting."
fi
if [ ! -x "${HIFIASM_PATH}" ]; then
    echo "Error: hifiasm not found or not executable at ${HIFIASM_PATH}. Exiting."
fi
if ! command -v awk &> /dev/null; then
    echo "Error: awk command could not be found. Exiting."
fi
if [ ! -f "${STATS_FILE}" ] || [ ! -r "${STATS_FILE}" ]; then
    echo "Error: Stats file ${STATS_FILE} not found or not readable. Exiting."
fi
if [ ! -f "${PHASED_VCF_GZ}" ] || [ ! -r "${PHASED_VCF_GZ}" ]; then
    echo "Error: Phased VCF ${PHASED_VCF_GZ} not found or not readable. Exiting."
fi
if [ ! -f "${NORMAL_BAM}" ] || [ ! -r "${NORMAL_BAM}" ]; then
    echo "Error: Normal BAM ${NORMAL_BAM} not found or not readable. Exiting."
fi
if [ ! -f "${ONT_FASTQ_GZ}" ] || [ ! -r "${ONT_FASTQ_GZ}" ]; then
    echo "Error: ONT FASTQ file ${ONT_FASTQ_GZ} not found or not readable. Exiting."
fi

# Validate that PHASED_VCF_GZ ends with .gz
if [[ "${PHASED_VCF_GZ}" != *.gz ]]; then
    echo "Error: Phased VCF file must be gzipped and have a .gz extension."
    echo "Received: ${PHASED_VCF_GZ}"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# --- Step 1: Index Phased VCF ---
echo "Step 1: Indexing phased VCF (${PHASED_VCF_GZ})..."
# The index file (.tbi) will be created in the same directory as PHASED_VCF_GZ.
if [ ! -f "${PHASED_VCF_GZ}.tbi" ]; then
    echo "Index file ${PHASED_VCF_GZ}.tbi not found. Creating index..."
    
    tabix -p vcf "${PHASED_VCF_GZ}"
    
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create index for ${PHASED_VCF_GZ}. Exiting."
        exit 1
    else
        echo "Successfully created index for ${PHASED_VCF_GZ}."
    fi
else
    echo "Index file ${PHASED_VCF_GZ}.tbi already exists. Skipping indexing."
fi
echo "Step 1 finished."

# --- Step 2: Haplotagging Normal BAM (Whole Genome) ---
echo "Step 2: Haplotagging Normal BAM file (Whole Genome)..."

NORMAL_HAPLOTAG_BAM="${OUTPUT_DIR}/normal.re_haplotag.bam"
NORMAL_HAPLOTAG_TSV="${OUTPUT_DIR}/normal.re_haplotag.haplotype.tsv.gz"

echo "Haplotagging Normal BAM: ${NORMAL_BAM}"
echo "Output BAM: ${NORMAL_HAPLOTAG_BAM}"
echo "Output TSV: ${NORMAL_HAPLOTAG_TSV}"

singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${WHATSHAP_CONTAINER_SIF_PATH}" \
  whatshap haplotag \
    "${PHASED_VCF_GZ}" \
    "${NORMAL_BAM}" \
    -r "${REF_FASTA_PATH}" \
    -o "${NORMAL_HAPLOTAG_BAM}" \
    --output-threads "${NSLOTS}" \
    --output-haplotag-list "${NORMAL_HAPLOTAG_TSV}" \
    --ignore-read-groups \
    --tag-supplementary \
    --skip-missing-contigs

if [ $? -ne 0 ]; then
    echo "Error: whatshap haplotag failed for Normal BAM. Exiting."
    exit 1
else
    echo "Successfully haplotagged Normal BAM: ${NORMAL_HAPLOTAG_BAM}"
fi

echo "Step 2 finished."

# --- Step 3: Create BED file from Stats File ---
echo "Step 3: Creating BED file from ${STATS_FILE}..."
REGIONS_BED="${OUTPUT_DIR}/regions_to_extract.bed"

# Assumes STATS_FILE is tab-delimited, has a header, 'diff' is the 9th column.
# Assumes chr_num, start_block_pos, end_block_pos are 1st, 2nd, 3rd columns.
awk 'BEGIN{FS="\t"; OFS="\t"} NR>1 && $9=="True" {print $1, $2, $3}' "${STATS_FILE}" > "${REGIONS_BED}"

if [ $? -ne 0 ]; then
    echo "Error: Failed to create BED file ${REGIONS_BED} from ${STATS_FILE}. Exiting."
    exit 1
fi
if [ ! -s "${REGIONS_BED}" ]; then
    echo "Warning: BED file ${REGIONS_BED} is empty. This may be expected if no rows in ${STATS_FILE} have diff=True."
    echo "Continuing, but subsequent filtering steps might yield no reads."
else
    echo "Successfully created BED file: ${REGIONS_BED}"
fi
echo "Step 3 finished."

# --- Step 4: Filter BAM by HP tag and BED, then convert to FASTQ ---
echo "Step 4: Filtering BAM and converting to FASTQ..."

HP1_FILTERED_BAM="${OUTPUT_DIR}/HP1_filtered.bam"
HP2_FILTERED_BAM="${OUTPUT_DIR}/HP2_filtered.bam"
HP1_FILTERED_FQ_GZ="${OUTPUT_DIR}/HP1_filtered.fq.gz"
HP2_FILTERED_FQ_GZ="${OUTPUT_DIR}/HP2_filtered.fq.gz"

if [ ! -f "${NORMAL_HAPLOTAG_BAM}" ]; then
    echo "Error: Haplotagged BAM ${NORMAL_HAPLOTAG_BAM} not found. Cannot proceed with filtering. Exiting."
    exit 1
fi

# Proceed with filtering only if the BED file has content.
if [ ! -s "${REGIONS_BED}" ]; then
    echo "Info: BED file ${REGIONS_BED} is empty. Skipping samtools view and fastq steps for HP1/HP2."
    # Create empty gzipped FASTQ files if they might be expected by downstream processes
    echo "Creating empty ${HP1_FILTERED_FQ_GZ} and ${HP2_FILTERED_FQ_GZ}."
    gzip -c /dev/null > "${HP1_FILTERED_FQ_GZ}"
    gzip -c /dev/null > "${HP2_FILTERED_FQ_GZ}"
    # Optionally create empty BAMs too if needed, though samtools fastq on empty BAM is fine
    touch "${HP1_FILTERED_BAM}" "${HP2_FILTERED_BAM}"

else
    echo "Filtering HP:1 reads from ${NORMAL_HAPLOTAG_BAM} using ${REGIONS_BED}..."
    "${SAMTOOLS_PATH}" view -@ "${SAMTOOLS_THREADS}" -bh -d HP:1 -L "${REGIONS_BED}" "${NORMAL_HAPLOTAG_BAM}" > "${HP1_FILTERED_BAM}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter HP:1 reads into ${HP1_FILTERED_BAM}. Exiting."
        exit 1
    else
        echo "Successfully created ${HP1_FILTERED_BAM}."
    fi

    echo "Filtering HP:2 reads from ${NORMAL_HAPLOTAG_BAM} using ${REGIONS_BED}..."
    "${SAMTOOLS_PATH}" view -@ "${SAMTOOLS_THREADS}" -bh -d HP:2 -L "${REGIONS_BED}" "${NORMAL_HAPLOTAG_BAM}" > "${HP2_FILTERED_BAM}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter HP:2 reads into ${HP2_FILTERED_BAM}. Exiting."
        exit 1
    else
        echo "Successfully created ${HP2_FILTERED_BAM}."
    fi

    echo "Converting ${HP1_FILTERED_BAM} to FASTQ (${HP1_FILTERED_FQ_GZ})..."
    if [ -s "${HP1_FILTERED_BAM}" ]; then
        "${SAMTOOLS_PATH}" fastq -@ "${SAMTOOLS_THREADS}" "${HP1_FILTERED_BAM}" | gzip > "${HP1_FILTERED_FQ_GZ}"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to convert ${HP1_FILTERED_BAM} to ${HP1_FILTERED_FQ_GZ}. Exiting."
            exit 1
        else
            echo "Successfully created ${HP1_FILTERED_FQ_GZ}."
        fi
    else
        echo "Info: ${HP1_FILTERED_BAM} is empty. Creating empty ${HP1_FILTERED_FQ_GZ}."
        gzip -c /dev/null > "${HP1_FILTERED_FQ_GZ}"
    fi

    echo "Converting ${HP2_FILTERED_BAM} to FASTQ (${HP2_FILTERED_FQ_GZ})..."
    if [ -s "${HP2_FILTERED_BAM}" ]; then
        "${SAMTOOLS_PATH}" fastq -@ "${SAMTOOLS_THREADS}" "${HP2_FILTERED_BAM}" | gzip > "${HP2_FILTERED_FQ_GZ}"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to convert ${HP2_FILTERED_BAM} to ${HP2_FILTERED_FQ_GZ}. Exiting."
            exit 1
        else
            echo "Successfully created ${HP2_FILTERED_FQ_GZ}."
        fi
    else
        echo "Info: ${HP2_FILTERED_BAM} is empty. Creating empty ${HP2_FILTERED_FQ_GZ}."
        gzip -c /dev/null > "${HP2_FILTERED_FQ_GZ}"
    fi
fi
echo "Step 4 finished."

echo "--- Script finished successfully ---"
echo "All output files are located in: ${OUTPUT_DIR}"
echo "Key outputs include:"
echo "  - Haplotagged BAM: ${NORMAL_HAPLOTAG_BAM}"
echo "  - Regions BED: ${REGIONS_BED}"
# Only list filtered files if the BED file was not empty and thus filtering was attempted
if [ -s "${REGIONS_BED}" ]; then
    echo "  - Filtered HP1 BAM: ${HP1_FILTERED_BAM}"
    echo "  - Filtered HP2 BAM: ${HP2_FILTERED_BAM}"
    echo "  - Filtered HP1 FASTQ: ${HP1_FILTERED_FQ_GZ}"
    echo "  - Filtered HP2 FASTQ: ${HP2_FILTERED_FQ_GZ}"
    if [ -f "${OUTPUT_DIR}/hp1.yak" ]; then # Check if yak files were created
        echo "  - HP1 Yak File: ${OUTPUT_DIR}/hp1.yak"
        echo "  - HP2 Yak File: ${OUTPUT_DIR}/hp2.yak"
        echo "  - Hifiasm Assembly Log: ${OUTPUT_DIR}/hifiasm_assembly.log"
        echo "  - Hifiasm Assembly GFA (example): ${OUTPUT_DIR}/hifiasm_assembly.bp.p_ctg.gfa"
    fi
fi

# --- Step 5: Run yak count ---
echo "Step 5: Running yak count..."
HP1_YAK_FILE="${OUTPUT_DIR}/hp1.yak"
HP2_YAK_FILE="${OUTPUT_DIR}/hp2.yak"
YAK_RUN_SUCCESSFUL=false

if [ -s "${HP1_FILTERED_FQ_GZ}" ]; then
    echo "Running yak count for HP1 reads: ${HP1_FILTERED_FQ_GZ} -> ${HP1_YAK_FILE}"
    "${YAK_PATH}" count -b37 -t"${YAK_THREADS}" -o "${HP1_YAK_FILE}" "${HP1_FILTERED_FQ_GZ}"
    if [ $? -ne 0 ]; then
        echo "Error: yak count failed for ${HP1_FILTERED_FQ_GZ}. Exiting."
        exit 1
    else
        echo "Successfully created ${HP1_YAK_FILE}."
    fi
else
    echo "Info: ${HP1_FILTERED_FQ_GZ} is empty. Skipping yak count for HP1."
    # Create an empty placeholder if hifiasm requires the file to exist, or handle in hifiasm step
    touch "${HP1_YAK_FILE}" # hifiasm might fail if file doesn't exist, even if empty
fi

if [ -s "${HP2_FILTERED_FQ_GZ}" ]; then
    echo "Running yak count for HP2 reads: ${HP2_FILTERED_FQ_GZ} -> ${HP2_YAK_FILE}"
    "${YAK_PATH}" count -b37 -t"${YAK_THREADS}" -o "${HP2_YAK_FILE}" "${HP2_FILTERED_FQ_GZ}"
    if [ $? -ne 0 ]; then
        echo "Error: yak count failed for ${HP2_FILTERED_FQ_GZ}. Exiting."
        exit 1
    else
        echo "Successfully created ${HP2_YAK_FILE}."
        YAK_RUN_SUCCESSFUL=true # Set to true only if both could potentially run and HP2 succeeded
    fi
else
    echo "Info: ${HP2_FILTERED_FQ_GZ} is empty. Skipping yak count for HP2."
    touch "${HP2_YAK_FILE}"
fi

# Check if both yak files are non-empty if we expect them for hifiasm
if [ -s "${HP1_YAK_FILE}" ] && [ -s "${HP2_YAK_FILE}" ]; then
    YAK_RUN_SUCCESSFUL=true
elif [ ! -s "${HP1_FILTERED_FQ_GZ}" ] && [ ! -s "${HP2_FILTERED_FQ_GZ}" ]; then
    echo "Info: Both HP1 and HP2 filtered FASTQ files were empty. Yak files will be empty."
    # YAK_RUN_SUCCESSFUL remains false, hifiasm will be skipped.
else
    echo "Warning: One or both yak files are empty or were not generated from non-empty FASTQ. Hifiasm might not run as expected."
    YAK_RUN_SUCCESSFUL=false # Ensure hifiasm is skipped if inputs are problematic
fi
echo "Step 5 finished."


# --- Step 6: Run hifiasm ---
echo "Step 6: Running hifiasm..."
HIFIASM_OUTPUT_PREFIX="${OUTPUT_DIR}/hifiasm_assembly"
HIFIASM_LOG="${OUTPUT_DIR}/hifiasm_assembly.log"

if [ "$YAK_RUN_SUCCESSFUL" = true ]; then
    echo "Running hifiasm with:"
    echo "  ONT FASTQ: ${ONT_FASTQ_GZ}"
    echo "  HP1 Yak: ${HP1_YAK_FILE}"
    echo "  HP2 Yak: ${HP2_YAK_FILE}"
    echo "  Output Prefix: ${HIFIASM_OUTPUT_PREFIX}"
    echo "  Threads: ${HIFIASM_THREADS}"
    echo "  Log: ${HIFIASM_LOG}"
    
    # HIFIASM_THREADS is now aligned with NSLOTS.
    # Added --ont -i options as requested.
    "${HIFIASM_PATH}" --ont -i -l1 --dual-scaf -t"${HIFIASM_THREADS}" \
        -1 "${HP1_YAK_FILE}" -2 "${HP2_YAK_FILE}" \
        -o "${HIFIASM_OUTPUT_PREFIX}" \
        "${ONT_FASTQ_GZ}" > "${HIFIASM_LOG}" 2>&1
    
    if [ $? -ne 0 ]; then
        echo "Error: hifiasm failed. Check log ${HIFIASM_LOG}. Exiting."
        # Optionally, display last few lines of log: tail "${HIFIASM_LOG}"
        exit 1
    else
        echo "Successfully ran hifiasm. Assembly outputs prefixed with ${HIFIASM_OUTPUT_PREFIX}."
        echo "Hifiasm log is at ${HIFIASM_LOG}."
    fi
else
    echo "Info: Skipping hifiasm because yak files were not successfully generated (likely due to empty input FASTQs)."
fi
echo "Step 6 finished."
