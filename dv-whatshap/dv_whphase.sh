#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J dv_whphase
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
# Usage: ./dv_whphase.sh --type <ont|hifi> [--reference hg38|chm13] [--strict-filter] <input.bam>

# Start time
START_TIME=$(date +%s)

# Parse arguments
SEQ_TYPE=""
INPUT_BAM=""
REFERENCE_TYPE="hg38"
STRICT_FILTER=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --type)
            if [ "$2" = "ont" ] || [ "$2" = "hifi" ]; then
                SEQ_TYPE="$2"
            else
                echo "Error: --type must be 'ont' or 'hifi'"
                exit 1
            fi
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
        --strict-filter)
            STRICT_FILTER=true
            shift
            ;;
        *)
            INPUT_BAM="$1"
            shift
            ;;
    esac
done

# Check required arguments
if [ -z "$SEQ_TYPE" ] || [ -z "$INPUT_BAM" ]; then
    echo "Error: --type and input BAM file are required"
    exit 1
fi

# Convert input BAM to absolute path
INPUT_BAM=$(realpath "$INPUT_BAM")

# Thread definition
THREADS=${SLURM_CPUS_PER_TASK:-16}

# Set DeepVariant model type
if [ "$SEQ_TYPE" = "ont" ]; then
    DV_MODEL="ONT_R104"
elif [ "$SEQ_TYPE" = "hifi" ]; then
    DV_MODEL="PACBIO"
fi

# Set reference genome path
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

# Container paths
DEEPVARIANT_SIF="/home/itoyu8/singularity/deepvariant_1.9.0.sif"
WHATSHAP_SIF="/home/itoyu8/singularity/scarpia-python_0.2.0.sif"

# BCFtools path
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"

# [ARCHIVED] gnomAD VCF for filtering (currently unused):
# GNOMAD_VCF="/home/itoyu8/database/reference/gnomAD_4.1/gnomad.genomes.v4.1.sites.merged.light.vcf.bgz"

# Get directory from input BAM (now absolute path)
OUTPUT_DIR=$(dirname "$INPUT_BAM")

# Output files in same directory as input BAM
DV_OUTPUT="${OUTPUT_DIR}/dv.vcf.gz"
FILTERED_OUTPUT="${OUTPUT_DIR}/dv.filtered.vcf.gz"
PHASED_OUTPUT="${OUTPUT_DIR}/phased.vcf.gz"

# DeepVariant intermediate directory (avoid /tmp capacity issues in Singularity)
DV_TEMP_DIR="${OUTPUT_DIR}/dv_intermediate"
mkdir -p "${DV_TEMP_DIR}"

mkdir -p log

# Step 1: Run DeepVariant
echo "Running DeepVariant with model: $DV_MODEL"
DV_START=$(date +%s)

singularity run --nv \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1/ \
    "${DEEPVARIANT_SIF}" run_deepvariant \
    --model_type "${DV_MODEL}" \
    --ref "${REFERENCE_GENOME_PATH}" \
    --reads "${INPUT_BAM}" \
    --output_vcf "${DV_OUTPUT}" \
    --intermediate_results_dir "${DV_TEMP_DIR}" \
    --num_shards "${THREADS}"

DV_END=$(date +%s)
DV_ELAPSED=$((DV_END - DV_START))
echo "DeepVariant completed: ${DV_OUTPUT}"
echo "DeepVariant time: ${DV_ELAPSED} seconds"

# Step 2: Pass DV output to WhatsHap (with optional strict filter)
#
# [ARCHIVED] gnomAD filtering commands - can be restored if needed:
# ----------------------------------------------------------------
# # Filter for known SNPs only (hg38):
# "${BCFTOOLS}" isec -n=2 -w2 -c snps -O z -o "${UNPHASED_OUTPUT}" "${GNOMAD_VCF}" "${DV_OUTPUT}"
#
# # Filter for known variants including indels (hg38):
# "${BCFTOOLS}" isec -n=2 -w2 -O z -o "${UNPHASED_OUTPUT}" "${GNOMAD_VCF}" "${DV_OUTPUT}"
#
# # gnomAD + GQ + VAF filter (hg38):
# GNOMAD_TEMP="${OUTPUT_DIR}/temp_gnomad.vcf.gz"
# "${BCFTOOLS}" isec -n=2 -w2 -O z -o "${GNOMAD_TEMP}" "${GNOMAD_VCF}" "${DV_OUTPUT}"
# tabix -p vcf "${GNOMAD_TEMP}"
# "${BCFTOOLS}" view -f PASS -m2 -M2 --genotype het "${GNOMAD_TEMP}" | \
# "${BCFTOOLS}" filter -i "FORMAT/GQ >= 20 && FORMAT/VAF >= 0.3 && FORMAT/VAF <= 0.7" -O z -o "${FILTERED_OUTPUT}"
# rm -f "${GNOMAD_TEMP}" "${GNOMAD_TEMP}.tbi"
# ----------------------------------------------------------------

STRICT_FILTER_ELAPSED=0
WHATSHAP_INPUT="${DV_OUTPUT}"

if [ "$STRICT_FILTER" = true ]; then
    echo "Applying strict filter (GQ + VAF)..."
    STRICT_FILTER_START=$(date +%s)

    MIN_GQ=20
    MIN_VAF=0.3
    MAX_VAF=0.7

    # GQ + VAF filter for heterozygous biallelic variants
    "${BCFTOOLS}" view \
        -f PASS \
        -m2 -M2 \
        --genotype het \
        "${DV_OUTPUT}" | \
    "${BCFTOOLS}" filter \
        -i "FORMAT/GQ >= ${MIN_GQ} && FORMAT/VAF >= ${MIN_VAF} && FORMAT/VAF <= ${MAX_VAF}" \
        -O z -o "${FILTERED_OUTPUT}"

    tabix -p vcf "${FILTERED_OUTPUT}"
    WHATSHAP_INPUT="${FILTERED_OUTPUT}"

    STRICT_FILTER_END=$(date +%s)
    STRICT_FILTER_ELAPSED=$((STRICT_FILTER_END - STRICT_FILTER_START))
    echo "Strict filter completed: ${FILTERED_OUTPUT}"
    echo "Strict filter time: ${STRICT_FILTER_ELAPSED} seconds"
else
    echo "No filtering applied - passing DeepVariant output directly to WhatsHap"
fi

# Step 3: Run Whatshap
echo "Running Whatshap for phasing..."
WHATSHAP_START=$(date +%s)

singularity exec --nv \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1/ \
    "${WHATSHAP_SIF}" whatshap phase \
    --reference "${REFERENCE_GENOME_PATH}" \
    --ignore-read-groups \
    --distrust-genotypes \
    -o "${PHASED_OUTPUT}" \
    "${WHATSHAP_INPUT}" \
    "${INPUT_BAM}"

WHATSHAP_END=$(date +%s)
WHATSHAP_ELAPSED=$((WHATSHAP_END - WHATSHAP_START))
echo "Whatshap completed: ${PHASED_OUTPUT}"
echo "Whatshap time: ${WHATSHAP_ELAPSED} seconds"

# Step 3: Index phased VCF with tabix
echo "Indexing phased VCF with tabix..."
tabix -p vcf "${PHASED_OUTPUT}"
echo "Indexing completed"

# End time
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=== Processing Summary ==="
echo "DeepVariant: ${DV_ELAPSED} seconds"
if [ "$STRICT_FILTER" = true ]; then
    echo "Strict filter: ${STRICT_FILTER_ELAPSED} seconds"
fi
echo "WhatsHap: ${WHATSHAP_ELAPSED} seconds"
echo "Total time: ${TOTAL_ELAPSED} seconds"
echo ""
echo "Pipeline completed successfully"
echo "DeepVariant VCF: ${DV_OUTPUT}"
if [ "$STRICT_FILTER" = true ]; then
    echo "Filtered VCF: ${FILTERED_OUTPUT}"
fi
echo "Phased VCF: ${PHASED_OUTPUT}"
