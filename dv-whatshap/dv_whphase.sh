#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J dv_whphase
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
# Usage: ./dv_whatshap.sh --type <ont|hifi> [--reference hg38|chm13] <input.bam>

# Start time
START_TIME=$(date +%s)

# Parse arguments
SEQ_TYPE=""
INPUT_BAM=""
REFERENCE_TYPE="hg38"

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

# BCFtools and gnomAD (only for hg38)
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"
GNOMAD_VCF="/home/itoyu8/database/reference/gnomAD_4.1/gnomad.genomes.v4.1.sites.merged.light.vcf.bgz"

# Get base name and directory from input BAM
BASENAME=$(basename "$INPUT_BAM" .bam)
OUTPUT_DIR=$(dirname "$INPUT_BAM")

# Output files in same directory as input BAM
DV_OUTPUT="${OUTPUT_DIR}/normal_dv.vcf.gz"
UNPHASED_OUTPUT="${OUTPUT_DIR}/snp.unphased.vcf.gz"
PHASED_OUTPUT="${OUTPUT_DIR}/snp.phased.vcf.gz"

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
    --num_shards "${THREADS}"

DV_END=$(date +%s)
DV_ELAPSED=$((DV_END - DV_START))
echo "DeepVariant completed: ${DV_OUTPUT}"
echo "DeepVariant time: ${DV_ELAPSED} seconds"

# Step 2: Filter with bcftools (only for hg38)
BCFTOOLS_START=$(date +%s)

if [ "$REFERENCE_TYPE" = "hg38" ]; then
    echo "Filtering VCF with bcftools for known SNPs..."
    "${BCFTOOLS}" isec -n=2 -w2 -c snps -O z -o "${UNPHASED_OUTPUT}" "${GNOMAD_VCF}" "${DV_OUTPUT}"
    echo "Filtered VCF: ${UNPHASED_OUTPUT}"
else
    echo "Skipping bcftools filtering for CHM13"
    mv "${DV_OUTPUT}" "${UNPHASED_OUTPUT}"
fi

BCFTOOLS_END=$(date +%s)
BCFTOOLS_ELAPSED=$((BCFTOOLS_END - BCFTOOLS_START))
echo "BCFtools filtering time: ${BCFTOOLS_ELAPSED} seconds"

# Total time for DeepVariant + bcftools
DV_BCFTOOLS_TOTAL=$((DV_ELAPSED + BCFTOOLS_ELAPSED))
echo "DeepVariant + BCFtools total: ${DV_BCFTOOLS_TOTAL} seconds"

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
    "${UNPHASED_OUTPUT}" \
    "${INPUT_BAM}"

WHATSHAP_END=$(date +%s)
WHATSHAP_ELAPSED=$((WHATSHAP_END - WHATSHAP_START))
echo "Whatshap completed: ${PHASED_OUTPUT}"
echo "Whatshap time: ${WHATSHAP_ELAPSED} seconds"

# Step 4: Index with tabix
echo "Indexing VCF files with tabix..."
tabix -p vcf "${UNPHASED_OUTPUT}"
tabix -p vcf "${PHASED_OUTPUT}"

echo "Indexing completed"

# Cleanup intermediate file if bcftools was used
if [ "$REFERENCE_TYPE" = "hg38" ]; then
    rm -f "${DV_OUTPUT}" "${DV_OUTPUT}.tbi"
fi

# End time
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=== Processing Summary ==="
echo "DeepVariant + BCFtools: ${DV_BCFTOOLS_TOTAL} seconds"
echo "Whatshap: ${WHATSHAP_ELAPSED} seconds"
echo "Total time: ${TOTAL_ELAPSED} seconds"
echo ""
echo "Pipeline completed successfully"
echo "Unphased VCF: ${UNPHASED_OUTPUT}"
echo "Phased VCF: ${PHASED_OUTPUT}"
