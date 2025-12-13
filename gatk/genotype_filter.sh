#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J gatk_genotype
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 4
# Usage: sbatch genotype_filter.sh [--reference hg38|chm13] -V <input1.g.vcf.gz> [-V <input2.g.vcf.gz> ...] -o <output_dir>

set -e

# Parse arguments
REFERENCE_TYPE="hg38"
INPUT_GVCFS=()
OUTPUT_DIR=""

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
        -V)
            INPUT_GVCFS+=("$2")
            shift 2
            ;;
        -o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown option $1"
            echo "Usage: sbatch $0 [--reference hg38|chm13] -V <input1.g.vcf.gz> [-V <input2.g.vcf.gz> ...] -o <output_dir>"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ ${#INPUT_GVCFS[@]} -lt 1 ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: sbatch $0 [--reference hg38|chm13] -V <input1.g.vcf.gz> [-V <input2.g.vcf.gz> ...] -o <output_dir>"
    exit 1
fi

# Check input GVCF format
for GVCF in "${INPUT_GVCFS[@]}"; do
    if [[ ! "${GVCF}" == *.g.vcf.gz ]]; then
        echo "Error: Input GVCF file '${GVCF}' does not end with .g.vcf.gz"
        exit 1
    fi
done

# Set reference-specific paths
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REF_FASTA="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REF_FASTA="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

# Tool paths
GATK_CONTAINER="/home/itoyu8/singularity/compat_parabricks-0.2.2.sif"
GATK_JAR="/tools/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar"
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"

# Create output directory
OUTPUT_BASE=$(basename "$OUTPUT_DIR")
mkdir -p "$OUTPUT_DIR"
mkdir -p log

# Define output files
RAW_VCF="${OUTPUT_DIR}/${OUTPUT_BASE}.vcf.gz"
SNP_VCF="${OUTPUT_DIR}/${OUTPUT_BASE}.snp.vcf.gz"
FILT_VCF="${OUTPUT_DIR}/${OUTPUT_BASE}.snp.filter.vcf.gz"

# Check and index input GVCFs if needed
for GVCF in "${INPUT_GVCFS[@]}"; do
    GVCF_INDEX="${GVCF}.tbi"
    if [ ! -f "${GVCF_INDEX}" ]; then
        singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${GATK_CONTAINER}" tabix -p vcf "${GVCF}"
    fi
done

# Build -V arguments for GenotypeGVCFs
GVCF_ARGS=""
for GVCF in "${INPUT_GVCFS[@]}"; do
    GVCF_ARGS="${GVCF_ARGS} -V \"${GVCF}\""
done

# Step 1: GenotypeGVCFs
eval singularity exec --bind /home/itoyu8/:/home/itoyu8/ \"${GATK_CONTAINER}\" /usr/bin/java \
    -Xmx6G -jar \"${GATK_JAR}\" GenotypeGVCFs \
    -R \"${REF_FASTA}\" \
    ${GVCF_ARGS} \
    -O \"${RAW_VCF}\"

# Step 2: SelectVariants (SNP only)
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${GATK_CONTAINER}" /usr/bin/java \
    -Xmx6G -jar "${GATK_JAR}" SelectVariants \
    -V "${RAW_VCF}" \
    -select-type SNP \
    -O "${SNP_VCF}"

# Step 3: VariantFiltration
TEMP_FILT_VCF="${OUTPUT_DIR}/temp.filt.vcf.gz"
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${GATK_CONTAINER}" /usr/bin/java \
    -Xmx6G -jar "${GATK_JAR}" VariantFiltration \
    -V "${SNP_VCF}" \
    -O "${TEMP_FILT_VCF}" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

# Step 4: Filter for PASS variants only
"${BCFTOOLS}" view -f PASS -O z "${TEMP_FILT_VCF}" > "${FILT_VCF}"

# Clean up temporary file
rm -f "${TEMP_FILT_VCF}" "${TEMP_FILT_VCF}.tbi"

# Index final VCF
tabix -p vcf "${FILT_VCF}"
