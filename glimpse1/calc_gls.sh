#!/bin/bash
#SBATCH -J calc_gls
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=32G
#SBATCH -c 1

# Usage: sbatch calc_gls.sh /path/to/sample.bam [output_folder_name]

BAM=$1
OUTPUT_FOLDER_NAME=${2:-"glimpse1_gl"}

# Set up directories
BAM_DIR=$(dirname "$BAM")
OUTPUT_BASE="${BAM_DIR}/${OUTPUT_FOLDER_NAME}"
REFPANEL_DIR="/home/itoyu8/database/tools/glimpse1/reference_panel"
REFGEN="/home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta"
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"

mkdir -p "${OUTPUT_BASE}"
mkdir -p log

# Process chromosomes 1-22
for CHR in {1..22}; do
    chr="chr${CHR}"
    
    # Compute genotype likelihoods using pre-created site files
    VCF="${REFPANEL_DIR}/1000GP.chr${CHR}.sites.vcf.gz"
    TSV="${REFPANEL_DIR}/1000GP.chr${CHR}.sites.tsv.gz"
    OUT="${OUTPUT_BASE}/sample.chr${CHR}.vcf.gz"
    
    "${BCFTOOLS}" mpileup -f "${REFGEN}" -I -E -a 'FORMAT/DP' -T "${VCF}" -r "${chr}" "${BAM}" -Ou |
    "${BCFTOOLS}" call -Aim -C alleles -T "${TSV}" -Oz -o "${OUT}"
    
    "${BCFTOOLS}" index -f "${OUT}"
done

# Process chromosome X
VCF="${REFPANEL_DIR}/1000GP.chrX.sites.vcf.gz"
TSV="${REFPANEL_DIR}/1000GP.chrX.sites.tsv.gz"
OUT="${OUTPUT_BASE}/sample.chrX.vcf.gz"

"${BCFTOOLS}" mpileup -f "${REFGEN}" -I -E -a 'FORMAT/DP' -T "${VCF}" -r "chrX" "${BAM}" -Ou |
"${BCFTOOLS}" call -Aim -C alleles -T "${TSV}" -Oz -o "${OUT}"

"${BCFTOOLS}" index -f "${OUT}"