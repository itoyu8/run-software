#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J ncgm_pure
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 4

CHR=$1

if [ -z "$CHR" ]; then
    echo "Usage: sbatch $0 <chr>"
    exit 1
fi

INPUT_DIR="/home/itoyu8/database/1000genomes/ncbn_window40_Map"
OUTPUT_DIR="/home/itoyu8/database/1000genomes/ncbn_pure"
SAMPLE_LIST="/home/kechiba/work_ncbn_211130/database/1000GenomesReferencePanelSamples.txt"
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"

mkdir -p "$OUTPUT_DIR"

INPUT_VCF="${INPUT_DIR}/NCBN_genotype.beagle_${CHR}.vcf.gz"
OUTPUT_VCF="${OUTPUT_DIR}/NCBN_genotype.beagle_${CHR}.vcf.gz"
"$BCFTOOLS" view -S ^"${SAMPLE_LIST}" --force-samples --min-ac 1 --threads 4 -O z -o "${OUTPUT_VCF}" "${INPUT_VCF}"
tabix -p vcf "${OUTPUT_VCF}"
