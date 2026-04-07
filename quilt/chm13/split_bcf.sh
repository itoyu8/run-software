#!/bin/bash
### 全ゲノム BCF を染色体ごとの VCF.gz に分割する (1回だけ実行)

#SBATCH -p rjobs,mjobs
#SBATCH -J split_bcf_chm13
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
# Usage: sbatch split_bcf.sh

set -euxo pipefail

BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"
INPUT_BCF="/home/itoyu8/database/1000genomes/chm13/whole_genome/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf.gz"
OUTPUT_DIR="/home/itoyu8/database/tools/quilt/chm13/per_chr_vcf"

mkdir -p "${OUTPUT_DIR}"
mkdir -p log

time {
    for chr_num in {1..22}; do
        chr="chr${chr_num}"
        "${BCFTOOLS}" view -r "${chr}" "${INPUT_BCF}" -Oz -o "${OUTPUT_DIR}/${chr}.vcf.gz"
        "${BCFTOOLS}" index -t "${OUTPUT_DIR}/${chr}.vcf.gz"
    done
}

echo "Exit status: $?"
