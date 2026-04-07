#!/bin/bash
### QUILT2_prepare_reference を指定された染色体の全チャンクに対して実行する
### 事前に split_bcf.sh で per_chr_vcf/ を作成しておくこと

#SBATCH -p rjobs,mjobs
#SBATCH -J quilt_prep_ref_chm13
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=32G
#SBATCH -c 1
# Usage: sbatch prepare_reference.sh <chr1|chr2|...|chr22>
# All chromosomes: for i in {1..22}; do sbatch prepare_reference.sh chr$i; done

set -euxo pipefail

CHR=$1

QUILT_CONTAINER_SIF_PATH="/home/itoyu8/singularity/quilt_v0.1.0.sif"
VCF_DIR="/home/itoyu8/database/tools/quilt/chm13/per_chr_vcf"
GENETIC_MAP_DIR="/home/itoyu8/database/tools/quilt/chm13/maps"
CHUNK_DIR="/home/itoyu8/database/tools/quilt/chm13/chunk_output"
OUTPUT_DIR="/home/itoyu8/database/tools/quilt/chm13/prepared_reference"
NGEN=100
BUFFER=500000

mkdir -p "${OUTPUT_DIR}"
mkdir -p log

chunk_file="${CHUNK_DIR}/chunks_${CHR}.txt"
genetic_map_file="${GENETIC_MAP_DIR}/CEU-${CHR}-final.chm13.txt.gz"
vcf_file="${VCF_DIR}/${CHR}.vcf.gz"

time {
    tail -n +2 "$chunk_file" | while IFS=$'\t' read -r chunk_id chr_name region; do
        region_start=$(echo $region | sed 's/.*://' | sed 's/-.*//')
        region_end=$(echo $region | sed 's/.*-//')

        singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ "${QUILT_CONTAINER_SIF_PATH}" \
            /bin/QUILT2_prepare_reference.R \
            --genetic_map_file="$genetic_map_file" \
            --reference_vcf_file="$vcf_file" \
            --chr="$CHR" \
            --regionStart="$region_start" \
            --regionEnd="$region_end" \
            --nGen="$NGEN" \
            --buffer="$BUFFER" \
            --outputdir="$OUTPUT_DIR"
    done
}

echo "Exit status: $?"
