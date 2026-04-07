#!/bin/bash
### QUILT2_prepare_reference を SHAPEIT4 パネルで全染色体・全チャンクに対して実行するスクリプト

#SBATCH -p rjobs,mjobs
#SBATCH -J quilt_prep_ref_shapeit4
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=32G
#SBATCH -c 1
# Usage: sbatch prepare_reference_shapeit4.sh

set -euxo pipefail

QUILT_CONTAINER_SIF_PATH="/home/itoyu8/singularity/quilt_v0.1.0.sif"

CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")

GENETIC_MAP_DIR="/home/itoyu8/database/tools/quilt/hg38/maps"
VCF_DIR="/home/itoyu8/database/1000genomes/hg38/shapeit4_phased"
CHUNK_DIR="/home/itoyu8/database/tools/quilt/hg38/chunk_output"
OUTPUT_DIR="/home/itoyu8/database/tools/quilt/hg38/prepared_reference_shapeit4"
NGEN=100
BUFFER=500000

mkdir -p $OUTPUT_DIR
mkdir -p log

echo "Starting QUILT2_prepare_reference (SHAPEIT4) for all chromosomes and chunks..."
echo "Output directory: $OUTPUT_DIR"
echo "=========================="

for chr_num in "${CHROMOSOMES[@]}"; do
    chr="chr${chr_num}"
    chunk_file="${CHUNK_DIR}/chunks_${chr}.txt"
    genetic_map_file="${GENETIC_MAP_DIR}/CEU-${chr}-final.b38.txt.gz"
    vcf_file="${VCF_DIR}/1kGP_high_coverage_Illumina.${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

    echo "Processing $chr..."

    if [ ! -f "$chunk_file" ]; then
        echo "Warning: Chunk file not found: $chunk_file"
        continue
    fi

    if [ ! -f "$genetic_map_file" ]; then
        echo "Warning: Genetic map file not found: $genetic_map_file"
        continue
    fi

    if [ ! -f "$vcf_file" ]; then
        echo "Warning: VCF file not found: $vcf_file"
        continue
    fi

    tail -n +2 "$chunk_file" | while IFS=$'\t' read -r chunk_id chr_name region; do
        echo "  Processing chunk $chunk_id: $region"

        region_start=$(echo $region | sed 's/.*://' | sed 's/-.*//')
        region_end=$(echo $region | sed 's/.*-//')

        singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ "${QUILT_CONTAINER_SIF_PATH}" \
            /bin/QUILT2_prepare_reference.R \
            --genetic_map_file="$genetic_map_file" \
            --reference_vcf_file="$vcf_file" \
            --chr="$chr" \
            --regionStart="$region_start" \
            --regionEnd="$region_end" \
            --nGen="$NGEN" \
            --buffer="$BUFFER" \
            --outputdir="$OUTPUT_DIR"

        if [ $? -eq 0 ]; then
            echo "    Successfully processed chunk $chunk_id ($region)"
        else
            echo "    Error processing chunk $chunk_id ($region)"
        fi
    done

    echo "Completed $chr"
    echo "------------------------"
done

echo "=========================="
echo "QUILT2_prepare_reference (SHAPEIT4) completed."
echo "Output saved to: $OUTPUT_DIR"
echo "Exit status: $?"
