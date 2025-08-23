#!/bin/bash
### QUILT2_prepare_referenceを全染色体・全チャンクに対して実行するスクリプト

#SBATCH -J quilt_prepare_reference
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=32G
#SBATCH -c 1

# Set the path to your QUILT Singularity container
QUILT_CONTAINER_SIF_PATH="/home/itoyu8/singularity/quilt_v0.1.0.sif"

# Define chromosomes to process
CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")

# Set directories and parameters
GENETIC_MAP_DIR="/home/itoyu8/database/tools/quilt/maps/hg38"
VCF_DIR="/home/itoyu8/database/1000genomes/vcf"
CHUNK_DIR="chunk_output"
OUTPUT_DIR="/home/itoyu8/database/tools/quilt/output"
NGEN=100
BUFFER=500000

# Create output directory
mkdir -p $OUTPUT_DIR

# Create log directory
mkdir -p log

echo "Starting QUILT2_prepare_reference for all chromosomes and chunks..."
echo "Output directory: $OUTPUT_DIR"
echo "=========================="

# Process each chromosome
for chr_num in "${CHROMOSOMES[@]}"; do
    chr="chr${chr_num}"
    chunk_file="${CHUNK_DIR}/chunks_${chr}.txt"
    genetic_map_file="${GENETIC_MAP_DIR}/CEU-${chr}-final.b38.txt.gz"
    vcf_file="${VCF_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
    
    echo "Processing $chr..."
    
    # Check if required files exist
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
    
    # Read chunks file and process each chunk (skip header)
    tail -n +2 "$chunk_file" | while IFS=$'\t' read -r chunk_id chr_name region; do
        echo "  Processing chunk $chunk_id: $region"
        
        # Extract regionStart and regionEnd from region (format: chr:start-end)
        region_start=$(echo $region | sed 's/.*://' | sed 's/-.*//')
        region_end=$(echo $region | sed 's/.*-//')
        
        # Run QUILT2_prepare_reference
        singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${QUILT_CONTAINER_SIF_PATH}" \
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
echo "QUILT2_prepare_reference completed for all chromosomes and chunks."
echo "Output saved to: $OUTPUT_DIR"