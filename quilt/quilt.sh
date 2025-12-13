#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J quilt_diploid
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=64G
#SBATCH -c 1

# Usage: sbatch run_quilt.sh /path/to/sample.bam [output_folder_name]
# Note: output_folder_name can include subdirectories (e.g., "results/quilt_analysis")

CONTAINER="/home/itoyu8/singularity/quilt_v0.1.0.sif"
BAM=$1
OUTPUT_FOLDER_NAME=${2:-"quilt_output"}

# Set up directories
BAM_DIR=$(dirname "$BAM")
OUTPUT_BASE="${BAM_DIR}/${OUTPUT_FOLDER_NAME}"
PREPARED_REFERENCE_DIR="/home/itoyu8/database/tools/quilt/output/RData"
CHUNK_DIR="/home/itoyu8/database/tools/quilt/chunk_output"
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.22/bcftools"

# Parameters
NGEN=100
BUFFER=500000

# Create output directories
mkdir -p "${OUTPUT_BASE}/impute"
mkdir -p "${OUTPUT_BASE}/ligate"
mkdir -p log

# Create bamlist and sample names files for QUILT
sample_name=$(basename "$BAM" .bam | sed 's/\..*$//')
bamlist_file="${OUTPUT_BASE}/bamlist.txt"
samplenames_file="${OUTPUT_BASE}/sample_names.txt"

echo "$(realpath "$BAM")" > "$bamlist_file"
echo "$sample_name" > "$samplenames_file"

# Process chromosomes 1-22
for chr_num in {1..22}; do
    chr="chr${chr_num}"
    chunk_file="${CHUNK_DIR}/chunks_${chr}.txt"
    
    [ ! -f "$chunk_file" ] && continue
    
    # Process each chunk (skip header)
    tail -n +2 "$chunk_file" | while IFS=$'\t' read -r chunk_id chr_name region; do
        # Extract regionStart and regionEnd from region (format: chr:start-end)
        region_start=$(echo $region | sed 's/.*://' | sed 's/-.*//')
        region_end=$(echo $region | sed 's/.*-//')
        
        # Construct prepared reference filename
        prepared_reference="${PREPARED_REFERENCE_DIR}/QUILT_prepared_reference.${chr}.${region_start}.${region_end}.RData"
        
        [ ! -f "$prepared_reference" ] && continue
        
        # Output filename
        output_vcf="${OUTPUT_BASE}/impute/quilt2.diploid.${chr}.chunk${chunk_id}.vcf.gz"
        
        # Run QUILT2 diploid imputation
        singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ "${CONTAINER}" \
            /bin/QUILT2.R \
            --prepared_reference_filename="$prepared_reference" \
            --bamlist="$(realpath "$bamlist_file")" \
            --sampleNames_file="$(realpath "$samplenames_file")" \
            --method=diploid \
            --chr="$chr" \
            --regionStart="$region_start" \
            --regionEnd="$region_end" \
            --nGen="$NGEN" \
            --buffer="$BUFFER" \
            --output_filename="$output_vcf"
    done
    
    # Create list of VCF files for this chromosome
    vcf_list="${OUTPUT_BASE}/ligate/list.${chr}.txt"
    ls -1v "${OUTPUT_BASE}/impute/quilt2.diploid.${chr}.chunk"*.vcf.gz > "$vcf_list" 2>/dev/null
    
    [ ! -s "$vcf_list" ] && continue
    
    # Ligate chunks for this chromosome
    ligated_vcf="${OUTPUT_BASE}/ligate/quilt2.diploid.${chr}.ligated.vcf.gz"
    "$BCFTOOLS" concat --ligate --output-type z --output "$ligated_vcf" --file-list "$vcf_list"
done

# Concatenate all chromosomes
all_vcf_list="${OUTPUT_BASE}/ligate/all_chromosomes.txt"
ls -1v "${OUTPUT_BASE}/ligate/quilt2.diploid.chr"*.ligated.vcf.gz > "$all_vcf_list" 2>/dev/null

final_vcf="${OUTPUT_BASE}/quilt.phased.vcf.gz"
[ -s "$all_vcf_list" ] && "$BCFTOOLS" concat --output-type z --output "$final_vcf" --file-list "$all_vcf_list"
[ -f "$final_vcf" ] && "$BCFTOOLS" index -f "$final_vcf"