#!/bin/bash
### genetic mapのファイルから染色体ごとにchunkを作ってtxtファイルで書き出すスクリプト

#SBATCH -J quilt_chunk_mapping
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 1

# Set the path to your QUILT Singularity container
QUILT_CONTAINER_SIF_PATH="/home/itoyu8/singularity/quilt_v0.1.0.sif"

# Set genetic map directory
GENETIC_MAP_DIR="/home/itoyu8/database/tools/quilt/maps/hg38"

# Set output directory
OUTDIR="chunk_output"
mkdir -p $OUTDIR

# Create R script to generate chunks for all chromosomes
cat << 'EOF' > temp_chunk_script.R
library(QUILT)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
genetic_map_dir <- args[1]
output_dir <- args[2]

# Process all chromosomes
for (chr_num in 1:22) {
  chr <- paste0("chr", chr_num)
  genetic_map_file <- file.path(genetic_map_dir, paste0("CEU-", chr, "-final.b38.txt.gz"))
  output_file <- file.path(output_dir, paste0("chunks_", chr, ".txt"))
  
  if (file.exists(genetic_map_file)) {
    cat("Processing", chr, "...\n")
    
    # Generate chunk mapping
    dat <- QUILT::quilt_chunk_map(chr, genetic_map_file)
    
    # Write to file
    write.table(dat, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("Saved", nrow(dat), "chunks for", chr, "to", output_file, "\n")
  } else {
    cat("Warning: Genetic map file not found for", chr, ":", genetic_map_file, "\n")
  }
}

cat("Chunk mapping completed for all chromosomes.\n")
EOF

# Run the R script inside the container
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "${QUILT_CONTAINER_SIF_PATH}" \
  Rscript temp_chunk_script.R "${GENETIC_MAP_DIR}" "${OUTDIR}"

# Clean up temporary script
rm temp_chunk_script.R

echo "Chunk mapping completed. Results saved to: ${OUTDIR}/"
echo "List of chunk files:"
ls -la ${OUTDIR}/chunks_*.txt