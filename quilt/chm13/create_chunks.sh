#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J quilt_chunk_chm13
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 1
# Usage: sbatch create_chunks.sh

set -euxo pipefail

QUILT_CONTAINER_SIF_PATH="/home/itoyu8/singularity/quilt_v0.1.0.sif"
GENETIC_MAP_DIR="/home/itoyu8/database/tools/quilt/chm13/maps"
OUTDIR="/home/itoyu8/database/tools/quilt/chm13/chunk_output"

mkdir -p "${OUTDIR}"
mkdir -p log

cat << 'EOF' > temp_chunk_script.R
library(QUILT)

args <- commandArgs(trailingOnly = TRUE)
genetic_map_dir <- args[1]
output_dir <- args[2]

for (chr_num in 1:22) {
  chr <- paste0("chr", chr_num)
  genetic_map_file <- file.path(genetic_map_dir, paste0("CEU-", chr, "-final.chm13.txt.gz"))
  output_file <- file.path(output_dir, paste0("chunks_", chr, ".txt"))

  if (file.exists(genetic_map_file)) {
    cat("Processing", chr, "...\n")
    dat <- QUILT::quilt_chunk_map(chr, genetic_map_file)
    write.table(dat, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Saved", nrow(dat), "chunks for", chr, "to", output_file, "\n")
  } else {
    cat("Warning: Genetic map file not found for", chr, ":", genetic_map_file, "\n")
  }
}

cat("Chunk mapping completed for all chromosomes.\n")
EOF

singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ "${QUILT_CONTAINER_SIF_PATH}" \
  Rscript temp_chunk_script.R "${GENETIC_MAP_DIR}" "${OUTDIR}"

rm temp_chunk_script.R

echo "Chunk mapping completed. Results saved to: ${OUTDIR}/"
ls -la ${OUTDIR}/chunks_*.txt

echo "Exit status: $?"
