#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J qc_dist
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=2G
#SBATCH -c 8
# Usage: sbatch qc_dist.sh <input.fastq.gz> [output.tsv]

INPUT_FILE="$1"
OUTPUT_FILE="${2:-$(dirname "$INPUT_FILE")/$(basename "$INPUT_FILE" .fastq.gz).qc_dist.tsv}"

cat "$INPUT_FILE" | singularity run /home/itoyu8/singularity/fastq-qc-counter_0.2.0.sif > "$OUTPUT_FILE"
