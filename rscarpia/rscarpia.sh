#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J rscarpia_pipeline
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
# Usage: ./rscarpia.sh <vcf_file> <tumor_fastq> <normal_fastq>

# Start time
START_TIME=$(date +%s)

echo "=== RScarpia Pipeline started at: $(date) ==="

# Parse arguments
if [ "$#" -ne 3 ]; then
    echo "Error: Requires 3 arguments"
    echo "Usage: ./rscarpia.sh <vcf_file> <tumor_fastq> <normal_fastq>"
    exit 1
fi

VCF_FILE="$1"
TUMOR_FASTQ="$2"
NORMAL_FASTQ="$3"

# Check if files exist
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: VCF file not found: $VCF_FILE"
    exit 1
fi

if [ ! -f "$TUMOR_FASTQ" ]; then
    echo "Error: Tumor FASTQ file not found: $TUMOR_FASTQ"
    exit 1
fi

if [ ! -f "$NORMAL_FASTQ" ]; then
    echo "Error: Normal FASTQ file not found: $NORMAL_FASTQ"
    exit 1
fi

# Thread definition
THREADS=${SLURM_CPUS_PER_TASK:-16}

# Container and reference paths
CONTAINER="/home/itoyu8/singularity/rscarpia_0.1.0.sif"
REFERENCE_FASTA="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"

echo "VCF file: $VCF_FILE"
echo "Tumor FASTQ: $TUMOR_FASTQ"
echo "Normal FASTQ: $NORMAL_FASTQ"
echo "Reference: $REFERENCE_FASTA"
echo "Threads: $THREADS"

# Create log directory
mkdir -p log

# Run pipeline binary from container
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "$CONTAINER" \
  /workspace/rscarpia/target/release/pipeline \
  --reference "$REFERENCE_FASTA" \
  --vcf "$VCF_FILE" \
  --tumor "$TUMOR_FASTQ" \
  --normal "$NORMAL_FASTQ"

# End time
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=== RScarpia Pipeline completed at: $(date) ==="
echo "Total time: ${TOTAL_ELAPSED} seconds"
