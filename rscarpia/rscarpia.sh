#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J rscarpia_pipeline
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
# Usage: ./rscarpia.sh [--reference hg38|chm13] <vcf_file> <tumor_fastq> <normal_fastq>

# Start time
START_TIME=$(date +%s)

echo "=== RScarpia Pipeline started at: $(date) ==="

# Parse arguments
INPUT_FILES=()
REFERENCE_TYPE="hg38"

while [[ $# -gt 0 ]]; do
    case $1 in
        --reference)
            if [ "$2" = "chm13" ]; then
                REFERENCE_TYPE="chm13"
            elif [ "$2" = "hg38" ]; then
                REFERENCE_TYPE="hg38"
            else
                echo "Error: --reference must be 'hg38' or 'chm13'"
                exit 1
            fi
            shift 2
            ;;
        *)
            INPUT_FILES+=("$1")
            shift
            ;;
    esac
done

if [ "${#INPUT_FILES[@]}" -ne 3 ]; then
    echo "Error: Requires 3 arguments"
    echo "Usage: ./rscarpia.sh [--reference hg38|chm13] <vcf_file> <tumor_fastq> <normal_fastq>"
    exit 1
fi

VCF_FILE="${INPUT_FILES[0]}"
TUMOR_FASTQ="${INPUT_FILES[1]}"
NORMAL_FASTQ="${INPUT_FILES[2]}"

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

# Container path
CONTAINER="/home/itoyu8/singularity/rscarpia_0.1.0.sif"

# Reference genome selection
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_FASTA="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_FASTA="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

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
