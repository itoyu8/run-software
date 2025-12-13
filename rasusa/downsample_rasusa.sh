#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J rasusa_downsample
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 16
# Usage: sbatch downsample_rasusa.sh --coverage 30 [--reference hg38|chm13] [-o output.fastq.gz] <input.fastq.gz>

# Parse arguments
COVERAGE=""
REFERENCE_TYPE="hg38"
OUTPUT_FILE=""
INPUT_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --coverage|-c)
            COVERAGE="$2"
            shift 2
            ;;
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
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        *)
            INPUT_FILE="$1"
            shift
            ;;
    esac
done

if [ -z "$COVERAGE" ] || [ -z "$INPUT_FILE" ]; then
    echo "Usage: sbatch $0 --coverage 30 [--reference hg38|chm13] [-o output.fastq.gz] <input.fastq.gz>"
    exit 1
fi

# Convert to absolute paths for Singularity compatibility
INPUT_FILE=$(realpath "$INPUT_FILE")

# Set genome size based on reference type
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    GENOME_SIZE="3.055gb"
else
    GENOME_SIZE="2.91gb"
fi

# Set output file name
if [ -z "$OUTPUT_FILE" ]; then
    BASENAME=$(basename "$INPUT_FILE")
    if [[ "$BASENAME" == *.fastq.gz ]] || [[ "$BASENAME" == *.fq.gz ]]; then
        BASENAME="${BASENAME%.gz}"
        BASENAME="${BASENAME%.fastq}"
        BASENAME="${BASENAME%.fq}"
    elif [[ "$BASENAME" == *.fastq ]] || [[ "$BASENAME" == *.fq ]]; then
        BASENAME="${BASENAME%.fastq}"
        BASENAME="${BASENAME%.fq}"
    fi
    OUTPUT_FILE="${BASENAME}.downsampled_${COVERAGE}x.fastq.gz"
fi

# Create log directory
mkdir -p log

# Run rasusa
START_TIME=$(date +%s)

singularity exec /home/itoyu8/singularity/rasusa_0.1.0.sif \
    rasusa reads \
    --coverage "$COVERAGE" \
    --genome-size "$GENOME_SIZE" \
    --output "$OUTPUT_FILE" \
    --seed 42 \
    "$INPUT_FILE"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "Downsampled to ${COVERAGE}x: ${ELAPSED} seconds"
echo "Output: ${OUTPUT_FILE}"
