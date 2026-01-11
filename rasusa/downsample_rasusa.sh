#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J rasusa_downsample
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 16
# Usage: sbatch downsample_rasusa.sh --coverage 30 [--reference hg38|chm13] [-d output_dir] [-o output_name] <input.fastq.gz>
# Output: <output_dir>/<output_name>.fastq.gz

set -euxo pipefail

# Parse arguments
COVERAGE=""
REFERENCE_TYPE="hg38"
OUTPUT_DIR="."
OUTPUT_NAME=""
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
        -d)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -o)
            OUTPUT_NAME="$2"
            shift 2
            ;;
        -*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            if [ -z "$INPUT_FILE" ]; then
                INPUT_FILE="$1"
            else
                echo "Too many arguments"
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$COVERAGE" ] || [ -z "$INPUT_FILE" ]; then
    echo "Usage: $0 --coverage 30 [--reference hg38|chm13] [-d output_dir] [-o output_name] <input.fastq.gz>"
    exit 1
fi

INPUT_FILE=$(realpath "$INPUT_FILE")
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")

# Set default output name from input filename
if [ -z "$OUTPUT_NAME" ]; then
    BASENAME=$(basename "$INPUT_FILE")
    BASENAME="${BASENAME%.gz}"
    BASENAME="${BASENAME%.fastq}"
    BASENAME="${BASENAME%.fq}"
    OUTPUT_NAME="${BASENAME}.${COVERAGE}x"
fi

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    GENOME_SIZE="3.055gb"
else
    GENOME_SIZE="2.91gb"
fi

OUTPUT_FILE="${OUTPUT_DIR}/${OUTPUT_NAME}.fastq.gz"

mkdir -p ./log

time singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ /home/itoyu8/singularity/rasusa_0.1.0.sif \
    rasusa reads \
    --coverage "$COVERAGE" \
    --genome-size "$GENOME_SIZE" \
    --output "$OUTPUT_FILE" \
    --seed 42 \
    "$INPUT_FILE"

echo "Exit status: $?"
