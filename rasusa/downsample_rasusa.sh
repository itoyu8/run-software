#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J rasusa_downsample
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 16
# Usage: ./downsample_rasusa.sh --coverage 30 [--reference hg38|chm13] <input.fastq.gz> [<input2.fastq.gz> ...]

# Start time
START_TIME=$(date +%s)

# Thread definition
THREADS=${SLURM_CPUS_PER_TASK:-16}

# Parse arguments
COVERAGE=""
REFERENCE_TYPE="hg38"
INPUT_FILES=()

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
        *)
            INPUT_FILES+=("$1")
            shift
            ;;
    esac
done

# Check if coverage is provided
if [ -z "$COVERAGE" ]; then
    echo "Error: --coverage option is required"
    echo "Usage: ./downsample_rasusa.sh --coverage 30 [--reference hg38|chm13] <input.fastq.gz> [<input2.fastq.gz> ...]"
    exit 1
fi

# Check if input files are provided
if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "Error: No input files provided"
    echo "Usage: ./downsample_rasusa.sh --coverage 30 [--reference hg38|chm13] <input.fastq.gz> [<input2.fastq.gz> ...]"
    exit 1
fi

# Set genome size based on reference type
# hg38/GRCh38: Effective mappable size ~2.91 Gb
# CHM13 T2T: ~3.055 Gb (gaps are filled, so physical/effective size are similar)
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    GENOME_SIZE="3.055gb"
else
    # Use effective genome size for hg38 to match stats calculation
    GENOME_SIZE="2.91gb"
fi

# Singularity container path
CONTAINER="/home/itoyu8/singularity/rasusa_0.1.0.sif"

# Random seed for reproducibility
SEED=42

# Compression level
COMPRESSION_LEVEL=6

# Process each input file
for INPUT_FILE in "${INPUT_FILES[@]}"; do
    if [ ! -f "$INPUT_FILE" ]; then
        echo "Error: File not found: $INPUT_FILE"
        continue
    fi

    # Get base name without extension
    BASENAME=$(basename "$INPUT_FILE")
    if [[ "$BASENAME" == *.fastq.gz ]] || [[ "$BASENAME" == *.fq.gz ]]; then
        BASENAME="${BASENAME%.gz}"
        BASENAME="${BASENAME%.fastq}"
        BASENAME="${BASENAME%.fq}"
    elif [[ "$BASENAME" == *.fastq ]] || [[ "$BASENAME" == *.fq ]]; then
        BASENAME="${BASENAME%.fastq}"
        BASENAME="${BASENAME%.fq}"
    fi

    # Output file in current directory
    OUTPUT_FILE="${BASENAME}.downsampled_${COVERAGE}x.fastq.gz"

    echo "Processing: $INPUT_FILE"
    echo "Reference: $REFERENCE_TYPE (Genome size: $GENOME_SIZE)"
    echo "Coverage: ${COVERAGE}x"
    echo "Output: $OUTPUT_FILE"

    # Start time for this file
    FILE_START=$(date +%s)

    # Run rasusa via Singularity
    singularity exec "$CONTAINER" \
        rasusa reads \
        --coverage "$COVERAGE" \
        --genome-size "$GENOME_SIZE" \
        --output "$OUTPUT_FILE" \
        --seed "$SEED" \
        "$INPUT_FILE"

    # End time for this file
    FILE_END=$(date +%s)
    FILE_ELAPSED=$((FILE_END - FILE_START))

    echo "Completed: $OUTPUT_FILE"
    echo "Time taken: ${FILE_ELAPSED} seconds"
    echo ""
done

# End time
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))

echo "All files processed successfully"
echo "Total time taken: ${TOTAL_ELAPSED} seconds"
