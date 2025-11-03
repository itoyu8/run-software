#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J downsample_seqtk
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
# Usage: ./downsample_seqtk.sh --fraction 0.1 <input.fastq.gz> [<input2.fastq.gz> ...]

# Start time
START_TIME=$(date +%s)

# Parse arguments
FRACTION=""
INPUT_FILES=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --fraction)
            FRACTION="$2"
            shift 2
            ;;
        *)
            INPUT_FILES+=("$1")
            shift
            ;;
    esac
done

# Check if fraction is provided
if [ -z "$FRACTION" ]; then
    echo "Error: --fraction option is required"
    echo "Usage: ./downsample.sh --fraction 0.1 <input.fastq.gz> [<input2.fastq.gz> ...]"
    exit 1
fi

# Check if input files are provided
if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "Error: No input files provided"
    echo "Usage: ./downsample.sh --fraction 0.1 <input.fastq.gz> [<input2.fastq.gz> ...]"
    exit 1
fi

# seqtk path
SEQTK="/home/itoyu8/bin/seqtk/seqtk-1.5/seqtk"

# Random seed for reproducibility
SEED=42

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
    OUTPUT_FILE="${BASENAME}.downsampled_${FRACTION}.fastq.gz"

    echo "Processing: $INPUT_FILE"
    echo "Output: $OUTPUT_FILE"
    echo "Fraction: $FRACTION"

    # Start time for this file
    FILE_START=$(date +%s)

    # Run seqtk sample
    ${SEQTK} sample -s${SEED} "$INPUT_FILE" "$FRACTION" | gzip > "$OUTPUT_FILE"

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
