#!/bin/bash
#SBATCH -J sam_downsample
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=64G
#SBATCH -c 1

# Usage: sbatch sam_downsample.sh input.bam -r 0.1 --output_name output_prefix
# ~/bin/samtools/samtools-1.19/samtools stats input.bam で，入力ファイルのcoverageを把握してrを決定すること．
# Default downsampling rate: 0.1 (10%)

BAM_FILE=$1
DOWNSAMPLE_RATE=0.1
OUTPUT_NAME=""

# Parse options
while [[ $# -gt 1 ]]; do
    case $2 in
        -r)
            DOWNSAMPLE_RATE="$3"
            shift 2
            ;;
        --output_name)
            OUTPUT_NAME="$3"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

# Calculate -s value (downsample rate + 42)
S_VALUE=$(echo "42 + $DOWNSAMPLE_RATE" | bc -l)

# Output r value to stdout
echo "r=${DOWNSAMPLE_RATE}"

# Create output filename
BASE_NAME=$(basename "$BAM_FILE" .bam)
DIR_NAME=$(dirname "$BAM_FILE")
if [ -n "$OUTPUT_NAME" ]; then
    OUTPUT_BAM="${DIR_NAME}/${OUTPUT_NAME}.bam"
else
    OUTPUT_BAM="${DIR_NAME}/${BASE_NAME}.downs.bam"
fi

# Downsample BAM
/home/itoyu8/bin/samtools/samtools-1.19/samtools view -s $S_VALUE -b -o "$OUTPUT_BAM" "$BAM_FILE"

# Index the downsampled BAM
/home/itoyu8/bin/samtools/samtools-1.19/samtools index "$OUTPUT_BAM"