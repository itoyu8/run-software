#!/bin/bash
#SBATCH -J sam_stats
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 32

# Usage: ./sam_stats.sh [--reference hg38|chm13] <input_bam_file1> [input_bam_file2] [...]

INPUT_BAMS=()
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
            INPUT_BAMS+=("$1")
            shift
            ;;
    esac
done

if [ ${#INPUT_BAMS[@]} -eq 0 ]; then
    echo "Error: At least one BAM file argument is required"
    echo "Usage: $0 [--reference hg38|chm13] <input_bam_file1> [input_bam_file2] [...]"
    exit 1
fi

THREADS=${SLURM_CPUS_PER_TASK:-32}

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

for INPUT_BAM in "${INPUT_BAMS[@]}"; do
    INPUT_DIR=$(dirname "$INPUT_BAM")
    BASE_NAME=$(basename "$INPUT_BAM" .bam)
    OUTPUT_STATS="${INPUT_DIR}/${BASE_NAME}.stats.txt"
    
    /home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ ${THREADS} --reference "${REFERENCE_GENOME_PATH}" "${INPUT_BAM}" > "${OUTPUT_STATS}"
done