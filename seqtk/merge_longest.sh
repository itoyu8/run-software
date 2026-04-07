#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J merge_longest
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 4
# Usage: sbatch seqtk/merge_longest.sh --coverage 60 [-d output_dir] [-o output_name] <ontul.fastq.gz> <ont.fastq.gz>
# Output: <output_dir>/<output_name>.fastq.gz
#
# Merges ONT-UL (all reads) + ONT (longest reads first) to reach target coverage.
# Genome size: 3.1 Gb (fixed)

set -euxo pipefail

SEQTK="/home/itoyu8/bin/seqtk/seqtk-1.5/seqtk"
GENOME_SIZE=3100000000  # 3.1 Gb

# Parse arguments
COVERAGE=""
OUTPUT_DIR="."
OUTPUT_NAME=""
ONTUL_FILE=""
ONT_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --coverage)
            COVERAGE="$2"
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
            if [ -z "$ONTUL_FILE" ]; then
                ONTUL_FILE="$1"
            elif [ -z "$ONT_FILE" ]; then
                ONT_FILE="$1"
            else
                echo "Error: Too many arguments"
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$COVERAGE" ] || [ -z "$ONTUL_FILE" ] || [ -z "$ONT_FILE" ]; then
    echo "Usage: $0 --coverage 60 [-d output_dir] [-o output_name] <ontul.fastq.gz> <ont.fastq.gz>"
    exit 1
fi

ONTUL_FILE=$(realpath "$ONTUL_FILE")
ONT_FILE=$(realpath "$ONT_FILE")
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")

if [ -z "$OUTPUT_NAME" ]; then
    BASENAME=$(basename "$ONTUL_FILE")
    BASENAME="${BASENAME%.gz}"
    BASENAME="${BASENAME%.fastq}"
    BASENAME="${BASENAME%.fq}"
    OUTPUT_NAME="${BASENAME}.${COVERAGE}x"
fi

OUTPUT_FILE="${OUTPUT_DIR}/${OUTPUT_NAME}.fastq.gz"
TARGET_BASES=$(awk "BEGIN {printf \"%.0f\", ${GENOME_SIZE} * ${COVERAGE}}")

mkdir -p ./log

# Step 1: Count ONT-UL total bases
ONTUL_BASES=$(${SEQTK} comp "$ONTUL_FILE" | awk '{sum += $2} END {printf "%.0f", sum}')

# Step 2: Check if ONT-UL already meets target
if [ "$ONTUL_BASES" -ge "$TARGET_BASES" ]; then
    echo "Error: ONT-UL bases (${ONTUL_BASES}) already >= target bases (${TARGET_BASES})"
    echo "ONT-UL coverage alone: $(awk "BEGIN {printf \"%.1f\", ${ONTUL_BASES} / ${GENOME_SIZE}")"
    exit 1
fi

REMAINING_BASES=$((TARGET_BASES - ONTUL_BASES))

# Step 3: Get ONT read names sorted by length (descending), select until target reached
TMPDIR=$(mktemp -d "${OUTPUT_DIR}/merge_longest.XXXXXX")
ONT_READNAMES="${TMPDIR}/ont_selected_reads.txt"
ONT_SUBSET="${TMPDIR}/ont_subset.fastq.gz"

time ${SEQTK} comp "$ONT_FILE" \
    | sort -k2,2nr -S 2G --parallel=4 \
    | awk -v target="$REMAINING_BASES" 'BEGIN {sum=0; done=0} !done {sum += $2; print $1; if (sum >= target) done=1}' \
    > "$ONT_READNAMES"

# Step 4: Extract selected reads from ONT
time ${SEQTK} subseq "$ONT_FILE" "$ONT_READNAMES" | gzip > "$ONT_SUBSET"

# Step 5: Concatenate ONT-UL (all) + ONT (selected) into output
cat "$ONTUL_FILE" "$ONT_SUBSET" > "$OUTPUT_FILE"

# Step 6: Report stats
ONT_SELECTED_BASES=$(${SEQTK} comp "$ONT_SUBSET" | awk '{sum += $2} END {printf "%.0f", sum}')
TOTAL_BASES=$((ONTUL_BASES + ONT_SELECTED_BASES))
ONT_SELECTED_READS=$(wc -l < "$ONT_READNAMES")
ACTUAL_COVERAGE=$(awk "BEGIN {printf \"%.1f\", ${TOTAL_BASES} / ${GENOME_SIZE}}")

echo "=== Summary ==="
echo "ONT-UL bases: ${ONTUL_BASES}"
echo "ONT selected reads: ${ONT_SELECTED_READS}"
echo "ONT selected bases: ${ONT_SELECTED_BASES}"
echo "Total bases: ${TOTAL_BASES}"
echo "Target coverage: ${COVERAGE}x"
echo "Actual coverage: ${ACTUAL_COVERAGE}x"
echo "Output: ${OUTPUT_FILE}"

# Cleanup
rm -rf "$TMPDIR"

echo "Exit status: $?"
