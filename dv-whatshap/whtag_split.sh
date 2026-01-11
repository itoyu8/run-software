#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J whtag_split
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 4
# Usage: sbatch whtag_split.sh [--reference hg38|chm13] [-d output_dir] <phased.vcf.gz> <input.bam|cram>
# Output: <output_dir>/<basename>.hptag.bam, <output_dir>/<basename>.h1.bam, <output_dir>/<basename>.h2.bam

# whatshap does not allow multi-threading system (except polyphase)
# CRAM files are converted to BAM first using samtools for faster processing

set -euxo pipefail

# Parse arguments
REFERENCE_TYPE="hg38"
OUTPUT_DIR=""
PHASED_VCF=""
INPUT_FILE=""

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
        -d)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            if [ -z "$PHASED_VCF" ]; then
                PHASED_VCF="$1"
            elif [ -z "$INPUT_FILE" ]; then
                INPUT_FILE="$1"
            else
                echo "Too many arguments"
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$PHASED_VCF" ] || [ -z "$INPUT_FILE" ]; then
    echo "Usage: $0 [--reference hg38|chm13] [-d output_dir] <phased.vcf.gz> <input.bam|cram>"
    exit 1
fi

if [[ "$INPUT_FILE" == *.cram ]]; then
    FILE_BASENAME=$(basename "$INPUT_FILE" .cram)
    FILE_EXT="cram"
elif [[ "$INPUT_FILE" == *.bam ]]; then
    FILE_BASENAME=$(basename "$INPUT_FILE" .bam)
    FILE_EXT="bam"
else
    echo "Error: Input file must be .bam or .cram"
    exit 1
fi

INPUT_DIR=$(dirname "$INPUT_FILE")
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="${INPUT_DIR}"
fi
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REF_FASTA_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REF_FASTA_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

DV_WHATSHAP_SIF="/home/itoyu8/singularity/dv-whatshap_0.1.0.sif"
SAMTOOLS="/home/itoyu8/bin/samtools/samtools-1.19/samtools"
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p ./log

HAPLOTAG_BAM="${OUTPUT_DIR}/${FILE_BASENAME}.hptag.bam"
HAPLOTAG_TSV="${OUTPUT_DIR}/${FILE_BASENAME}.haplotype.tsv.gz"
H1_BAM="${OUTPUT_DIR}/${FILE_BASENAME}.h1.bam"
H2_BAM="${OUTPUT_DIR}/${FILE_BASENAME}.h2.bam"

# Always recreate index to avoid stale .tbi issues
tabix -f -p vcf "${PHASED_VCF}"

# Step 0: Convert CRAM to BAM if input is CRAM
if [ "$FILE_EXT" = "cram" ]; then
    CONVERTED_BAM="${INPUT_DIR}/${FILE_BASENAME}.bam"

    time "${SAMTOOLS}" view -@ "${THREADS}" -b -T "${REF_FASTA_PATH}" -o "${CONVERTED_BAM}" "${INPUT_FILE}"
    "${SAMTOOLS}" index -@ "${THREADS}" "${CONVERTED_BAM}"

    PROCESSING_FILE="${CONVERTED_BAM}"
else
    PROCESSING_FILE="${INPUT_FILE}"
fi

# Step 1: Haplotagging
time singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ "${DV_WHATSHAP_SIF}" \
    whatshap haplotag \
    "${PHASED_VCF}" \
    "${PROCESSING_FILE}" \
    -r "${REF_FASTA_PATH}" \
    -o "${HAPLOTAG_BAM}" \
    --output-threads ${THREADS} \
    --output-haplotag-list "${HAPLOTAG_TSV}" \
    --ignore-read-groups \
    --tag-supplementary \
    --skip-missing-contigs

"${SAMTOOLS}" index "${HAPLOTAG_BAM}"

# Step 2: Decompress TSV for splitting
gunzip -f "${HAPLOTAG_TSV}"
HAPLOTAG_TSV_UNZIPPED="${HAPLOTAG_TSV%.gz}"

# Step 3: Split BAM by haplotype
time singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ "${DV_WHATSHAP_SIF}" \
    whatshap split \
    --output-h1 "${H1_BAM}" \
    --output-h2 "${H2_BAM}" \
    "${HAPLOTAG_BAM}" \
    "${HAPLOTAG_TSV_UNZIPPED}"

"${SAMTOOLS}" index "${H1_BAM}"
"${SAMTOOLS}" index "${H2_BAM}"

rm -f "${HAPLOTAG_TSV_UNZIPPED}"

echo "Exit status: $?"
