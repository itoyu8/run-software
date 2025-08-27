#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J bwa_samsort
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 32

# Usage: ./bwa_samsort.sh [--reference hg38|chm13] <input_fastq_R1> <output_base_name>

# Parse arguments
INPUT_FASTQ_R1=""
OUTPUT_BASE_NAME=""
REFERENCE_TYPE="hg38"

while [[ $# -gt 0 ]]; do
    case $1 in
        --reference)
            if [ "$2" = "chm13" ]; then
                REFERENCE_TYPE="chm13"
            elif [ "$2" = "hg38" ]; then
                REFERENCE_TYPE="hg38"
            else
                echo "Error: --reference must be either 'hg38' or 'chm13'"
                echo "Usage: $0 [--reference hg38|chm13] <input_fastq_R1> <output_base_name>"
                exit 1
            fi
            shift 2
            ;;
        -*)
            echo "Unknown option $1"
            echo "Usage: $0 [--reference hg38|chm13] <input_fastq_R1> <output_base_name>"
            exit 1
            ;;
        *)
            if [ -z "$INPUT_FASTQ_R1" ]; then
                INPUT_FASTQ_R1="$1"
            elif [ -z "$OUTPUT_BASE_NAME" ]; then
                OUTPUT_BASE_NAME="$1"
            else
                echo "Too many arguments"
                echo "Usage: $0 [--reference hg38|chm13] <input_fastq_R1> <output_base_name>"
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$INPUT_FASTQ_R1" ] || [ -z "$OUTPUT_BASE_NAME" ]; then
    echo "Usage: $0 [--reference hg38|chm13] <input_fastq_R1> <output_base_name>"
    exit 1
fi
FINAL_OUTPUT_BAM="${OUTPUT_BASE_NAME}.bam"
THREADS=${SLURM_CPUS_PER_TASK:-32}

INPUT_DIR=$(dirname "$INPUT_FASTQ_R1")

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi
CONTAINER_PATH="/home/itoyu8/singularity/compat_parabricks-0.2.2.sif"

BAM_FILE="${INPUT_DIR}/${OUTPUT_BASE_NAME}.temp.bam"
SORTED_BAM_FILE="${INPUT_DIR}/${OUTPUT_BASE_NAME}.sorted.bam"
SORTED_BAM_INDEX="${INPUT_DIR}/${OUTPUT_BASE_NAME}.sorted.bam.bai"
MARKDUP_BAM_FILE="${INPUT_DIR}/${FINAL_OUTPUT_BAM}"
MARKDUP_METRICS="${INPUT_DIR}/${OUTPUT_BASE_NAME}.metrics.txt"
MARKDUP_BAM_INDEX="${INPUT_DIR}/${FINAL_OUTPUT_BAM}.bai"

mkdir -p ./log
/home/itoyu8/bin/bwa/bwa-0.7.19/bwa mem -t ${THREADS} \
        -R "@RG\tID:${OUTPUT_BASE_NAME}\tSM:${OUTPUT_BASE_NAME}\tPL:ILLUMINA\tLB:${OUTPUT_BASE_NAME}\tPU:${OUTPUT_BASE_NAME}" \
        ${REFERENCE_GENOME_PATH} ${INPUT_FASTQ_R1} | \
/home/itoyu8/bin/samtools/samtools-1.19/samtools view -bh -@ ${THREADS} - > ${BAM_FILE} \
    || { echo "BWA alignment or SAM to BAM conversion failed"; exit 1; }

/home/itoyu8/bin/samtools/samtools-1.19/samtools sort -@ ${THREADS} -o ${SORTED_BAM_FILE} ${BAM_FILE} \
        || { echo "BAM sorting failed"; exit 1; }

/home/itoyu8/bin/samtools/samtools-1.19/samtools index ${SORTED_BAM_FILE} ${SORTED_BAM_INDEX} \
        || { echo "Sorted BAM indexing failed"; exit 1; }

singularity exec \
        --bind /home/itoyu8/:/home/itoyu8/ \
	--bind /lustre1:/lustre1 \
        ${CONTAINER_PATH} /usr/bin/java \
	-Xmx64G -jar /tools/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar MarkDuplicates \
	--INPUT ${SORTED_BAM_FILE} \
	--OUTPUT ${MARKDUP_BAM_FILE} \
       	--METRICS_FILE ${MARKDUP_METRICS} \
	--CREATE_INDEX true \
	|| { echo "Marking duplicates failed"; exit 1; }

rm ${BAM_FILE}
rm ${SORTED_BAM_FILE}
rm ${SORTED_BAM_INDEX}
