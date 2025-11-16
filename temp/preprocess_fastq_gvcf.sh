#!/bin/bash
#SBATCH -J script_job
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G # Request 8GB RAM per core
#SBATCH -c 16             # Request 16 cores

# Check if three arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fastq_R1> <output_base_name>"
    exit 1
fi

# Assign arguments to variables
# INPUT_FASTQ_R1=$1 # Using absolute paths below
OUTPUT_BASE_NAME=$2
OUTPUT_GVCF="${OUTPUT_BASE_NAME}.g.vcf.gz" # Append .g.vcf.gz to the base name

# Construct absolute paths based on the script's execution directory
EXEC_DIR=$(pwd)
INPUT_FASTQ_R1="${EXEC_DIR}/${1}"

# Define the container path
CONTAINER_PATH="/home/itoyu8/singularity/compat_parabricks-0.2.2.sif"

# Define intermediate filenames using the provided base name
# SAM_FILE is no longer needed as bwa output is piped directly to samtools view
BAM_FILE="./${OUTPUT_BASE_NAME}.bam"
SORTED_BAM_FILE="./${OUTPUT_BASE_NAME}.sorted.bam"
SORTED_BAM_INDEX="./${OUTPUT_BASE_NAME}.sorted.bam.bai"
MARKDUP_BAM_FILE="./${OUTPUT_BASE_NAME}.markdup.bam"
MARKDUP_METRICS="./${OUTPUT_BASE_NAME}.markdup.metrics.txt"
MARKDUP_BAM_INDEX="./${OUTPUT_BASE_NAME}.markdup.bam.bai" # GATK MarkDuplicates creates this by default with CREATE_INDEX=true

echo "Starting BWA alignment and SAM to BAM conversion (piped)..."
# Alignment to the reference genome by bwa, piped directly to samtools view to create BAM
# Using 16 threads for bwa mem (-t 16) and samtools view (-@ 16)
singularity exec \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1 \
    ${CONTAINER_PATH} /bin/bash -c "\
        /tools/bwa-0.7.15/bwa mem -t 16 \
        -R \"@RG\\tID:${OUTPUT_BASE_NAME}\\tSM:${OUTPUT_BASE_NAME}\\tPL:ILLUMINA\\tLB:${OUTPUT_BASE_NAME}\\tPU:${OUTPUT_BASE_NAME}\" \
        /home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa \
        ${INPUT_FASTQ_R1} | \
        /tools/samtools-1.14/samtools view -bh -@ 16 - \
    " > ${BAM_FILE}
echo "BWA alignment and SAM to BAM conversion finished."

echo "Sorting BAM file..."
# samtools sort using 16 threads (-@ 16)
singularity exec \
        --bind /home/itoyu8/:/home/itoyu8/ \
	--bind /lustre1:/lustre1 \
        ${CONTAINER_PATH} \
        /tools/samtools-1.14/samtools sort -@ 16 -o \
        ${SORTED_BAM_FILE} \
        ${BAM_FILE}
echo "BAM sorting finished."

echo "Indexing sorted BAM file..."
# samtools index using 16 threads (-@ 16)
singularity exec \
        --bind /home/itoyu8/:/home/itoyu8/ \
	--bind /lustre1:/lustre1 \
        ${CONTAINER_PATH} \
        /tools/samtools-1.14/samtools index \
        ${SORTED_BAM_FILE} \
        ${SORTED_BAM_INDEX}
echo "Sorted BAM indexing finished."

echo "Marking duplicates..."
# gatk markduplicates
# Note: GATK MarkDuplicates creates an index file (.bai) automatically when CREATE_INDEX=true
singularity exec \
        --bind /home/itoyu8/:/home/itoyu8/ \
	--bind /lustre1:/lustre1 \
        ${CONTAINER_PATH} /usr/bin/java \
	-Xmx64G -jar /tools/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar MarkDuplicates \
	--INPUT ${SORTED_BAM_FILE} \
	--OUTPUT ${MARKDUP_BAM_FILE} \
       	--METRICS_FILE ${MARKDUP_METRICS} \
	--CREATE_INDEX true
echo "Marking duplicates finished."

echo "Running HaplotypeCaller..."
# gatk making gvcf file
# Added --native-pair-hmm-threads 16 for HaplotypeCaller parallelization
singularity exec \
        --bind /home/itoyu8/:/home/itoyu8/ \
	--bind /lustre1:/lustre1 \
        ${CONTAINER_PATH} /usr/bin/java \
        -Xmx64G -jar /tools/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar HaplotypeCaller \
        --reference /home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta \
	--emit-ref-confidence GVCF \
	--input ${MARKDUP_BAM_FILE} \
        --sample-name ${OUTPUT_BASE_NAME} \
      	--output ./${OUTPUT_GVCF} \
        --native-pair-hmm-threads 16 # Use 16 threads for PairHMM
echo "HaplotypeCaller finished."

echo "Indexing GVCF file with tabix..."
# Index the output gvcf file using tabix
singularity exec \
        --bind /home/itoyu8/:/home/itoyu8/ \
	--bind /lustre1:/lustre1 \
        ${CONTAINER_PATH} \
        /tools/htslib-1.14/tabix -p vcf \
        ./${OUTPUT_GVCF}
echo "Tabix indexing finished."

echo "Cleaning up intermediate files..."
# Remove intermediate files
# rm ${SAM_FILE} # SAM file is no longer created
rm ${BAM_FILE}
rm ${SORTED_BAM_FILE}
rm ${SORTED_BAM_INDEX}
# rm ${MARKDUP_BAM_FILE}
# rm ${MARKDUP_METRICS}
# rm ${MARKDUP_BAM_INDEX} # Also remove the index created by MarkDuplicates

echo "Script finished successfully. Output file: ${OUTPUT_GVCF}"
