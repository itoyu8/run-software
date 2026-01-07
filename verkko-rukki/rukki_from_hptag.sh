#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J rukki_hptag
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=7G
#SBATCH -c 56

# Usage: sbatch rukki_from_hptag.sh -b <hptag.bam> -g <p_utg.gfa> -d <output_dir> [-o <output_name>] [--try-fill-bubbles]
# Output: <output_dir>/<output_name>_paths.tsv, <output_dir>/<output_name>_assign.tsv,
#         <output_dir>/<output_name>_MAT.fa, <output_dir>/<output_name>_PAT.fa

set -euxo pipefail

# Default values
HPTAG_BAM=""
INPUT_GFA=""
OUTPUT_DIR=""
OUTPUT_NAME="sample"
TRY_FILL_BUBBLES=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--bam)
            HPTAG_BAM="$2"
            shift 2
            ;;
        -g|--gfa)
            INPUT_GFA="$2"
            shift 2
            ;;
        -d|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-name)
            OUTPUT_NAME="$2"
            shift 2
            ;;
        --try-fill-bubbles)
            TRY_FILL_BUBBLES=true
            shift
            ;;
        *)
            echo "Error: Unknown option $1"
            echo "Usage: $0 -b <hptag.bam> -g <p_utg.gfa> -d <output_dir> [-o <output_name>] [--try-fill-bubbles]"
            exit 1
            ;;
    esac
done

# Validate
if [ -z "$HPTAG_BAM" ] || [ -z "$INPUT_GFA" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Required arguments: -b <hptag.bam> -g <p_utg.gfa> -d <output_dir>"
    exit 1
fi

# Tools
SAMTOOLS="/home/itoyu8/bin/samtools/samtools-1.19/samtools"
YAK="/home/itoyu8/bin/yak/yak-0.1/yak"
RUKKI_CONTAINER="/home/itoyu8/singularity/rukki_0.4.0.sif"
THREADS=${SLURM_CPUS_PER_TASK:-56}

# Setup paths
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
HPTAG_BAM=$(realpath "${HPTAG_BAM}")
INPUT_GFA=$(realpath "${INPUT_GFA}")

# Intermediate files
H1_FASTQ="${OUTPUT_DIR}/${OUTPUT_NAME}_H1.fastq.gz"
H2_FASTQ="${OUTPUT_DIR}/${OUTPUT_NAME}_H2.fastq.gz"
H1_YAK="${OUTPUT_DIR}/${OUTPUT_NAME}_H1.yak"
H2_YAK="${OUTPUT_DIR}/${OUTPUT_NAME}_H2.yak"
NODES_FA="${OUTPUT_DIR}/${OUTPUT_NAME}_nodes.fa"
TRIOEVAL_OUT="${OUTPUT_DIR}/${OUTPUT_NAME}_trioeval.txt"
MARKER_TSV="${OUTPUT_DIR}/${OUTPUT_NAME}_marker_cnts.tsv"
OUT_PATHS="${OUTPUT_DIR}/${OUTPUT_NAME}_paths.tsv"
OUT_ASSIGN="${OUTPUT_DIR}/${OUTPUT_NAME}_assign.tsv"
OUT_FASTA="${OUTPUT_DIR}/${OUTPUT_NAME}_haplotypes.fa"
MAT_FASTA="${OUTPUT_DIR}/${OUTPUT_NAME}_MAT.fa"
PAT_FASTA="${OUTPUT_DIR}/${OUTPUT_NAME}_PAT.fa"

# # Step 1: Extract HP-tagged reads to FASTQ
# time "${SAMTOOLS}" view -@ "${THREADS}" -d HP:1 -u "${HPTAG_BAM}" | \
#     "${SAMTOOLS}" fastq -@ "${THREADS}" - | \
#     gzip -c > "${H1_FASTQ}"
#
# time "${SAMTOOLS}" view -@ "${THREADS}" -d HP:2 -u "${HPTAG_BAM}" | \
#     "${SAMTOOLS}" fastq -@ "${THREADS}" - | \
#     gzip -c > "${H2_FASTQ}"
#
# # Step 2: Build yak k-mer databases
# time "${YAK}" count -b37 -t "${THREADS}" -o "${H1_YAK}" "${H1_FASTQ}"
# time "${YAK}" count -b37 -t "${THREADS}" -o "${H2_YAK}" "${H2_FASTQ}"
#
# # Step 3: Extract node sequences from GFA
# time awk '/^S\t/ { print ">"$2; print $3 }' "${INPUT_GFA}" > "${NODES_FA}"
#
# # Step 4: Run yak trioeval to count markers per node
# time "${YAK}" trioeval -t "${THREADS}" "${H1_YAK}" "${H2_YAK}" "${NODES_FA}" > "${TRIOEVAL_OUT}"
#
# # Step 5: Convert trioeval output to marker counts TSV
# awk '/^S\t/ { print $2"\t"$3"\t"$4 }' "${TRIOEVAL_OUT}" > "${MARKER_TSV}"
#
# # Step 6: Run rukki trio
# RUKKI_ARGS=""
# if [ "$TRY_FILL_BUBBLES" = true ]; then
#     RUKKI_ARGS="--try-fill-bubbles"
# fi
#
# time singularity exec \
#     --bind /home/itoyu8/:/home/itoyu8/ \
#     --bind /lustre1:/lustre1/ \
#     "${RUKKI_CONTAINER}" rukki trio \
#     -g "${INPUT_GFA}" \
#     -m "${MARKER_TSV}" \
#     -p "${OUT_PATHS}" \
#     --final-assign "${OUT_ASSIGN}" \
#     ${RUKKI_ARGS}

# Step 7: Convert rukki paths to FASTA using verkko's get_paths_fasta.py
time singularity exec \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1/ \
    "${RUKKI_CONTAINER}" python3 /opt/verkko_scripts/get_paths_fasta.py \
    "${OUT_PATHS}" \
    "${NODES_FA}" \
    "${INPUT_GFA}" \
    "${OUT_FASTA}"

# Step 8: Split by haplotype (MAT/PAT)
awk '/^>/{p=0} /^>mat_/{p=1} p' "${OUT_FASTA}" > "${MAT_FASTA}"
awk '/^>/{p=0} /^>pat_/{p=1} p' "${OUT_FASTA}" > "${PAT_FASTA}"

echo "Exit status: $?"
