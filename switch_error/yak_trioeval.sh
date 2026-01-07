#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J yak_trioeval
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 32
# Usage: bash yak_trioeval.sh -d <output_dir> -o <output_name> <query.fa> <truth_hap1.fa> <truth_hap2.fa>
# Output: <output_dir>/<output_name>.trioeval.txt

set -euxo pipefail

YAK="/home/itoyu8/bin/yak/yak-0.1/yak"
THREADS=${SLURM_CPUS_PER_TASK:-32}

OUTPUT_DIR="."
OUTPUT_NAME="output"

while [[ $# -gt 0 ]]; do
    case $1 in
        -d) OUTPUT_DIR="$2"; shift 2 ;;
        -o) OUTPUT_NAME="$2"; shift 2 ;;
        *)  break ;;
    esac
done

if [[ $# -lt 3 ]]; then
    echo "Error: Requires <query.fa> <truth_hap1.fa> <truth_hap2.fa>"
    exit 1
fi

QUERY_FA="$1"
TRUTH_HAP1_FA="$2"
TRUTH_HAP2_FA="$3"

OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
mkdir -p "${OUTPUT_DIR}"

HAP1_YAK="${OUTPUT_DIR}/${OUTPUT_NAME}_hap1.yak"
HAP2_YAK="${OUTPUT_DIR}/${OUTPUT_NAME}_hap2.yak"
TRIOEVAL_FILE="${OUTPUT_DIR}/${OUTPUT_NAME}.trioeval.txt"

# Count k-mers for truth haplotypes (k=31, bloom filter 37 bits)
time "${YAK}" count -b37 -t "${THREADS}" -o "${HAP1_YAK}" "${TRUTH_HAP1_FA}"
time "${YAK}" count -b37 -t "${THREADS}" -o "${HAP2_YAK}" "${TRUTH_HAP2_FA}"

# Evaluate switch error
time "${YAK}" trioeval -t "${THREADS}" "${HAP1_YAK}" "${HAP2_YAK}" "${QUERY_FA}" > "${TRIOEVAL_FILE}"

echo "Exit status: $?"
