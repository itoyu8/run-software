#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J quilt_multi
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
# Usage: sbatch quilt_multi.sh [-d output_dir] <input.bam>
# Output: <output_dir>/quilt.phased.vcf.gz

set -euxo pipefail

# Parse arguments
OUTPUT_DIR="."
INPUT_BAM=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -d)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            if [ -z "$INPUT_BAM" ]; then
                INPUT_BAM="$1"
            else
                echo "Too many arguments"
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$INPUT_BAM" ]; then
    echo "Usage: $0 [-d output_dir] <input.bam>"
    exit 1
fi

INPUT_BAM=$(realpath "$INPUT_BAM")
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")

THREADS=${SLURM_CPUS_PER_TASK:-16}

CONTAINER="/home/itoyu8/singularity/quilt_v0.1.0.sif"
PREPARED_REFERENCE_DIR="/home/itoyu8/database/tools/quilt/output/RData"
CHUNK_DIR="/home/itoyu8/database/tools/quilt/chunk_output"
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.22/bcftools"

NGEN=100
BUFFER=500000

mkdir -p "${OUTPUT_DIR}/impute"
mkdir -p "${OUTPUT_DIR}/ligate"
mkdir -p ./log

sample_name=$(basename "$INPUT_BAM" .bam | sed 's/\..*$//')
bamlist_file="${OUTPUT_DIR}/bamlist.txt"
samplenames_file="${OUTPUT_DIR}/sample_names.txt"

echo "$(realpath "$INPUT_BAM")" > "$bamlist_file"
echo "$sample_name" > "$samplenames_file"

process_chromosome() {
    local chr=$1
    local chunk_file="${CHUNK_DIR}/chunks_${chr}.txt"

    [ ! -f "$chunk_file" ] && return

    tail -n +2 "$chunk_file" | while IFS=$'\t' read -r chunk_id chr_name region; do
        region_start=$(echo $region | sed 's/.*://' | sed 's/-.*//')
        region_end=$(echo $region | sed 's/.*-//')

        prepared_reference="${PREPARED_REFERENCE_DIR}/QUILT_prepared_reference.${chr}.${region_start}.${region_end}.RData"

        [ ! -f "$prepared_reference" ] && continue

        output_vcf="${OUTPUT_DIR}/impute/quilt2.diploid.${chr}.chunk${chunk_id}.vcf.gz"

        singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ "${CONTAINER}" \
            /bin/QUILT2.R \
            --prepared_reference_filename="$prepared_reference" \
            --bamlist="$(realpath "$bamlist_file")" \
            --sampleNames_file="$(realpath "$samplenames_file")" \
            --method=diploid \
            --chr="$chr" \
            --regionStart="$region_start" \
            --regionEnd="$region_end" \
            --nGen="$NGEN" \
            --buffer="$BUFFER" \
            --output_filename="$output_vcf"
    done

    vcf_list="${OUTPUT_DIR}/ligate/list.${chr}.txt"
    ls -1v "${OUTPUT_DIR}/impute/quilt2.diploid.${chr}.chunk"*.vcf.gz > "$vcf_list" 2>/dev/null

    [ ! -s "$vcf_list" ] && return

    ligated_vcf="${OUTPUT_DIR}/ligate/quilt2.diploid.${chr}.ligated.vcf.gz"
    "$BCFTOOLS" concat --ligate --output-type z --output "$ligated_vcf" --file-list "$vcf_list"
}

export -f process_chromosome
export CONTAINER PREPARED_REFERENCE_DIR OUTPUT_DIR CHUNK_DIR NGEN BUFFER BCFTOOLS bamlist_file samplenames_file

time {
    for chr_num in {1..22}; do
        chr="chr${chr_num}"

        process_chromosome "$chr" &

        while (( $(jobs -r | wc -l) >= THREADS )); do
            sleep 1
        done
    done

    wait

    all_vcf_list="${OUTPUT_DIR}/ligate/all_chromosomes.txt"
    ls -1v "${OUTPUT_DIR}/ligate/quilt2.diploid.chr"*.ligated.vcf.gz > "$all_vcf_list" 2>/dev/null

    final_vcf="${OUTPUT_DIR}/quilt.phased.vcf.gz"
    [ -s "$all_vcf_list" ] && "$BCFTOOLS" concat --output-type z --output "$final_vcf" --file-list "$all_vcf_list"
    [ -f "$final_vcf" ] && "$BCFTOOLS" index -f "$final_vcf"
}

echo "Exit status: $?"
