#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J dv_rescue
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
# Usage: bash scripts/dv_rescue.sh --log <stdout_logfile> --original-vcf <original.dv.vcf.gz> --type <ont|hifi> [--reference hg38|chm13] [--strict-filter] -d <output_dir> <input.bam>
#
# This script rescues DeepVariant runs that failed with "invalid allele index" errors (chm13).
# It parses the stdout log file to find error positions and shard boundaries,
# re-runs DeepVariant for affected regions, and merges the results with the original VCF.
#
# Workflow:
#   1. Parse log for "which is invalid" errors and "Processing region" shard info
#   2. Determine affected chromosomes (error chr + collateral damage chr in same shard)
#   3. Re-run DeepVariant with --regions for each affected region
#   4. Merge rescued VCFs with unaffected chromosomes from original VCF
#   5. Apply strict filter (if --strict-filter)
#   6. Run WhatsHap phasing
#   7. Backup original files (.broken suffix) and copy rescued files to original names
#
# Output:
#   <output_dir>/rescued.dv.vcf.gz          - Merged raw DeepVariant output
#   <output_dir>/rescued.dv.filtered.vcf.gz - Filtered VCF (if --strict-filter)
#   <output_dir>/rescued.phased.vcf.gz      - WhatsHap phased VCF
#   Original files backed up as normal.dv.broken.vcf.gz, etc.

set -euxo pipefail

# Parse arguments
SEQ_TYPE=""
INPUT_BAM=""
OUTPUT_DIR=""
REFERENCE_TYPE="hg38"
STRICT_FILTER=false
LOG_FILE=""
ORIGINAL_VCF=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --log)
            LOG_FILE="$2"
            shift 2
            ;;
        --original-vcf)
            ORIGINAL_VCF="$2"
            shift 2
            ;;
        --type)
            if [ "$2" = "ont" ] || [ "$2" = "hifi" ]; then
                SEQ_TYPE="$2"
            else
                echo "Error: --type must be 'ont' or 'hifi'"
                exit 1
            fi
            shift 2
            ;;
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
        --strict-filter)
            STRICT_FILTER=true
            shift
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

if [ -z "$SEQ_TYPE" ] || [ -z "$INPUT_BAM" ] || [ -z "$LOG_FILE" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$ORIGINAL_VCF" ]; then
    echo "Usage: $0 --log <stdout_logfile> --original-vcf <original.dv.vcf.gz> --type <ont|hifi> [--reference hg38|chm13] [--strict-filter] -d <output_dir> <input.bam>"
    exit 1
fi

INPUT_BAM=$(realpath "$INPUT_BAM")
LOG_FILE=$(realpath "$LOG_FILE")
ORIGINAL_VCF=$(realpath "$ORIGINAL_VCF")
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")

THREADS=${SLURM_CPUS_PER_TASK:-16}

if [ "$SEQ_TYPE" = "ont" ]; then
    DV_MODEL="ONT_R104"
elif [ "$SEQ_TYPE" = "hifi" ]; then
    DV_MODEL="PACBIO"
fi

if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi

DV_WHATSHAP_SIF="/home/itoyu8/singularity/dv-whatshap_0.1.0.sif"
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools"

if [ ! -f "${ORIGINAL_VCF}" ]; then
    echo "Error: Original VCF not found at ${ORIGINAL_VCF}"
    exit 1
fi

mkdir -p ./log

# Step 1: Parse log file for error positions and shard boundaries
echo "=== Parsing log file for errors ==="

# Extract error positions (chr:pos in 0-based)
ERROR_POSITIONS=$(grep -B 3 "which is invalid" "${LOG_FILE}" | \
    awk '/^reference_name:/ {gsub(/"/, ""); chr=$2} /^start:/ {print chr ":" $2}')

if [ -z "${ERROR_POSITIONS}" ]; then
    echo "No 'invalid allele index' errors found in log file."
    exit 0
fi

echo "Error positions found:"
echo "${ERROR_POSITIONS}"

# Extract shard boundaries
SHARD_INFO=$(grep "Processing region" "${LOG_FILE}" | sed 's/.*Processing region //')

# Step 2: Determine regions to re-run
echo "=== Determining regions to re-run ==="

REGIONS_FILE="${OUTPUT_DIR}/rescue_regions.txt"
AFFECTED_CHRS_FILE="${OUTPUT_DIR}/affected_chrs.txt"
> "${REGIONS_FILE}"
> "${AFFECTED_CHRS_FILE}"

# Get chromosome lengths from reference .fai file
REFERENCE_FAI="${REFERENCE_GENOME_PATH}.fai"
if [ ! -f "${REFERENCE_FAI}" ]; then
    echo "Error: Reference index not found at ${REFERENCE_FAI}"
    exit 1
fi

echo "${ERROR_POSITIONS}" | while read err_line; do
    err_chr=$(echo "$err_line" | cut -d: -f1)
    err_pos=$(echo "$err_line" | cut -d: -f2)

    # Find the shard containing this chromosome
    shard=$(echo "${SHARD_INFO}" | grep "^${err_chr}:" | head -1)

    if [ -n "$shard" ]; then
        # Parse shard: chrA:start-chrB:end
        start_chr=$(echo "$shard" | cut -d: -f1)
        end_part=$(echo "$shard" | cut -d- -f2)
        end_chr=$(echo "$end_part" | cut -d: -f1)

        # Get chromosome lengths from .fai file
        err_chr_len=$(awk -v chr="${err_chr}" '$1 == chr {print $2}' "${REFERENCE_FAI}")
        end_chr_len=$(awk -v chr="${end_chr}" '$1 == chr {print $2}' "${REFERENCE_FAI}")

        # Region 1: Error chromosome from (error_pos + 2) in 1-based to end
        # DeepVariant --regions uses 1-based coordinates
        region_start=$((err_pos + 2))
        echo "${err_chr}:${region_start}-${err_chr_len}" >> "${REGIONS_FILE}"
        echo "${err_chr}" >> "${AFFECTED_CHRS_FILE}"

        # Region 2: Collateral damage chromosome (if different)
        if [ "$start_chr" != "$end_chr" ]; then
            echo "${end_chr}:1-${end_chr_len}" >> "${REGIONS_FILE}"
            echo "${end_chr}" >> "${AFFECTED_CHRS_FILE}"
        fi
    fi
done

# Remove duplicates
sort -u "${AFFECTED_CHRS_FILE}" -o "${AFFECTED_CHRS_FILE}"
sort -u "${REGIONS_FILE}" -o "${REGIONS_FILE}"

echo "Regions to re-run:"
cat "${REGIONS_FILE}"

echo "Affected chromosomes:"
cat "${AFFECTED_CHRS_FILE}"

# Step 3: Run DeepVariant for each region
echo "=== Running DeepVariant for rescue regions ==="

RESCUE_VCFS_FILE="${OUTPUT_DIR}/rescue_vcfs.txt"
> "${RESCUE_VCFS_FILE}"

while read region; do
    # Sanitize region name for filename (chr10:17- -> chr10_17)
    region_name=$(echo "$region" | sed 's/:/_/g' | sed 's/-$//')

    DV_RESCUE_OUTPUT="${OUTPUT_DIR}/rescue_${region_name}.vcf.gz"
    DV_TEMP_DIR="${OUTPUT_DIR}/rescue_${region_name}_intermediate"
    mkdir -p "${DV_TEMP_DIR}"

    echo "Running DeepVariant for region: ${region}"

    time singularity exec --nv \
        --bind /home/itoyu8/:/home/itoyu8/ \
        --bind /lustre1:/lustre1/ \
        --bind "${DV_TEMP_DIR}:/tmp" \
        --env TMPDIR=/tmp \
        "${DV_WHATSHAP_SIF}" run_deepvariant \
        --model_type "${DV_MODEL}" \
        --ref "${REFERENCE_GENOME_PATH}" \
        --reads "${INPUT_BAM}" \
        --output_vcf "${DV_RESCUE_OUTPUT}" \
        --intermediate_results_dir "${DV_TEMP_DIR}" \
        --num_shards "${THREADS}" \
        --regions "${region}"

    tabix -f -p vcf "${DV_RESCUE_OUTPUT}"
    echo "${DV_RESCUE_OUTPUT}" >> "${RESCUE_VCFS_FILE}"
done < "${REGIONS_FILE}"

RESCUE_VCFS=$(cat "${RESCUE_VCFS_FILE}" | tr '\n' ' ')

# Step 4: Concatenate and sort all VCFs
# Note: Original VCF contains data up to error positions, rescue VCFs contain data after.
# Simply concatenate and sort - no overlap expected.
echo "=== Concatenating and sorting VCFs ==="

MERGED_VCF="${OUTPUT_DIR}/rescued.dv.vcf.gz"

"${BCFTOOLS}" concat -a "${ORIGINAL_VCF}" ${RESCUE_VCFS} | \
    "${BCFTOOLS}" sort -m 8G -O z -o "${MERGED_VCF}"

tabix -f -p vcf "${MERGED_VCF}"

# Step 6: Optional strict filter
echo "=== Applying filters ==="

FILTERED_OUTPUT="${OUTPUT_DIR}/rescued.dv.filtered.vcf.gz"
WHATSHAP_INPUT="${MERGED_VCF}"

if [ "$STRICT_FILTER" = true ]; then
    MIN_GQ=20
    MIN_VAF=0.3
    MAX_VAF=0.7

    time "${BCFTOOLS}" view \
        -f PASS \
        -m2 -M2 \
        --genotype het \
        "${MERGED_VCF}" | \
    "${BCFTOOLS}" filter \
        -i "FORMAT/GQ >= ${MIN_GQ} && FORMAT/VAF >= ${MIN_VAF} && FORMAT/VAF <= ${MAX_VAF}" \
        -O z -o "${FILTERED_OUTPUT}"

    tabix -f -p vcf "${FILTERED_OUTPUT}"
    WHATSHAP_INPUT="${FILTERED_OUTPUT}"
fi

# Step 7: Run WhatsHap
echo "=== Running WhatsHap ==="

PHASED_OUTPUT="${OUTPUT_DIR}/rescued.phased.vcf.gz"

time singularity exec --nv \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1/ \
    "${DV_WHATSHAP_SIF}" whatshap phase \
    --reference "${REFERENCE_GENOME_PATH}" \
    --ignore-read-groups \
    --distrust-genotypes \
    -o "${PHASED_OUTPUT}" \
    "${WHATSHAP_INPUT}" \
    "${INPUT_BAM}"

tabix -f -p vcf "${PHASED_OUTPUT}"

# Step 8: Verify output
echo "=== Verification ==="
echo "Chromosomes in rescued VCF:"
"${BCFTOOLS}" index -s "${MERGED_VCF}" | cut -f1 | sort -V

echo ""
echo "Output files:"
ls -lh "${OUTPUT_DIR}"/rescued.*

# Step 9: Backup original files and rename rescued files
echo "=== Backing up original files and renaming rescued files ==="

ORIGINAL_DIR=$(dirname "${ORIGINAL_VCF}")

# Backup original files (add .broken suffix)
for ext in vcf.gz vcf.gz.tbi vcf.gz.csi; do
    if [ -f "${ORIGINAL_DIR}/normal.dv.${ext}" ]; then
        mv "${ORIGINAL_DIR}/normal.dv.${ext}" "${ORIGINAL_DIR}/normal.dv.broken.${ext}"
    fi
    if [ -f "${ORIGINAL_DIR}/normal.dv.filtered.${ext}" ]; then
        mv "${ORIGINAL_DIR}/normal.dv.filtered.${ext}" "${ORIGINAL_DIR}/normal.dv.filtered.broken.${ext}"
    fi
    if [ -f "${ORIGINAL_DIR}/normal.phased.${ext}" ]; then
        mv "${ORIGINAL_DIR}/normal.phased.${ext}" "${ORIGINAL_DIR}/normal.phased.broken.${ext}"
    fi
done

# Copy rescued files to original names
cp "${MERGED_VCF}" "${ORIGINAL_DIR}/normal.dv.vcf.gz"
cp "${MERGED_VCF}.tbi" "${ORIGINAL_DIR}/normal.dv.vcf.gz.tbi"

if [ "$STRICT_FILTER" = true ]; then
    cp "${FILTERED_OUTPUT}" "${ORIGINAL_DIR}/normal.dv.filtered.vcf.gz"
    cp "${FILTERED_OUTPUT}.tbi" "${ORIGINAL_DIR}/normal.dv.filtered.vcf.gz.tbi"
fi

cp "${PHASED_OUTPUT}" "${ORIGINAL_DIR}/normal.phased.vcf.gz"
cp "${PHASED_OUTPUT}.tbi" "${ORIGINAL_DIR}/normal.phased.vcf.gz.tbi"

echo ""
echo "Original files backed up with .broken suffix"
echo "Rescued files copied to original names:"
ls -lh "${ORIGINAL_DIR}"/normal.*.vcf.gz

echo "Exit status: $?"
