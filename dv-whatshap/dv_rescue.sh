#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J dv_rescue
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
# Usage: bash scripts/dv_rescue.sh --regions "chr10:17,chr11,chr15" --original-vcf <original.dv.vcf.gz> --type <ont|hifi> [--reference hg38|chm13] [--strict-filter] -d <output_dir> <input.bam>
#
# This script rescues DeepVariant runs that failed with "invalid allele index" errors (chm13).
# Regions are specified manually via --regions argument.
#
# Region format (comma-separated):
#   chr10:17      - chr10 from position 17 to end
#   chr11         - entire chr11 (position 1 to end)
#   chr14:101161446  - chr14 from position 101161446 to end
#
# Workflow:
#   1. Parse --regions argument
#   2. Re-run DeepVariant with --regions for each affected region
#   3. Merge rescued VCFs with original VCF
#   4. Apply strict filter (if --strict-filter)
#   5. Run WhatsHap phasing
#   6. Backup original files (.broken suffix) and copy rescued files to original names
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
REGIONS=""
ORIGINAL_VCF=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --regions)
            REGIONS="$2"
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

if [ -z "$SEQ_TYPE" ] || [ -z "$INPUT_BAM" ] || [ -z "$REGIONS" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$ORIGINAL_VCF" ]; then
    echo "Usage: $0 --regions \"chr10:17,chr11,chr15\" --original-vcf <original.dv.vcf.gz> --type <ont|hifi> [--reference hg38|chm13] [--strict-filter] -d <output_dir> <input.bam>"
    echo ""
    echo "Region format (comma-separated):"
    echo "  chr10:17      - chr10 from position 17 to end"
    echo "  chr11         - entire chr11 (position 1 to end)"
    exit 1
fi

INPUT_BAM=$(realpath "$INPUT_BAM")
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

REFERENCE_FAI="${REFERENCE_GENOME_PATH}.fai"
if [ ! -f "${REFERENCE_FAI}" ]; then
    echo "Error: Reference index not found at ${REFERENCE_FAI}"
    exit 1
fi

if [ ! -f "${ORIGINAL_VCF}" ]; then
    echo "Error: Original VCF not found at ${ORIGINAL_VCF}"
    exit 1
fi

mkdir -p ./log

# Step 1: Parse regions argument and build DeepVariant region strings
echo "=== Parsing regions ==="

REGIONS_FILE="${OUTPUT_DIR}/rescue_regions.txt"
AFFECTED_CHRS_FILE="${OUTPUT_DIR}/affected_chrs.txt"
> "${REGIONS_FILE}"
> "${AFFECTED_CHRS_FILE}"

# Parse comma-separated regions
IFS=',' read -ra REGION_ARRAY <<< "$REGIONS"

for region_spec in "${REGION_ARRAY[@]}"; do
    # Trim whitespace
    region_spec=$(echo "$region_spec" | xargs)

    if [[ "$region_spec" == *":"* ]]; then
        # Format: chr10:17 (from position 17 to end)
        chr=$(echo "$region_spec" | cut -d: -f1)
        start_pos=$(echo "$region_spec" | cut -d: -f2)
    else
        # Format: chr11 (entire chromosome)
        chr="$region_spec"
        start_pos=1
    fi

    # Get chromosome length from .fai
    chr_len=$(awk -v chr="${chr}" '$1 == chr {print $2}' "${REFERENCE_FAI}")

    if [ -z "$chr_len" ]; then
        echo "Error: Chromosome ${chr} not found in reference index"
        exit 1
    fi

    # Build region string: chr:start-end
    echo "${chr}:${start_pos}-${chr_len}" >> "${REGIONS_FILE}"
    echo "${chr}" >> "${AFFECTED_CHRS_FILE}"
    echo "  ${chr}:${start_pos}-${chr_len}"
done

echo ""
echo "Regions to rescue:"
cat "${REGIONS_FILE}"
echo ""
echo "Affected chromosomes:"
cat "${AFFECTED_CHRS_FILE}"

# Step 2: Run DeepVariant for each region
echo "=== Running DeepVariant for rescue regions ==="

RESCUE_VCFS_FILE="${OUTPUT_DIR}/rescue_vcfs.txt"
> "${RESCUE_VCFS_FILE}"

while read region; do
    # Sanitize region name for filename (chr10:17-134758134 -> chr10_17-134758134)
    region_name=$(echo "$region" | sed 's/:/_/g')

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

# Step 3: Concatenate and sort all VCFs
echo "=== Concatenating and sorting VCFs ==="

MERGED_VCF="${OUTPUT_DIR}/rescued.dv.vcf.gz"

"${BCFTOOLS}" concat -a "${ORIGINAL_VCF}" ${RESCUE_VCFS} | \
    "${BCFTOOLS}" sort -m 8G -O z -o "${MERGED_VCF}"

tabix -f -p vcf "${MERGED_VCF}"

# Step 4: Optional strict filter
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

# Step 5: Run WhatsHap
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

# Step 6: Verify output
echo "=== Verification ==="
echo "Chromosomes in rescued VCF:"
"${BCFTOOLS}" index -s "${MERGED_VCF}"

echo ""
echo "Output files:"
ls -lh "${OUTPUT_DIR}"/rescued.*

# Step 7: Backup original files and rename rescued files
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
