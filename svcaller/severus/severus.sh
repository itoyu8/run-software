#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J severus
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16
# Usage: ./severus.sh --tumor <tumor.bam> [--normal <normal.bam>] --phased-vcf <phased.vcf> --out-dir <output_dir> [--reference hg38|chm13]

# Start time
START_TIME=$(date +%s)

# Parse arguments
TUMOR_BAM=""
NORMAL_BAM=""
PHASED_VCF=""
OUT_DIR=""
REFERENCE_TYPE="hg38"

while [[ $# -gt 0 ]]; do
    case $1 in
        --tumor)
            TUMOR_BAM="$2"
            shift 2
            ;;
        --normal)
            NORMAL_BAM="$2"
            shift 2
            ;;
        --phased-vcf)
            PHASED_VCF="$2"
            shift 2
            ;;
        --out-dir)
            OUT_DIR="$2"
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
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check required arguments
if [ -z "$TUMOR_BAM" ] || [ -z "$PHASED_VCF" ] || [ -z "$OUT_DIR" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: ./severus.sh --tumor <tumor.bam> [--normal <normal.bam>] --phased-vcf <phased.vcf> --out-dir <output_dir> [--reference hg38|chm13]"
    exit 1
fi

# Check if files exist
if [ ! -f "$TUMOR_BAM" ]; then
    echo "Error: Tumor BAM file not found: $TUMOR_BAM"
    exit 1
fi

if [ -n "$NORMAL_BAM" ] && [ ! -f "$NORMAL_BAM" ]; then
    echo "Error: Normal BAM file not found: $NORMAL_BAM"
    exit 1
fi

if [ ! -f "$PHASED_VCF" ]; then
    echo "Error: Phased VCF file not found: $PHASED_VCF"
    exit 1
fi

# Thread definition
THREADS=${SLURM_CPUS_PER_TASK:-16}

# Container path
CONTAINER="/home/itoyu8/singularity/severus_0.1.0.sif"

# Set VNTR bed and PoN paths based on reference type
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    VNTR_BED="/opt/conda/envs/severus_env/lib/python3.10/site-packages/severus-1.6-py3.10.egg/vntrs/chm13v2.0_maskedY_rCRS.trf.bed"
    PON_FILE="/opt/conda/envs/severus_env/lib/python3.10/site-packages/severus-1.6-py3.10.egg/pon/PoN_1000G_chm13.tsv.gz"
else
    VNTR_BED="/opt/conda/envs/severus_env/lib/python3.10/site-packages/severus-1.6-py3.10.egg/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed"
    PON_FILE="/opt/conda/envs/severus_env/lib/python3.10/site-packages/severus-1.6-py3.10.egg/pon/PoN_1000G_hg38.tsv.gz"
fi

# Create output directory
mkdir -p "$OUT_DIR"
mkdir -p log

echo "=== Severus SV calling started at: $(date) ==="
echo "Tumor BAM: $TUMOR_BAM"
if [ -n "$NORMAL_BAM" ]; then
    echo "Normal BAM: $NORMAL_BAM"
    echo "Mode: Tumor/Normal pair"
else
    echo "Mode: Tumor-only (using PoN)"
fi
echo "Phased VCF: $PHASED_VCF"
echo "Output directory: $OUT_DIR"
echo "Reference: $REFERENCE_TYPE"
echo "VNTR bed: $VNTR_BED"
echo "Threads: $THREADS"

# Build severus command
SEVERUS_CMD="singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1:/lustre1 \"$CONTAINER\" severus"
SEVERUS_CMD="$SEVERUS_CMD --target-bam \"$TUMOR_BAM\""

if [ -n "$NORMAL_BAM" ]; then
    # Tumor/Normal mode
    SEVERUS_CMD="$SEVERUS_CMD --control-bam \"$NORMAL_BAM\""
else
    # Tumor-only mode (requires PoN)
    echo "PoN file: $PON_FILE"
    SEVERUS_CMD="$SEVERUS_CMD --PON \"$PON_FILE\""
fi

SEVERUS_CMD="$SEVERUS_CMD --out-dir \"$OUT_DIR\""
SEVERUS_CMD="$SEVERUS_CMD -t $THREADS"
SEVERUS_CMD="$SEVERUS_CMD --phasing-vcf \"$PHASED_VCF\""
SEVERUS_CMD="$SEVERUS_CMD --vntr-bed \"$VNTR_BED\""

# Execute severus
eval $SEVERUS_CMD

# End time
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=== Severus SV calling completed at: $(date) ==="
echo "Total time: ${TOTAL_ELAPSED} seconds"
