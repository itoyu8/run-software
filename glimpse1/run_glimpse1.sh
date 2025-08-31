#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J glimpse1
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=64G
#SBATCH -c 1

# Usage: sbatch run_glimpse1.sh --type <ONT|illumina|hifi> /path/to/sample.bam [output_base_name]
# Note: output_base_name can include subdirectories (e.g., "results/glimpse_analysis")

CONTAINER="/home/itoyu8/singularity/glimpse1_0.1.0.sif"

# Parse arguments
SEQ_TYPE=""
BAM=""
OUTPUT_BASE_NAME=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --type)
            if [[ "$2" =~ ^(ONT|illumina|hifi)$ ]]; then
                SEQ_TYPE="$2"
            else
                echo "Error: --type must be 'ONT', 'illumina', or 'hifi'"
                exit 1
            fi
            shift 2
            ;;
        *)
            if [ -z "$BAM" ]; then
                BAM="$1"
            elif [ -z "$OUTPUT_BASE_NAME" ]; then
                OUTPUT_BASE_NAME="$1"
            else
                echo "Too many arguments"
                echo "Usage: $0 --type <ONT|illumina|hifi> /path/to/sample.bam [output_base_name]"
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$SEQ_TYPE" ] || [ -z "$BAM" ]; then
    echo "Error: --type and BAM file arguments are required"
    echo "Usage: $0 --type <ONT|illumina|hifi> /path/to/sample.bam [output_base_name]"
    exit 1
fi

OUTPUT_BASE_NAME=${OUTPUT_BASE_NAME:-"glimpse1_output"}

# Set up output directories in the same directory as BAM file
BAM_DIR=$(dirname "$BAM")
OUTPUT_BASE="${BAM_DIR}/${OUTPUT_BASE_NAME}"

mkdir -p "${OUTPUT_BASE}/gl_files"
mkdir -p "${OUTPUT_BASE}/glimpse_impute"
mkdir -p "${OUTPUT_BASE}/glimpse_ligate"
mkdir -p "${OUTPUT_BASE}/glimpse_sample"
mkdir -p log

# Fixed paths
REFPANEL_DIR="/home/itoyu8/database/tools/glimpse1/reference_panel"
MAP_DIR="/home/itoyu8/database/tools/glimpse1/maps/genetic_maps.b38"
CHUNK_DIR="/home/itoyu8/database/tools/glimpse1/chunk_output"
REFGEN="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
BCFTOOLS="/home/itoyu8/bin/bcftools/bcftools-1.22/bcftools"

for CHR in {1..22} X; do
    chr="chr${CHR}"
    
    # Step 1: Calculate genotype likelihoods
    VCF="${REFPANEL_DIR}/1000GP.chr${CHR}.sites.vcf.gz"
    TSV="${REFPANEL_DIR}/1000GP.chr${CHR}.sites.tsv.gz"
    GL_OUT="${OUTPUT_BASE}/gl_files/sample.chr${CHR}.vcf.gz"
    
    if [ "$SEQ_TYPE" = "illumina" ]; then
        "${BCFTOOLS}" mpileup -f "${REFGEN}" -I -E -a 'FORMAT/DP' -T "${VCF}" -r "${chr}" "${BAM}" -Ou |
        "${BCFTOOLS}" call -Aim -C alleles -T "${TSV}" -Oz -o "${GL_OUT}"
    else
        "${BCFTOOLS}" mpileup -f "${REFGEN}" -I -B --max-depth 500 -a 'FORMAT/DP' -T "${VCF}" -r "${chr}" "${BAM}" -Ou |
        "${BCFTOOLS}" call -Aim -C alleles -T "${TSV}" -Oz -o "${GL_OUT}"
    fi
    
    "${BCFTOOLS}" index -f "${GL_OUT}"
    
    # Step 2: GLIMPSE phase imputation
    REF="${REFPANEL_DIR}/1000GP.chr${CHR}.bcf"
    MAP="${MAP_DIR}/chr${CHR}.b38.gmap.gz"
    
    while IFS="" read -r LINE || [ -n "$LINE" ]; 
    do   
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        OUT="${OUTPUT_BASE}/glimpse_impute/sample.chr${CHR}.${ID}.bcf"
        
        singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ "$CONTAINER" \
            GLIMPSE_phase \
            --input "${GL_OUT}" \
            --reference "${REF}" \
            --map "${MAP}" \
            --input-region "${IRG}" \
            --output-region "${ORG}" \
            --output "${OUT}"
        
        "${BCFTOOLS}" index -f "${OUT}"
    done < "${CHUNK_DIR}/chunks.chr${CHR}.txt"

    # Step 3: Ligate chunks
    LST="${OUTPUT_BASE}/glimpse_ligate/list.chr${CHR}.txt"
    ls -1v "${OUTPUT_BASE}/glimpse_impute/sample.chr${CHR}."*.bcf > ${LST}

    OUT="${OUTPUT_BASE}/glimpse_ligate/sample_chr${CHR}_ligated.bcf"
    singularity exec --bind /home/itoyu8/:/home/itoyu8/ "$CONTAINER" \
        GLIMPSE_ligate --input ${LST} --output $OUT
    
    "${BCFTOOLS}" index -f "${OUT}"
    
    # Step 4: Sample haplotypes
    PHASED_OUT="${OUTPUT_BASE}/glimpse_sample/sample_chr${CHR}_phased.bcf"
    singularity exec --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ "$CONTAINER" \
        GLIMPSE_sample --input "${OUT}" --solve --output "${PHASED_OUT}"
    
    "${BCFTOOLS}" index -f "${PHASED_OUT}"
done

# Concatenate all chromosomes
ALL_CHRS=""
for CHR in {1..22} X; do
    ALL_CHRS="$ALL_CHRS ${OUTPUT_BASE}/glimpse_sample/sample_chr${CHR}_phased.bcf"
done

"${BCFTOOLS}" concat $ALL_CHRS -Oz -o "${OUTPUT_BASE}/sample.all_chroms.vcf.gz"
"${BCFTOOLS}" index -f "${OUTPUT_BASE}/sample.all_chroms.vcf.gz"