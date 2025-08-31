#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J glimpse2
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=64G
#SBATCH -c 1

# Usage: sbatch run_glimpse.sh /path/to/sample.bam [output_base_name]
# Note: output_base_name can include subdirectories (e.g., "results/glimpse_analysis")

CONTAINER="/home/itoyu8/singularity/glimpse_v2.0.0-27-g0919952_20221207.sif"
BAM=$1
OUTPUT_BASE_NAME=${2:-"glimpse2_output"}

# Set up output directories in the same directory as BAM file
BAM_DIR=$(dirname "$BAM")
OUTPUT_BASE="${BAM_DIR}/${OUTPUT_BASE_NAME}"

mkdir -p "${OUTPUT_BASE}/glimpse_impute"
mkdir -p "${OUTPUT_BASE}/glimpse_ligate"

for CHR in {1..22} X; do
    echo "Processing chromosome $CHR"
    
    REF=/home/itoyu8/database/tools/glimpse2/reference_panel/split/1000GP.chr${CHR}
    
    while IFS="" read -r LINE || [ -n "$LINE" ]; 
    do   
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        CHR_NUM=$(echo ${LINE} | cut -d" " -f2)
        REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
        REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
        OUT="${OUTPUT_BASE}/glimpse_impute/sample_imputed"
        
        singularity exec --bind /home/itoyu8/:/home/itoyu8/ "$CONTAINER" \
            GLIMPSE2_phase \
            --bam-file ${BAM} \
            --reference ${REF}_${CHR_NUM}_${REGS}_${REGE}.bin \
            --output ${OUT}_${CHR_NUM}_${REGS}_${REGE}.bcf
    done < /home/itoyu8/database/tools/glimpse2/chunk_output/chunks.chr${CHR}.txt

    LST="${OUTPUT_BASE}/glimpse_ligate/list.chr${CHR}.txt"
    ls -1v "${OUTPUT_BASE}/glimpse_impute/sample_imputed_chr${CHR}_"*.bcf > ${LST}

    OUT="${OUTPUT_BASE}/glimpse_ligate/sample_chr${CHR}_ligated.bcf"
    singularity exec --bind /home/itoyu8/:/home/itoyu8/ "$CONTAINER" \
        GLIMPSE2_ligate --input ${LST} --output $OUT
    
    echo "Completed chromosome $CHR"
done

# Concatenate all chromosomes
ALL_CHRS=""
for CHR in {1..22} X; do
    ALL_CHRS="$ALL_CHRS ${OUTPUT_BASE}/glimpse_ligate/sample_chr${CHR}_ligated.bcf"
done

/home/itoyu8/bin/bcftools/bcftools-1.22/bcftools concat $ALL_CHRS -Oz -o "${OUTPUT_BASE}/sample.all_chroms.vcf.gz"
/home/itoyu8/bin/bcftools/bcftools-1.22/bcftools index -f "${OUTPUT_BASE}/sample.all_chroms.vcf.gz"