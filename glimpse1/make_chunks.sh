#!/bin/bash
#SBATCH -J make_chunks_g1
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1

CONTAINER="/home/itoyu8/singularity/glimpse1_0.1.0.sif"

mkdir -p /home/itoyu8/database/tools/glimpse1/chunk_output
mkdir -p log

# Process chromosomes 1-22
for CHR in {1..22}; do
    singularity exec --bind /home/itoyu8/:/home/itoyu8/ "$CONTAINER" \
        GLIMPSE_chunk \
        --input /home/itoyu8/database/tools/glimpse1/reference_panel/1000GP.chr${CHR}.sites.vcf.gz \
        --region chr${CHR} \
        --window-size 2000000 \
        --buffer-size 200000 \
        --output /home/itoyu8/database/tools/glimpse1/chunk_output/chunks.chr${CHR}.txt
done

# Process chromosome X separately
singularity exec --bind /home/itoyu8/:/home/itoyu8/ "$CONTAINER" \
    GLIMPSE_chunk \
    --input /home/itoyu8/database/tools/glimpse1/reference_panel/1000GP.chrX.sites.vcf.gz \
    --region chrX \
    --window-size 2000000 \
    --buffer-size 200000 \
    --output /home/itoyu8/database/tools/glimpse1/chunk_output/chunks.chrX.txt