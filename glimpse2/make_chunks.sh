#!/bin/bash
#SBATCH -J make_chunks
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=16G
#SBATCH -c 4

CONTAINER="/home/itoyu8/singularity/glimpse_v2.0.0-27-g0919952_20221207.sif"

mkdir -p /home/itoyu8/database/tools/glimpse2/chunk_output

# Process chromosomes 1-22
for CHR in {1..22}; do
    echo "Processing chromosome $CHR"
    
    singularity exec --bind /home/itoyu8/:/home/itoyu8/ "$CONTAINER" \
        GLIMPSE2_chunk \
        --input /home/itoyu8/database/tools/glimpse2/reference_panel/1000GP.chr${CHR}.sites.vcf.gz \
        --region chr${CHR} \
        --sequential \
        --output /home/itoyu8/database/tools/glimpse2/chunk_output/chunks.chr${CHR}.txt \
        --map /home/itoyu8/database/tools/glimpse2/maps/genetic_maps.b38/chr${CHR}.b38.gmap.gz
    
    echo "Completed chromosome $CHR"
done

# Process chromosome X separately
echo "Processing chromosome X"

singularity exec --bind /home/itoyu8/:/home/itoyu8/ "$CONTAINER" \
    GLIMPSE2_chunk \
    --input /home/itoyu8/database/tools/glimpse2/reference_panel/1000GP.chrX.sites.vcf.gz \
    --region chrX \
    --sequential \
    --output /home/itoyu8/database/tools/glimpse2/chunk_output/chunks.chrX.txt \
    --map /home/itoyu8/database/tools/glimpse2/maps/genetic_maps.b38/chrX.b38.gmap.gz

echo "Completed chromosome X"