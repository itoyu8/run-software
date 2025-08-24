#!/bin/bash
#SBATCH -J split_ref
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=16G
#SBATCH -c 4

CONTAINER="/home/itoyu8/singularity/glimpse_v2.0.0-27-g0919952_20221207.sif"

mkdir -p /home/itoyu8/database/tools/glimpse2/reference_panel/split

# Process chromosomes 1-22
for CHR in {1..22}; do
    echo "Processing chromosome $CHR"
    
    REF=/home/itoyu8/database/tools/glimpse2/reference_panel/1000GP.chr${CHR}.bcf
    MAP=/home/itoyu8/database/tools/glimpse2/maps/genetic_maps.b38/chr${CHR}.b38.gmap.gz
    
    while IFS="" read -r LINE || [ -n "$LINE" ]; do
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        
        singularity exec --bind /home/itoyu8/:/home/itoyu8/ "$CONTAINER" \
            GLIMPSE2_split_reference \
            --reference ${REF} \
            --map ${MAP} \
            --input-region ${IRG} \
            --output-region ${ORG} \
            --output /home/itoyu8/database/tools/glimpse2/reference_panel/split/1000GP.chr${CHR}
    done < /home/itoyu8/database/tools/glimpse2/chunk_output/chunks.chr${CHR}.txt
    
    echo "Completed chromosome $CHR"
done

# Process chromosome X separately
echo "Processing chromosome X"

REF=/home/itoyu8/database/tools/glimpse2/reference_panel/1000GP.chrX.bcf
MAP=/home/itoyu8/database/tools/glimpse2/maps/genetic_maps.b38/chrX.b38.gmap.gz

while IFS="" read -r LINE || [ -n "$LINE" ]; do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    
    singularity exec --bind /home/itoyu8/:/home/itoyu8/ "$CONTAINER" \
        GLIMPSE2_split_reference \
        --reference ${REF} \
        --map ${MAP} \
        --input-region ${IRG} \
        --output-region ${ORG} \
        --output /home/itoyu8/database/tools/glimpse2/reference_panel/split/1000GP.chrX
done < /home/itoyu8/database/tools/glimpse2/chunk_output/chunks.chrX.txt

echo "Completed chromosome X"