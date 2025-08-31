#!/bin/bash
#SBATCH -J gatk_dict
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1


CONTAINER_PATH="/home/itoyu8/singularity/compat_parabricks-0.2.2.sif"

singularity exec \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1 \
    ${CONTAINER_PATH} /usr/bin/java \
    -Xmx8G -jar /tools/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar CreateSequenceDictionary \
    -R /home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa \
    -O /home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.dict
