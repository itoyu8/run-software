#!/bin/bash
#SBATCH -J prepare_refpanel
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=16G
#SBATCH -c 4

mkdir -p reference_panel

# Process chromosomes 1-22
for CHR in {1..22}; do
    echo "Processing chromosome $CHR"
    
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools norm -m -any shapeit/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 |
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools view -m 2 -M 2 -v snps --threads 4 -Ob -o reference_panel/1000GP.chr${CHR}.bcf
    
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools index -f reference_panel/1000GP.chr${CHR}.bcf --threads 4
    
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools view -G -Oz -o reference_panel/1000GP.chr${CHR}.sites.vcf.gz reference_panel/1000GP.chr${CHR}.bcf
    
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools index -f reference_panel/1000GP.chr${CHR}.sites.vcf.gz
    
    echo "Completed chromosome $CHR"
done

# Process chromosome X separately (different filename format)
echo "Processing chromosome X"

/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools norm -m -any shapeit/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz -Ou --threads 4 |
/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools view -m 2 -M 2 -v snps --threads 4 -Ob -o reference_panel/1000GP.chrX.bcf

/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools index -f reference_panel/1000GP.chrX.bcf --threads 4

/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools view -G -Oz -o reference_panel/1000GP.chrX.sites.vcf.gz reference_panel/1000GP.chrX.bcf

/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools index -f reference_panel/1000GP.chrX.sites.vcf.gz

echo "Completed chromosome X"