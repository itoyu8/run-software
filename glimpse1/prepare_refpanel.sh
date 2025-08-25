#!/bin/bash
#SBATCH -J prepare_refpanel_g1
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=16G
#SBATCH -c 4

mkdir -p reference_panel
mkdir -p log

# Process chromosomes 1-22
for CHR in {1..22}; do
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools norm -m -any shapeit/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 |
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools view -m 2 -M 2 -v snps --threads 4 -Ob -o reference_panel/1000GP.chr${CHR}.bcf
    
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools index -f reference_panel/1000GP.chr${CHR}.bcf --threads 4
    
    # Create site files for GL calculation
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools view -G -m 2 -M 2 -v snps reference_panel/1000GP.chr${CHR}.bcf -Oz -o reference_panel/1000GP.chr${CHR}.sites.vcf.gz
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools index -f reference_panel/1000GP.chr${CHR}.sites.vcf.gz
    
    /home/itoyu8/bin/bcftools/bcftools-1.19/bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' reference_panel/1000GP.chr${CHR}.sites.vcf.gz | bgzip -c > reference_panel/1000GP.chr${CHR}.sites.tsv.gz
    tabix -s1 -b2 -e2 reference_panel/1000GP.chr${CHR}.sites.tsv.gz
done

# Process chromosome X separately (different filename format)
/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools norm -m -any shapeit/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz -Ou --threads 4 |
/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools view -m 2 -M 2 -v snps --threads 4 -Ob -o reference_panel/1000GP.chrX.bcf

/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools index -f reference_panel/1000GP.chrX.bcf --threads 4

# Create site files for chromosome X
/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools view -G -m 2 -M 2 -v snps reference_panel/1000GP.chrX.bcf -Oz -o reference_panel/1000GP.chrX.sites.vcf.gz
/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools index -f reference_panel/1000GP.chrX.sites.vcf.gz

/home/itoyu8/bin/bcftools/bcftools-1.19/bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' reference_panel/1000GP.chrX.sites.vcf.gz | bgzip -c > reference_panel/1000GP.chrX.sites.tsv.gz
tabix -s1 -b2 -e2 reference_panel/1000GP.chrX.sites.tsv.gz