# Dipcall - Diploid Variant Calling

## dipcall_male.sh
Phased haplotype assembliesからバリアントコール（男性サンプル、PAR領域対応）
**使用法**: `sbatch dipcall_male.sh --hap1 <hap1.fa> --hap2 <hap2.fa> -o <output_prefix> [--reference hg38|chm13]`
**出力**: `<prefix>.dip.vcf.gz`, `<prefix>.dip.bed`
**オプション**: `--reference hg38|chm13` (デフォルト: hg38)

### CHM13用PARファイルのインストール
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_PAR.bed