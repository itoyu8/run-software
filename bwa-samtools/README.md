# BWA/Minimap2 + Samtools パイプラインスクリプト

## bwa_samsort.sh
ショートリード向けBWAアライメント → samtools ソート/インデックス → GATK重複マーキングパイプライン  
**使用法**: `sbatch bwa_samsort.sh [--reference hg38|chm13] input_R1.fastq output_base`  
**出力**: `output_base.bam`, `output_base.bam.bai`, `output_base.metrics.txt`  
**オプション**: `--reference hg38` (デフォルト) または `--reference chm13`

## minimap2_samsort.sh
ONT/HiFiロングリード向けMinimap2アライメント → samtools ソート/インデックスパイプライン  
**使用法**: `sbatch minimap2_samsort.sh --type <ont|hifi> [--reference hg38|chm13] input.fastq output_base`  
**出力**: `output_base.bam`, `output_base.bam.bai`  
**オプション**: `--type ont|hifi` (必須), `--reference hg38|chm13` (任意、デフォルト: hg38)

## sam_stats.sh
samtools statsを使用したBAMファイル統計情報取得  
**使用法**: `sbatch sam_stats.sh <reference_fasta> <base_directory>`  
**出力**: `*.stats.txt` ファイル  
**オプション**: デフォルトhg38リファレンスとquiltディレクトリを使用

## sam_bedcov.sh
複数BAMファイルのBED領域カバレッジ計算  
**使用法**: `sbatch sam_bedcov.sh <BAM1> [BAM2] [BAM3] ...`  
**出力**: 最初のBAMディレクトリに `windows.bed`, `bedcov_results.txt`  
**オプション**: 10kb窓を自動作成

## sam_downsample.sh
samtools viewを使用したBAMファイルダウンサンプリング  
**使用法**: `sbatch sam_downsample.sh input.bam -r 0.1 --output_name output_prefix`  
**出力**: ダウンサンプリング済みBAMとインデックスファイル  
**オプション**: `-r` (サンプリング率), `--output_name` (プレフィックス)