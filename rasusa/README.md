# Rasusa - FASTQダウンサンプリング

## 使用方法

```bash
# 基本的な使用（hg38, 30x）
sbatch downsample_rasusa.sh --coverage 30 sample.fastq.gz

# CHM13リファレンスで50x
sbatch downsample_rasusa.sh -c 50 --reference chm13 sample.fastq.gz
```

出力: `{basename}.downsampled_{coverage}x.fastq.gz`
