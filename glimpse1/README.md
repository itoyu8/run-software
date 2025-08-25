# GLIMPSE1処理用スクリプト

GLIMPSE1を使った遺伝子型インピュテーションのための前処理スクリプト群です。

## スクリプト

### 1. prepare_refpanel.sh
1000ゲノムプロジェクトのVCFファイルからGLIMPSE1用リファレンスパネルを作成

- マルチアレリック変異をバイアレリックに分割
- バイアレリックSNPsのみを保持
- BCF形式で出力とインデックス作成
- サイト情報ファイル（VCF.gz + TSV.gz）も生成
- 全染色体（chr1-22, chrX）を処理

**入力**: `shapeit/`ディレクトリ内の1000ゲノムVCFファイル
**出力**: `reference_panel/1000GP.chr*.bcf`, `reference_panel/1000GP.chr*.sites.vcf.gz`, `reference_panel/1000GP.chr*.sites.tsv.gz`

### 2. make_chunks.sh
GLIMPSE_chunkを使ってゲノムを適切なチャンクに分割

- 各染色体を最適なサイズのチャンクに分割
- 最小ウィンドウサイズ: 2Mb、最小バッファサイズ: 200kb
- インピュテーション実行時の計算効率と精度のバランスを調整

**入力**: `reference_panel/1000GP.chr*.sites.vcf.gz`
**出力**: `/home/itoyu8/database/tools/glimpse1/chunk_output/chunks.chr*.txt`

### 3. calc_gls.sh  
BAMファイルからGenotype Likelihoods (GLs) を計算

- 事前作成されたサイト情報ファイルを使用
- bcftools mpileup + call でGL計算
- 各染色体で独立したVCFファイルを出力

**使用方法**:
```bash
sbatch calc_gls.sh /path/to/sample.bam [output_folder_name]
```

**出力**: `{BAMディレクトリ}/{フォルダ名}/sample.chr*.vcf.gz`
- デフォルトフォルダ名: `glimpse1_gl`

## 実行手順

1. `prepare_refpanel.sh` - リファレンスパネル作成
2. `make_chunks.sh` - ゲノムのチャンク分割
3. `calc_gls.sh` - BAMファイルからGL計算

全スクリプトSLURMジョブとして実行されます。

## 注意事項

GLIMPSE1では染色体ごとに独立して処理し、複数染色体の統合は推奨されません。