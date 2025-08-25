# GLIMPSE2用ファイル

このディレクトリには、GLIMPSE2を使った遺伝子型インピュテーション（欠損遺伝子型の推定）に関連するスクリプトが含まれています。

## スクリプト

### prepare_refpanel.sh

1000ゲノムプロジェクトのVCFファイルから、GLIMPSE2用の参照パネルを作成するSLURMジョブスクリプト。

**処理内容:**
- 染色体1-22とXの各染色体について処理
- マルチアレリック変異をバイアレリックに分割
- SNPsのみを抽出（インデル等を除外）
- バイアレリック変異のみを保持
- BCF形式とサイト情報のみのVCF.gz形式で出力
- 各ファイルのインデックスを作成

**入力:** `shapeit/`ディレクトリ内の1000ゲノムVCFファイル
**出力:** `reference_panel/`ディレクトリ内のBCFファイルとサイトファイル

**リソース要件:** CPU 4コア、メモリ16GB/コア

### make_chunks.sh

GLIMPSE2用の参照パネルからチャンクファイルを作成するSLURMジョブスクリプト。

**処理内容:**
- 染色体1-22とXの各染色体について処理
- GLIMPSE2_chunkを使用してゲノムを小さな領域に分割
- genetic mapファイルを使用して遺伝的距離を考慮
- 各染色体のチャンク情報をテキストファイルで出力

**入力:** `reference_panel/`内のサイトVCFファイルとgenetic mapファイル
**出力:** `chunk_output/`ディレクトリ内のchunks.chr*.txtファイル

**リソース要件:** CPU 4コア、メモリ16GB/コア

### split_reference.sh

参照パネルをGLIMPSE2用バイナリ形式に変換・分割するSLURMジョブスクリプト。

**処理内容:**
- 染色体1-22とXの各染色体について処理
- チャンクファイルの各領域に対してGLIMPSE2_split_referenceを実行
- BCF形式の参照パネルをバイナリ形式に変換
- 各チャンクごとに独立したバイナリファイルを作成
- genetic mapファイルを使用して遺伝的距離を考慮

**入力:** `reference_panel/`内のBCFファイル、チャンクファイル、genetic mapファイル
**出力:** `reference_panel/split/`ディレクトリ内の各チャンク用バイナリファイル

**リソース要件:** CPU 4コア、メモリ16GB/コア

### run_glimpse.sh

GLIMPSE2 phaseとligateを実行してBAMファイルから遺伝子型imputationを行うメインスクリプト。

**処理内容:**
- 全染色体（1-22, X）の全チャンクに対してGLIMPSE2_phaseを実行
- 各染色体内でチャンクをGLIMPSE2_ligateで結合
- 全染色体を統合して最終VCFファイルを作成

**使用方法:**
```bash
sbatch run_glimpse.sh /path/to/sample.bam [output_base_name]
```

- 出力ベース名は省略可能（デフォルト: glimpse_output）
- 最終結果: `{BAMファイルと同じディレクトリ}/{出力ベース名}/glimpse_ligate/sample_all_chromosomes.vcf.gz`

**リソース要件:** CPU 1コア、メモリ64GB

## 実行手順

1. `prepare_refpanel.sh` - リファレンスパネル作成
2. `make_chunks.sh` - チャンク作成
3. `split_reference.sh` - リファレンスパネル分割
4. **`run_glimpse.sh` - BAMファイルからimputation実行**