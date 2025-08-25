# Samtools関連ユーティリティスクリプト

BAM/SAMファイルの解析と処理に関するsamtoolsベースのユーティリティ集です。

## スクリプト

### sam_stats.sh
BAMファイルの統計情報を取得

- samtools statsを使用してBAM/SAMファイルの詳細統計を計算
- 複数サンプルの統計情報を一括処理
- リファレンス配列を指定して正確な統計計算

**使用方法**:
```bash
sbatch sam_stats.sh <reference_fasta> <base_directory>
```

**デフォルト**:
- reference_fasta: `/home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta`
- base_directory: `/home/itoyu8/project/pelt_expansion/rawdata/castle/quilt`

**出力**: `*.stats.txt` - 各BAMファイルに対応する統計ファイル

### sam_bedcov.sh
BEDファイル領域のカバレッジ計算

- bedtools makewinfowsで10kb窓を作成
- samtools bedcovで複数BAMファイルのカバレッジを一括計算
- ゲノム全体の窓単位カバレッジマトリックス生成

**使用方法**:
```bash
sbatch sam_bedcov.sh <BAM1> [BAM2] [BAM3] ...
```

**出力**: 最初のBAMファイルと同じディレクトリに生成
- `windows.bed`: 10kb窓BEDファイル
- `bedcov_results.txt`: カバレッジマトリックス

### sam_downsample.sh
BAMファイルのダウンサンプリング

- samtools viewを使用した確率的ダウンサンプリング
- ダウンサンプリング率と出力ファイル名を指定可能
- 自動的にインデックスファイルも生成

**使用方法**:
```bash
sbatch sam_downsample.sh input.bam -r 0.1 --output_name output_prefix
```

**パラメータ**:
- `-r`: ダウンサンプリング率（デフォルト: 0.1）
- `--output_name`: 出力ファイルのプレフィックス名

**出力**: 
- ダウンサンプリングされたBAMファイルとインデックス
- 標準出力に使用したダウンサンプリング率を表示

## 共通仕様

- 全スクリプトSLURMジョブとして実行
- samtools 1.19を使用
- hg38リファレンス配列に対応