# DeepVariant Rescue Guide

CHM13リファレンスでDeepVariantがテロメア領域でクラッシュした際のrescue手順。

## 問題の概要

CHM13リファレンスではテロメア領域に特殊な配列があり、DeepVariantが"invalid allele index"エラーでクラッシュすることがある。クラッシュすると:
1. エラー発生位置以降の染色体データが消失
2. 同じシャード内で順次処理されるはずだった他の染色体も巻き添えで消失

## Step 1: 欠損染色体の確認

```bash
# VCFに含まれる染色体を確認
bcftools index -s <path>/normal.phased.vcf.gz

# 期待される染色体: chr1-22, chrX, chrY, chrM
# 欠損しているものをメモ
```

## Step 2: エラーログの解析

stdoutログ (.o<jobid>) を確認:

```bash
# シャード構造とエラー位置を抽出
grep -E "(Processing region|which is invalid)" /path/to/log/defacto_*.o<jobid>
```

**出力例:**
```
Processing region chr16:0-chr18:80542538
Processing region chr19:0-chr22:51324926
 is [[1]], which is invalid.
 is [[1]], which is invalid.
```

## Step 3: エラー位置の特定

```bash
# エラー周辺のコンテキストを確認
grep -B 10 "which is invalid" /path/to/log/defacto_*.o<jobid>
```

**出力例:**
```
}
end: 30
reference_name: "chr17"
start: 29
 is [[1]], which is invalid.
```

**座標変換:**
- `start: 29` は0-based → 1-basedに変換: **position 30**
- rescueは `chr17:30` から開始

## Step 4: Rescue領域の決定

シャード情報とエラー位置から、rescue対象を決定:

| シャード | エラー位置 | Rescue対象 |
|----------|------------|------------|
| chr16:0-chr18:80542538 | chr17:29 (0-based) | chr17:30, chr18 |
| chr19:0-chr22:51324926 | chr20末端 | chr21, chr22 |

**ルール:**
1. エラー発生染色体: `chrN:<position+1>` (1-based開始位置)
2. 同シャード内でエラー染色体の後にある染色体: 全体 (`chrM`)

## Step 5: dv_rescue.sh の実行

```bash
bash scripts/dv_rescue.sh \
    --regions "chr17:30,chr18,chr21,chr22" \
    --original-vcf <path>/dv_whphase/normal.dv.vcf.gz \
    --type ont \
    --reference chm13 \
    --strict-filter \
    -d <path>/dv_whphase/rescue \
    <path>/normal.bam
```

**--regions フォーマット:**
- `chr17:30` → chr17のposition 30から末端まで (末端は自動取得)
- `chr18` → chr18全体 (position 1から末端まで)

## 典型的なCHM13エラーパターン

よく見られるシャード構成:
- `chr10:0-chr11:135127769` → chr10テロメアでエラー → chr10:16, chr11
- `chr14:0-chr15:99753195` → chr14テロメアでエラー → chr14:101161446, chr15
- `chr16:0-chr18:80542538` → chr17テロメアでエラー → chr17:30, chr18
- `chr19:0-chr22:51324926` → chr20/21テロメアでエラー → chr21, chr22

## Runscript例

```bash
#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J dv_rescue_SAMPLE_chm13
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16

set -euxo pipefail
mkdir -p ./log

DATA_DIR="data/chm13/SAMPLE/ul"

bash scripts/dv_rescue.sh \
    --regions "chr17:30,chr18,chr21,chr22" \
    --original-vcf "${DATA_DIR}/dv_whphase/normal.dv.vcf.gz" \
    --type ont \
    --reference chm13 \
    --strict-filter \
    -d "${DATA_DIR}/dv_whphase/rescue" \
    "${DATA_DIR}/normal.bam"

echo "Exit status: $?"
```

## ファイル配置

- スクリプト本体: `scripts/dv_rescue.sh`
- Runscript: `scripts/runscript/dv_rescue/<SAMPLE>_<modality>_chm13.sh`

## 注意事項

1. **hg38では通常不要** - このエラーはCHM13テロメア特有
2. **座標は1-based** - ログのstart値に+1して指定
3. **シャード巻き添え** - エラー染色体だけでなく同シャード内の後続染色体も確認
