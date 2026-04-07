# QUILT2 hg38 (SHAPEIT2 phased panel)

## Overview

GRCh38 座標系で QUILT2 imputation を実行するためのスクリプト群。
リファレンスパネルには 1000GP SHAPEIT2 phased VCF (旧版) を使用。

## Pipeline

```
1. create_chunks.sh    遺伝的マップからチャンク定義を作成 (1回だけ)
2. prepare_reference.sh チャンク+VCF → RData に変換 (1回だけ)
3. quilt_multi.sh       BAM → imputed VCF (サンプルごとに実行)
```

- Step 1-2 は preparation (初回のみ)。出力は HPC に保存済み。
- Step 3 が実際の imputation run。`quilt.sh` は単一スレッド版 (旧版)。

## Source Data

### Reference VCF: 1KGP SHAPEIT2 Phased (3202 samples, 2020-10 release)

- **URL**: <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/>
- **Files**: `CCDG_14151_B01_GRM_WGS_2020-08-05_chr{1-22}.filtered.shapeit2-duohmm-phased.vcf.gz`
- **chrX**: `CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz`
- **Samples**: 3,202 (全サンプル)
- **Content**: SNV + INDEL + SV (phased)

### Genetic Maps: CEU recombination rate (GRCh38)

- **Source**: QUILT 公式リポジトリ同梱 (<https://github.com/rwdavies/QUILT>)
- **Origin**: CEU recombination rate (build 37) を hg19→hg38 に liftOver
  - Build 37 map: `ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20130507_omni_recombination_rates/CEU_omni_recombination_20130507.tar`
  - LiftOver script: <https://github.com/rwdavies/QUILT/blob/master/scripts/make_b38_recomb_map.R>
- **Files**: `CEU-chr{1-22}-final.b38.txt.gz`

## HPC File Locations

```
/home/itoyu8/database/tools/quilt/hg38/
├── maps/                  # 遺伝的マップ
├── chunk_output/          # チャンク定義 (chunks_chr{1-22}.txt)
├── prepared_reference/    # 前処理済み RData
├── create_chunks.sh       # Prep step 1
└── prepare_reference.sh   # Prep step 2

/home/itoyu8/database/1000genomes/hg38/
└── shapeit2_phased/       # ダウンロードした元 VCF
```

## Usage

```bash
# Imputation (サンプルごとに実行)
sbatch quilt_multi.sh [-d output_dir] <input.bam>
# → <output_dir>/quilt.phased.vcf.gz
```

## Note: SHAPEIT2 vs SHAPEIT4

- **SHAPEIT2** (2020-10): `20201028_3202_phased/` — QUILT2 で使用 (本パイプライン)
- **SHAPEIT4** (2022-04): `20220422_3202_phased_SNV_INDEL_SV/` — Scarpia で使用

QUILT2 は SHAPEIT2 版で prepare_reference 済みのため、そのまま運用。

## References

- Davies RW, et al. Rapid genotype imputation from sequence with reference panels. *Nat Genet*. 2021;53(7):1104-1111. doi: [10.1038/s41588-021-00877-0](https://doi.org/10.1038/s41588-021-00877-0)
- Byrska-Bishop M, et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell*. 2022;185(18):3426-3440.e19. doi: [10.1016/j.cell.2022.08.004](https://doi.org/10.1016/j.cell.2022.08.004)
