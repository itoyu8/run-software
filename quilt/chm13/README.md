# QUILT2 CHM13 (T2T-native maps + SHAPEIT5 phased panel)

## Overview

T2T-CHM13v2.0 座標系で QUILT2 imputation を実行するためのリソース準備スクリプト群。
遺伝的マップは T2T-native (liftOver ではない)、リファレンスパネルは SHAPEIT5 phased を使用。

## Pipeline

```
1. download_2504_bcf.sh    Mac で実行: 2504版 BCF を HPC にダウンロード
2. create_chunks.sh        HPC: 遺伝的マップからチャンク定義を作成 (1回だけ)
3. split_bcf.sh            HPC: 全ゲノム BCF → 染色体別 VCF に分割 (1回だけ)
4. prepare_reference.sh    HPC: 染色体別 VCF + チャンク → RData に変換 (染色体ごとに実行)
5. (quilt_multi.sh)        BAM → imputed VCF (サンプルごとに実行、別途作成)
```

- Step 1 は Mac から実行 (HPC にインターネット接続なし)
- Step 2-3 は preparation (初回のみ)
- `split_bcf.sh` で全ゲノム BCF を染色体別 VCF に事前分割し、`prepare_reference.sh` はそれを使用
- `prepare_reference.sh` は染色体を引数に取り、22ジョブ並列で投入可能

## Scripts

| Script | Description |
|--------|-------------|
| `download_2504_bcf.sh` | Mac から 2504版 BCF を HPC にダウンロード (~11GB) |
| `create_chunks.sh` | 遺伝的マップからチャンク定義を作成 |
| `split_bcf.sh` | 全ゲノム BCF → 染色体別 VCF に分割 |
| `prepare_reference.sh` | 染色体別 VCF + チャンク → RData に変換 (染色体単位) |

## Source Data

### Genetic Maps: T2T-native recombination maps (Lalli et al. 2025)

- **GitHub**: <https://github.com/JosephLalli/phasing_T2T/tree/main/resources/recombination_maps/t2t_native_scaled_maps>
- **Zenodo**: <https://zenodo.org/records/14891074>
- **Download URL**: `https://raw.githubusercontent.com/JosephLalli/phasing_T2T/main/resources/recombination_maps/t2t_native_scaled_maps/chr{1-22,X}.t2t.scaled.gmap.gz`
- **Files**: `CEU-chr{1-22,X}-final.chm13.txt.gz` (QUILT 形式に変換済み)

**重要**: これは **T2T-CHM13 座標上で直接推定された recombination map** であり、hg38 からの liftOver ではない。ファイル名が `CEU-` で始まるのは QUILT の命名規則に合わせたため。

フォーマット変換:
```bash
# Original: pos<tab>cM/Mb<tab>cM
# → QUILT format: position COMBINED_rate.cM.Mb. Genetic_Map.cM. (space-separated)
curl -sL "${URL}/chr${chr}.t2t.scaled.gmap.gz" | gunzip | \
    awk 'BEGIN{OFS=" "} NR==1{print "position","COMBINED_rate.cM.Mb.","Genetic_Map.cM."} NR>1{print $1,$2,$3}' | \
    gzip > "CEU-chr${chr}-final.chm13.txt.gz"
```

### Reference VCF: 1KGP SHAPEIT5 Phased (2504 unrelated samples)

- **AWS S3 Browse**: <https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/>
- **S3 Direct**: `s3://human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/`
- **File**: `1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf.gz` (~11GB)

**2504 vs 3202**: QUILT2 imputation では統計的独立性のため 2504 (unrelated only) を使用。Scarpia では 3202 (全サンプル) を使用。詳細は `scarpia-paper/scripts/prepare_1kgp/README.md` 参照。

## HPC File Locations

```
/home/itoyu8/database/tools/quilt/chm13/
├── maps/                      # 遺伝的マップ (T2T-native, 配置済み)
├── chunk_output/              # チャンク定義 (create_chunks.sh で作成)
├── per_chr_vcf/               # 染色体別 VCF (split_bcf.sh で作成)
├── prepared_reference/        # 前処理済み RData (prepare_reference.sh で作成)
├── create_chunks.sh
├── split_bcf.sh
└── prepare_reference.sh

/home/itoyu8/database/1000genomes/chm13/
├── whole_genome/              # 全ゲノム BCF (2504版 + 3202版)
└── scarpia_input/             # Scarpia 用 BCF (SNP-only, GT-only)
```

## Usage

```bash
# Step 1: Mac で 2504版 BCF をダウンロード
bash download_2504_bcf.sh

# Step 2: HPC でチャンク作成
sbatch create_chunks.sh

# Step 3: HPC で BCF を染色体別に分割
sbatch split_bcf.sh

# Step 4: HPC で RData 作成 (全染色体を並列投入)
for i in {1..22}; do sbatch prepare_reference.sh chr$i; done
```

## References

- Lalli J, et al. A T2T-CHM13 recombination map and globally diverse haplotype reference panel improves phasing and imputation. *bioRxiv*. 2025. doi: [10.1101/2025.02.24.639687](https://doi.org/10.1101/2025.02.24.639687)
- Davies RW, et al. Rapid genotype imputation from sequence with reference panels. *Nat Genet*. 2021;53(7):1104-1111. doi: [10.1038/s41588-021-00877-0](https://doi.org/10.1038/s41588-021-00877-0)
