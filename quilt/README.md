# QUILT処理用スクリプト

### 1. create_chunks.sh
genetic mapのファイルから染色体ごとにchunkを作成してtxtファイルで書き出すスクリプト

- QUILTのためにはまずgenetic mapからchunkを作成する必要があります
- quiltのgithubにもgenetic mapは含まれています
- このスクリプトは全染色体（chr1-22）に対してチャンクマッピングを実行し、結果を`chunk_output/chunks_chr*.txt`ファイルに保存します

### 2. prepare_reference.sh
QUILT2_prepare_referenceを全染色体・全チャンクに対して実行するスクリプト

- chunkのデータを使って、population phased vcfをRDataに書き換えて後続の処理をしやすくします
- 各チャンクのregion情報（例：`chr1:1-3997428`）からregionStartとregionEndを自動抽出してQUILT2_prepare_reference.Rに渡します
- 処理結果は`/home/itoyu8/database/tools/quilt/output`に保存されます

### 3. run_quilt.sh
QUILT2 diploid imputationを実行するメインスクリプト

- 全染色体（chr1-22）の全チャンクに対してQUILT2 diploid imputationを実行
- チャンクファイルから座標を読み取り、対応するRDataファイルを使用
- 各染色体内でチャンクを結合し、最終的に全染色体を統合したVCFファイルを生成

## 使用方法

1. まず `create_chunks.sh` を実行してチャンク情報を生成
2. 次に `prepare_reference.sh` を実行してリファレンスデータを準備
3. **`run_quilt.sh` を実行してBAMファイルからimputationを実行**

```bash
sbatch run_quilt.sh /path/to/sample.bam [output_folder_name]
```

- 出力フォルダ名は省略可能（デフォルト: quilt_output）
- 最終結果: `{BAMファイルと同じディレクトリ}/{出力フォルダ名}/quilt2.diploid.all_chromosomes.vcf.gz`

全スクリプトSLURMジョブとして実行されます。