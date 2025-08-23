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

## 使用方法

1. まず `create_chunks.sh` を実行してチャンク情報を生成
2. 次に `prepare_reference.sh` を実行してリファレンスデータを準備

両スクリプトともSLURMジョブとして実行されます。