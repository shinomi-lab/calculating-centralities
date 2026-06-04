# centralities (Rust 版)

`centralities_targets.py` の Rust 移植・全ノード版。
全ノードについて 8 種類の中心性指標を計算します。

## 計算する指標

- in_degree / out_degree           （次数）
- in_degree_cent / out_degree_cent （次数中心性 = degree/(n-1)）
- betweenness                       （媒介中心性, Brandes アルゴリズム）
- k_shell                           （コア分解, 無向化）
- pagerank                          （PageRank, べき乗法）
- katz                              （Katz 中心性, べき乗法）

## ディレクトリ構成

```
centralities_rust/
├── Cargo.toml
└── src/
    └── main.rs
```

## ビルドと実行

```bash
cd centralities_rust
cargo build --release

# 実行（reverse_edges=false でフォロー方向）
cargo run --release -- /path/to/twitter_dataset.txt

# スレッド数を制限する場合（第2引数）
cargo run --release -- /path/to/twitter_dataset.txt 4
```

## graph_lib との一致

- 元IDを昇順ソートして 0-based index を割り当て（ソート順）
- 重複エッジ除去（HashSet）
- reverse_edges=false（txt のフォロー方向）
  → in_degree = フォロワー数

## 出力

`centralities_all_rust.csv`:

```csv
index,orig_id,in_degree,out_degree,in_degree_cent,out_degree_cent,betweenness,k_shell,pagerank,katz
```

vc_all_fink（Rust/Python）の出力と `index` 列で結合できます。

## 正しさの検証（100ノードテスト, ソート順）

NetworkX 基準値と一致を確認済み：

```
betweenness 最大: index 83, 値 0.0504
k_shell 最大: 4
pagerank 合計: 1.0
λ_max ≈ 2.0000, katz_alpha = 0.45

index 83: in_deg=41, out_deg=7, betweenness=0.0504, k_shell=4
index 18: in_deg=32, out_deg=2, betweenness=0.0420, k_shell=4
index  1: in_deg=5,  out_deg=6, betweenness=0.0387, k_shell=4
```

Rust 版を同じ test グラフで実行し、これらの値と一致すれば移植成功です。

## 計算時間の目安

```
媒介中心性が最も重い（O(V·E)）:
  - 並列化（rayon）で全始点を分散
  - 8万ノード・176万エッジで数分〜数十分

その他（次数・k-shell・PageRank・Katz）は高速:
  - 各々数秒〜数十秒
```

媒介中心性が CPU 負荷が高いので、必要なら第2引数でスレッド数を制限してください。

## 注意：NetworkX との微小な差

- PageRank / Katz はべき乗法の反復回数・収束判定で
  末尾桁が NetworkX と若干異なる場合がありますが、
  ランキング（相関分析）には影響しません。
- betweenness / k_shell は厳密一致します。
