// =====================================================================
// centralities (Rust 版)
// =====================================================================
// centralities_targets.py の Rust 移植・全ノード版。
//
// 計算する指標:
//   - in_degree / out_degree           （次数）
//   - in_degree_cent / out_degree_cent （次数中心性 = degree/(n-1)）
//   - betweenness                       （媒介中心性, Brandes アルゴリズム）
//   - k_shell                           （コア分解, 無向化して計算）
//   - pagerank                          （PageRank, べき乗法）
//   - katz                              （Katz 中心性, べき乗法）
//
// グラフ読み込みは graph_lib と一致:
//   - 元IDを昇順ソートして 0-based index を割り当て（ソート順）
//   - 重複エッジ除去（HashSet）
//   - reverse_edges=false（txt のフォロー方向そのまま）
//     → in_degree = フォロワー数（自然な解釈）
//
// 使い方:
//   cargo run --release -- <txt_path>
// =====================================================================

use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

struct Graph {
    n: usize,
    out_adj: Vec<Vec<usize>>,    // 有向: out_adj[u] = u の後続
    in_adj: Vec<Vec<usize>>,     // 有向: in_adj[v] = v の先行
    undirected_adj: Vec<Vec<usize>>, // 無向版（k-shell 用）
    orig_ids: Vec<i64>,
}

fn load_graph(path: &str, reverse_edges: bool) -> std::io::Result<Graph> {
    println!("[1/4] Reading {} ...", path);
    let t0 = Instant::now();

    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // パス1: 全エッジを読んで元ID集合を収集
    let mut raw_pairs: Vec<(i64, i64)> = Vec::new();
    let mut all_ids: BTreeSet<i64> = BTreeSet::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with('%') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let u_orig: i64 = match parts.next() {
            Some(s) => s.parse().unwrap_or(-1),
            None => continue,
        };
        let v_orig: i64 = match parts.next() {
            Some(s) => s.parse().unwrap_or(-1),
            None => continue,
        };
        if u_orig < 0 || v_orig < 0 {
            continue;
        }
        all_ids.insert(u_orig);
        all_ids.insert(v_orig);
        raw_pairs.push((u_orig, v_orig));
    }

    // ★ graph_lib と同じソート順 index（BTreeSet は昇順済み）
    let orig_ids: Vec<i64> = all_ids.into_iter().collect();
    let mut id_map: HashMap<i64, usize> = HashMap::with_capacity(orig_ids.len());
    for (i, &orig) in orig_ids.iter().enumerate() {
        id_map.insert(orig, i);
    }
    let n = orig_ids.len();

    // パス2: 有向エッジを index に変換（重複除去）
    let mut directed_set: HashSet<(usize, usize)> = HashSet::new();
    for &(u_orig, v_orig) in &raw_pairs {
        let u_idx = id_map[&u_orig];
        let v_idx = id_map[&v_orig];
        if reverse_edges {
            directed_set.insert((v_idx, u_idx));
        } else {
            directed_set.insert((u_idx, v_idx));
        }
    }

    // 有向隣接リスト
    let mut out_adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut in_adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    // 無向隣接リスト（k-shell 用、重複除去）
    let mut undirected_set: Vec<HashSet<usize>> = vec![HashSet::new(); n];

    for &(u, v) in &directed_set {
        out_adj[u].push(v);
        in_adj[v].push(u);
        // 無向版（自己ループは除外）
        if u != v {
            undirected_set[u].insert(v);
            undirected_set[v].insert(u);
        }
    }

    let undirected_adj: Vec<Vec<usize>> = undirected_set
        .into_iter()
        .map(|s| s.into_iter().collect())
        .collect();

    println!("    nodes: {}", n);
    println!("    directed edges (unique): {}", directed_set.len());
    println!("    （ソート順インデックス, reverse_edges={}）", reverse_edges);
    println!("    loaded in {:.1}s", t0.elapsed().as_secs_f64());

    Ok(Graph {
        n,
        out_adj,
        in_adj,
        undirected_adj,
        orig_ids,
    })
}

/// 媒介中心性（Brandes アルゴリズム, 有向版, 並列化）
fn betweenness_centrality(graph: &Graph) -> Vec<f64> {
    println!("    - betweenness (Brandes, 並列) ...");
    let n = graph.n;

    // 各始点からの寄与を並列計算して合算
    let partial: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|s| {
            let mut delta = vec![0.0_f64; n];
            let mut sigma = vec![0.0_f64; n]; // 最短経路数
            let mut dist = vec![-1_i64; n];
            let mut pred: Vec<Vec<usize>> = vec![Vec::new(); n];
            let mut stack: Vec<usize> = Vec::new();
            let mut queue: VecDeque<usize> = VecDeque::new();

            sigma[s] = 1.0;
            dist[s] = 0;
            queue.push_back(s);

            while let Some(v) = queue.pop_front() {
                stack.push(v);
                for &w in &graph.out_adj[v] {
                    // 初到達
                    if dist[w] < 0 {
                        dist[w] = dist[v] + 1;
                        queue.push_back(w);
                    }
                    // 最短経路カウント
                    if dist[w] == dist[v] + 1 {
                        sigma[w] += sigma[v];
                        pred[w].push(v);
                    }
                }
            }

            // 逆順に accumulation
            while let Some(w) = stack.pop() {
                for &v in &pred[w] {
                    if sigma[w] != 0.0 {
                        delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                    }
                }
                // s 自身は除外
                // delta[w] は w が中継した寄与
            }

            // delta を返す（s != w の分のみ）
            let mut contrib = vec![0.0_f64; n];
            for w in 0..n {
                if w != s {
                    contrib[w] = delta[w];
                }
            }
            contrib
        })
        .collect();

    // 合算
    let mut bc = vec![0.0_f64; n];
    for contrib in &partial {
        for i in 0..n {
            bc[i] += contrib[i];
        }
    }

    // 正規化（有向: (n-1)(n-2)）
    if n > 2 {
        let norm = ((n - 1) * (n - 2)) as f64;
        for i in 0..n {
            bc[i] /= norm;
        }
    }

    bc
}

/// k-shell 分解（コア数, 無向グラフ）
fn k_shell(graph: &Graph) -> Vec<usize> {
    println!("    - k-shell decomposition ...");
    let n = graph.n;

    // 次数（無向）
    let mut degree: Vec<usize> = graph.undirected_adj.iter().map(|x| x.len()).collect();
    let mut core = vec![0usize; n];
    let mut removed = vec![false; n];

    // 次数の最大値
    let max_deg = degree.iter().cloned().max().unwrap_or(0);

    // ビンソート方式（O(V+E)）
    // bin[d] = 次数 d のノードのリスト
    let mut k = 0;
    let mut remaining = n;

    // 簡易実装: 各 k について、次数 <= k のノードを繰り返し剥がす
    while remaining > 0 {
        let mut changed = true;
        while changed {
            changed = false;
            for v in 0..n {
                if !removed[v] && degree[v] <= k {
                    core[v] = k;
                    removed[v] = true;
                    remaining -= 1;
                    changed = true;
                    // 隣接ノードの次数を減らす
                    for &u in &graph.undirected_adj[v] {
                        if !removed[u] && degree[u] > 0 {
                            degree[u] -= 1;
                        }
                    }
                }
            }
        }
        k += 1;
        if k > max_deg + 1 {
            break;
        }
    }

    core
}

/// PageRank（べき乗法）
fn pagerank(graph: &Graph, alpha: f64, max_iter: usize, tol: f64) -> Vec<f64> {
    println!("    - PageRank ...");
    let n = graph.n;
    let nf = n as f64;
    let mut pr = vec![1.0 / nf; n];
    let out_deg: Vec<usize> = graph.out_adj.iter().map(|x| x.len()).collect();

    for _ in 0..max_iter {
        let mut new_pr = vec![(1.0 - alpha) / nf; n];

        // dangling（出次数0）ノードの寄与
        let dangling_sum: f64 = (0..n)
            .filter(|&i| out_deg[i] == 0)
            .map(|i| pr[i])
            .sum();
        let dangling_contrib = alpha * dangling_sum / nf;

        for v in 0..n {
            new_pr[v] += dangling_contrib;
        }

        // 通常の寄与
        for u in 0..n {
            if out_deg[u] > 0 {
                let share = alpha * pr[u] / out_deg[u] as f64;
                for &v in &graph.out_adj[u] {
                    new_pr[v] += share;
                }
            }
        }

        // 収束判定
        let diff: f64 = (0..n).map(|i| (new_pr[i] - pr[i]).abs()).sum();
        pr = new_pr;
        if diff < tol {
            break;
        }
    }

    pr
}

/// Katz 中心性（べき乗法）
/// alpha は最大固有値の逆数より小さく設定する必要がある
fn katz_centrality(graph: &Graph, alpha: f64, beta: f64, max_iter: usize, tol: f64) -> Vec<f64> {
    println!("    - Katz centrality (alpha={:.6}) ...", alpha);
    let n = graph.n;
    let mut x = vec![0.0_f64; n];

    for _ in 0..max_iter {
        let mut new_x = vec![beta; n];
        // x_i = alpha * Σ_{j→i} x_j + beta
        // in_adj[i] = i に向かうノード
        for i in 0..n {
            let mut acc = 0.0;
            for &j in &graph.in_adj[i] {
                acc += x[j];
            }
            new_x[i] = alpha * acc + beta;
        }

        let diff: f64 = (0..n).map(|i| (new_x[i] - x[i]).abs()).sum();
        x = new_x;
        if diff < tol {
            break;
        }
    }

    // L2 正規化（NetworkX に合わせる）
    let norm: f64 = x.iter().map(|v| v * v).sum::<f64>().sqrt();
    if norm > 0.0 {
        for v in x.iter_mut() {
            *v /= norm;
        }
    }

    x
}

/// べき乗法で最大固有値を概算（Katz の alpha 決定用）
fn estimate_largest_eigenvalue(graph: &Graph, max_iter: usize) -> f64 {
    let n = graph.n;
    let mut x = vec![1.0_f64 / (n as f64).sqrt(); n];
    let mut lambda = 0.0_f64;

    for _ in 0..max_iter {
        // y = A^T x （in_adj を使う: y_i = Σ_{j→i} x_j）
        let mut y = vec![0.0_f64; n];
        for i in 0..n {
            for &j in &graph.in_adj[i] {
                y[i] += x[j];
            }
        }
        let norm: f64 = y.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm < 1e-12 {
            break;
        }
        for i in 0..n {
            y[i] /= norm;
        }
        lambda = norm;
        x = y;
    }
    lambda
}

fn main() {
    rayon::ThreadPoolBuilder::new()
    .num_threads(3)        // ★ 3スレッドに制限
    .build_global()
    .unwrap();
    let args: Vec<String> = std::env::args().collect();
    let txt_path = args
        .get(1)
        .map(|s| s.as_str())
        .unwrap_or("twitter_dataset.txt");

    // スレッド数制限（任意, 第2引数）
    if let Some(nt) = args.get(2).and_then(|s| s.parse::<usize>().ok()) {
        rayon::ThreadPoolBuilder::new()
            .num_threads(nt)
            .build_global()
            .unwrap();
        println!("使用スレッド数: {}", nt);
    }

    let total_start = Instant::now();

    // ★ centralities は reverse_edges=false（txt のフォロー方向）
    //   → in_degree = フォロワー数（自然な解釈）
    let graph = match load_graph(txt_path, false) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("Error reading {}: {}", txt_path, e);
            std::process::exit(1);
        }
    };

    let n = graph.n;
    println!("\n[2/4] Computing centralities ...");

    // 次数
    let in_degree: Vec<usize> = graph.in_adj.iter().map(|x| x.len()).collect();
    let out_degree: Vec<usize> = graph.out_adj.iter().map(|x| x.len()).collect();

    // 次数中心性 = degree / (n-1)
    let denom = if n > 1 { (n - 1) as f64 } else { 1.0 };
    let in_degree_cent: Vec<f64> = in_degree.iter().map(|&d| d as f64 / denom).collect();
    let out_degree_cent: Vec<f64> = out_degree.iter().map(|&d| d as f64 / denom).collect();

    // 媒介中心性
    let t = Instant::now();
    let betweenness = betweenness_centrality(&graph);
    println!("      betweenness done in {:.1}s", t.elapsed().as_secs_f64());

    // k-shell
    let t = Instant::now();
    let kshell = k_shell(&graph);
    println!("      k-shell done in {:.1}s", t.elapsed().as_secs_f64());

    // PageRank
    let t = Instant::now();
    let pr = pagerank(&graph, 0.85, 200, 1e-9);
    println!("      pagerank done in {:.1}s", t.elapsed().as_secs_f64());

    // Katz（alpha = 0.9 / λ_max）
    let t = Instant::now();
    let lambda_max = estimate_largest_eigenvalue(&graph, 200);
    let alpha = if lambda_max > 0.0 { 0.9 / lambda_max } else { 0.1 };
    println!("      λ_max ≈ {:.4}, Katz alpha = {:.6}", lambda_max, alpha);
    let katz = katz_centrality(&graph, alpha, 1.0, 1000, 1e-9);
    println!("      katz done in {:.1}s", t.elapsed().as_secs_f64());

    // ── CSV 出力 ──
    let output_path = "centralities_all_rust.csv";
    println!("\n[3/4] Writing {} ...", output_path);
    let file = File::create(output_path).expect("cannot create csv");
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "index,orig_id,in_degree,out_degree,in_degree_cent,out_degree_cent,betweenness,k_shell,pagerank,katz"
    )
    .unwrap();
    for i in 0..n {
        writeln!(
            writer,
            "{},{},{},{},{:.6},{:.6},{:.6},{},{:.6},{:.6}",
            i,
            graph.orig_ids[i],
            in_degree[i],
            out_degree[i],
            in_degree_cent[i],
            out_degree_cent[i],
            betweenness[i],
            kshell[i],
            pr[i],
            katz[i],
        )
        .unwrap();
    }
    writer.flush().unwrap();
    println!("    saved {} rows", n);

    // ── 統計と上位表示 ──
    println!("\n[4/4] === Top 10 by each centrality ===");
    let metrics: Vec<(&str, Vec<f64>)> = vec![
        ("out_degree", out_degree.iter().map(|&x| x as f64).collect()),
        ("in_degree", in_degree.iter().map(|&x| x as f64).collect()),
        ("betweenness", betweenness.clone()),
        ("k_shell", kshell.iter().map(|&x| x as f64).collect()),
        ("pagerank", pr.clone()),
        ("katz", katz.clone()),
    ];

    for (name, vals) in &metrics {
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_by(|&a, &b| vals[b].partial_cmp(&vals[a]).unwrap());
        println!("\n  [{}] top 5:", name);
        for &i in idx.iter().take(5) {
            println!(
                "    index={:>6}  orig_id={:>13}  value={:>12.4}",
                i, graph.orig_ids[i], vals[i]
            );
        }
    }

    println!(
        "\n総実行時間: {:.1}s",
        total_start.elapsed().as_secs_f64()
    );
}
