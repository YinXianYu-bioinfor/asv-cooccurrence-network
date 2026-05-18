#!/usr/bin/env python3
"""
ASV Co-occurrence Network Analysis
====================================
基于 Spearman 秩相关的微生物共现网络构建

用途:
  仅输入特征表（OTU/ASV 丰度表），计算 ASV 之间的两两 Spearman 相关性，
  经 FDR 校正和 LOO 稳定性检验，输出可供 Gephi 可视化的共现网络边列表。

方法原理概述:
  1. ASV 过滤 — 出现频次和总丰度筛选低质量 ASV
  2. CPM 归一化 — 消除测序深度差异（Counts Per Million）
  3. CLR 变换 — 对数比变换，缓解成分数据强迫性负偏差（始终启用）
  4. Spearman 秩相关矩阵 — 向量化计算全部 ASV-ASV 相关系数
     ═ 对原始值做秩变换后计算 Pearson 相关系数 ≡ Spearman ρ
     ═ 利用协方差矩阵一次性计算所有配对，避开 Python 级循环
  5. 提取上三角 — 仅保留 i < j 的半矩阵，跳过对角线自相关和对称冗余
  6. 全局 BH FDR 校正 — 对上三角所有 p 值统一 Benjamini-Hochberg 校正
  7. LOO 稳定性检验 — 逐样本留出、完整重算，评估边的稳健性
  8. 导出 Gephi 兼容边列表（仅保留统计显著 + 稳定性高 + 强相关的边）

输入 (place in ../data/):
  - otutab.txt          : OTUID | Sample1 | Sample2 | ...
  - taxonomy.txt        : OTUID | Phylum | Genus | ...

输出 (written to ../results/):
  - asv_cooccurrence_network_edges_gephi.csv  : 显著 + 稳定 + 强相关边（Gephi 导入用）
  - [可选] asv_cooccurrence_network_edges.csv : 全部 ASV-ASV 配对及统计量（文件较大，默认关闭）
"""

import pandas as pd
import numpy as np
from scipy.stats import t as t_dist
import time as _time
import sys, os

# ======================== 参数配置 ========================
# File paths (relative to project root via ../)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))

INPUT_FILE        = os.path.join(PROJECT_ROOT, "data", "otutab.txt")
TAXONOMY_FILE     = os.path.join(PROJECT_ROOT, "data", "taxonomy.txt")
OUTPUT_EDGES      = os.path.join(PROJECT_ROOT, "results", "asv_cooccurrence_network_edges.csv")
NODE_TAXONOMY_OUT = os.path.join(PROJECT_ROOT, "results", "asv_cooccurrence_network_node_taxonomy.txt")

MIN_OCCURRENCE    = 14       # ASV 至少在 N 个样本中出现
MIN_TOTAL_CPM     = 200      # 所有样本的总 CPM ≥ N

USE_CLR           = True    # 始终启用 Centered Log-Ratio 变换（成分数据偏差校正）
PSEUDOCOUNT       = 0.5     # CLR 中处理零值的伪计数

Q_THRESHOLD       = 0.05    # FDR q 值阈值

DO_LOO            = True    # 是否执行 Leave-One-Out 稳定性检验
LOO_STABILITY     = 0.90    # 边至少要在该比例的 LOO 轮次中保持显著

MIN_ABS_RHO       = 0.70    # Gephi 导出要求的最小 |Spearman ρ|
                            # 0.70 是微生物共现网络的主流阈值（Faust & Raes 2012）

MIN_DEGREE        = 3       # 节点最小度（degree），剔除悬挂节点（degree=1）
                            # 仅连接 1 个其他 ASV 的节点对网络结构贡献有限

TOP_N_NODES       = 30      # 保留连接度最高的 N 个节点用于 Gephi 可视化
                            # 原始网络节点数可能过大（~3000+），无法直接可视化。
                            # 仅保留 Top N 核心节点展示"核心子网络"，是微生物组
                            # 共现网络论文中的通行做法。

OUTPUT_FULL_CSV   = False   # 是否输出完整结果 CSV（7.56M+ 行，文件较大）
# ===============================================================


def log(msg):
    """带时间戳的日志输出"""
    print(f"[{_time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------
#  1. 数据读取
# --------------------------------------------------------------
def read_asv_table(filepath):
    """
    读取特征表（OTU/ASV 丰度表）。

    输入格式: 首列为 OTUID，后续列为各样本的计数值。
    返回:   counts — DataFrame, 行=样本, 列=ASV
            otu_ids — 原始 OTU/ASV 编号列表
    """
    log(f"读取特征表 {filepath} ...")
    otu_df = pd.read_csv(filepath, sep='\t', header=0)
    otu_df.columns = otu_df.columns.str.strip()

    sample_cols = [c for c in otu_df.columns if c != 'OTUID']
    log(f"样本数: {len(sample_cols)}, ASV 数: {otu_df.shape[0]}")

    # 转置为 样本 × ASV
    counts = otu_df[sample_cols].T
    counts.columns = otu_df['OTUID'].values
    counts.index.name = 'SampleID'

    return counts, otu_df['OTUID'].tolist()


# --------------------------------------------------------------
#  2. ASV 过滤
# --------------------------------------------------------------
def filter_asvs(asv_counts, min_occurrence, min_total_cpm):
    """
    基于出现频次和总丰度过滤低质量 ASV。

    与 ASV-ARG 脚本相同的过滤逻辑：先计算 CPM 归一化，再按频次和丰度筛选。
    """
    n_samp = asv_counts.shape[0]
    present = (asv_counts > 0).sum(axis=0)

    lib_sizes = asv_counts.sum(axis=1)
    cpm = asv_counts.div(lib_sizes, axis=0) * 1e6
    total_cpm = cpm.sum(axis=0)

    keep = (present >= min_occurrence) & (total_cpm >= min_total_cpm)
    n_before = len(keep)
    n_after  = keep.sum()
    log(f"ASV 过滤: {n_before} → {n_after} 保留  "
        f"(出现≥{min_occurrence}次, 总CPM≥{min_total_cpm})")
    return cpm.loc[:, keep], keep


# --------------------------------------------------------------
#  3. [可选] CLR 变换 — 缓解成分数据偏差
# --------------------------------------------------------------
def clr_transform(cpm_df, pseudocount=0.5):
    """
    Centered Log-Ratio (CLR) 变换。

    数学原理:
      微生物组测序数据本质上是成分数据（compositional data）——每个样本的
      丰度总和受文库大小约束，ASV 之间并非独立变化。这会导致 Spearman 相关
      中出现强迫性负偏差（compositional negative bias）：即使两个 ASV 在
      绝对意义上不相关，由于"总和固定"的限制，它们的相对丰度会呈现出人为的
      负相关。

      CLR 变换将数据从单形空间（simplex）映射到实数空间：
        CLR(x_i) = ln(x_i / g(x))
      其中 g(x) = (∏_{j=1}^p x_j)^{1/p} 为样本的几何均值。

      变换后的数据消除了总和固定的约束，使后续相关分析更接近真实的协变关系。

    注意事项:
      - 零值无法取对数，故先加伪计数（pseudocount）处理。
      - CLR 对 Spearman 秩相关的影响弱于对 Pearson 的影响（秩不变性），
        但仍有助于减轻成分数据偏差。
      - 对于深度过滤后的 ASV 矩阵（大部分样本有值），CLR 效果最为明显。
    """
    if pseudocount > 0:
        log(f"CLR 变换: 添加伪计数 {pseudocount} 处理零值 ...")
        x = cpm_df + pseudocount
    else:
        x = cpm_df.copy()

    # 逐样本计算几何均值：exp(mean(ln(x)))
    geom_mean = np.exp(np.log(x).mean(axis=1))

    # CLR: ln(x / g(x))
    clr = np.log(x.div(geom_mean, axis=0))
    log("CLR 变换完成。")
    return clr


# --------------------------------------------------------------
#  4. Spearman 相关矩阵（向量化）
# --------------------------------------------------------------
def compute_spearman_matrix(data):
    """
    向量化计算所有 ASV 两两间的 Spearman 秩相关系数。

    数学原理:
      ASV-ASV 共现网络需要计算 p×(p-1)/2 个相关系数。逐对循环效率极低，
      采用矩阵运算实现：
        设 X 为 n×p 的 ASV 秩矩阵（n 样本, p ASV）。

        1. 秩变换: 对每列独立排序并替换为秩（1..n）
        2. 中心化: X_c = X - mean(X)           (n×p)
        3. 协方差: Cov = X_cᵀ @ X_c / (n-1)   (p×p)
        4. 标准差: σ_i = sqrt(ΣX_c² / (n-1))   (p 维向量)
        5. 相关系数矩阵: ρ = D⁻¹ Cov D⁻¹,
           其中 D = diag(σ₀, σ₁, ..., σ_{p-1})   (p×p)
           等价于逐元素: ρ_ij = Cov_ij / (σ_i × σ_j)

        由于 ρ 矩阵对称且对角线为 1（ASV 与自身完全正相关），
        仅需提取上三角部分（i < j）。

    参数:
      data: DataFrame, 行=样本, 列=ASV, 值为丰度

    返回:
      rho_full : np.ndarray (p×p), 完整相关矩阵
      triu_i   : np.ndarray, 上三角行索引
      triu_j   : np.ndarray, 上三角列索引
      rho_vec  : np.ndarray, 上三角的 ρ 值
      p_vals   : np.ndarray, 上三角各配对对应的 p 值
      asv_names: list, ASV 名称列表（与矩阵列顺序一致）
    """
    n, p = data.shape
    total_pairs = p * (p - 1) // 2
    log(f"向量化 Spearman 相关: {p} 个 ASV → {total_pairs} 对 (n={n}) ...")

    # 秩变换：逐列排序 → 替换为秩
    ranks = data.rank(axis=0).values           # n×p

    # 中心化
    centered = ranks - ranks.mean(axis=0)       # n×p

    # 协方差矩阵 (p×p)
    cov = centered.T @ centered / (n - 1)       # p×p

    # 标准差向量 (p,)
    sd = np.sqrt(np.sum(centered ** 2, axis=0) / (n - 1))

    # ρ 矩阵: ρ_ij = Cov_ij / (σ_i × σ_j)
    rho_full = cov / np.outer(sd, sd)           # p×p
    rho_full = np.clip(rho_full, -1.0, 1.0)    # 数值稳定

    # 提取上三角索引（不含对角线）
    triu_i, triu_j = np.triu_indices(p, k=1)   # 均为 (total_pairs,) 数组

    # 上三角 ρ 值
    rho_vec = rho_full[triu_i, triu_j]         # (total_pairs,)

    # t 统计量 → 双尾 p 值
    t_stat = rho_vec * np.sqrt(
        (n - 2) / np.maximum(1 - rho_vec ** 2, 1e-300)
    )
    p_vals = 2 * t_dist.sf(np.abs(t_stat), df=n - 2)
    p_vals = np.maximum(p_vals, 1e-300)

    asv_names = data.columns.tolist()
    log(f"  完整矩阵: {p}×{p}, 上三角: {total_pairs} 对")
    return rho_full, triu_i, triu_j, rho_vec, p_vals, asv_names


# --------------------------------------------------------------
#  5. 全局 Benjamini–Hochberg FDR 校正
# --------------------------------------------------------------
def bh_fdr(pvals):
    """
    全局 Benjamini–Hochberg 错误发现率（FDR）校正。

    与 ASV-ARG 版本不同：ASV-ASV 没有按基因分组的天然依据，
    故对上三角所有 p 值统一校正。
    """
    m = len(pvals)
    if m == 0:
        return np.array([], dtype=float)

    order = np.argsort(pvals)
    sorted_p = pvals[order]

    # 计算每个排序位置的临界值: q_i = p_i × m / (i+1)
    ratios = m / np.arange(1, m + 1, dtype=float)
    q_sorted = sorted_p * ratios

    # 单调性约束：从后向前取最小值
    # q_i = min_{k >= i} (p_k × m / (k+1))
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)

    # 恢复原始顺序
    qvals = np.zeros(m, dtype=float)
    qvals[order] = q_sorted
    return qvals


# --------------------------------------------------------------
#  6. 导出
# --------------------------------------------------------------
def export_edges(rho_full, triu_i, triu_j, rho_vec, p_vals, q_vals,
                 asv_names, sig_mask, stability_arr=None,
                 out_path=OUTPUT_EDGES, stability_cutoff=0.80):
    """
    导出边列表 CSV（Gephi 兼容格式 + 可选完整版）。

    Gephi 网络的三层筛选:
      (a) FDR 显著 (q < Q_THRESHOLD)
      (b) LOO 稳定 (stability ≥ stability_cutoff)  [如启用 LOO]
      (c) 强相关   (|ρ| ≥ MIN_ABS_RHO)

    与 ASV-ARG 版本的关键区别:
      - Source 和 Target 均为 ASV（同一类节点）
      - 不输出 ARG 相关信息
    """
    p = len(asv_names)
    has_loo = stability_arr is not None

    # --- 构建最终筛选掩码 ---
    if has_loo:
        final_mask = sig_mask & (stability_arr >= stability_cutoff) & (np.abs(rho_vec) >= MIN_ABS_RHO)
    else:
        final_mask = sig_mask & (np.abs(rho_vec) >= MIN_ABS_RHO)

    n_sig = sig_mask.sum()
    n_final = final_mask.sum()

    # --- 度筛选：剔除低连接节点 ---
    final_i = triu_i[final_mask]
    final_j = triu_j[final_mask]
    final_rho = rho_vec[final_mask]
    final_q = q_vals[final_mask]

    # 构建临时 DataFrame 以统计节点度
    df = pd.DataFrame({
        'Source':         [asv_names[i] for i in final_i],
        'Target':         [asv_names[j] for j in final_j],
        'Type':           'Undirected',
        'Weight':         np.abs(final_rho),
        'Sign':           np.sign(final_rho).astype(int),
        'Spearman_rho':   final_rho,
        'q_value':        final_q,
    })
    if has_loo:
        df['LOO_stability'] = stability_arr[final_mask]

    # 统计每个节点的度（无向图，Source+Target 合计）
    node_degree = pd.concat([
        df['Source'], df['Target']
    ]).value_counts()
    keep_nodes = set(node_degree[node_degree >= MIN_DEGREE].index)

    n_before_deg = len(df)
    n_nodes_before = node_degree.shape[0]
    df = df[df['Source'].isin(keep_nodes) & df['Target'].isin(keep_nodes)]
    n_removed_edges = n_before_deg - len(df)
    n_removed_nodes = n_nodes_before - len(keep_nodes)
    if n_removed_edges > 0:
        log(f"  度筛选 (degree≥{MIN_DEGREE}): 剔除 {n_removed_nodes} 个低连接节点、"
            f"{n_removed_edges} 条关联边")

    # --- Top-N 核心子网络：保留连接度最高的节点 ---
    node_degree = pd.concat([
        df['Source'], df['Target']
    ]).value_counts()
    top_nodes = set(node_degree.nlargest(min(TOP_N_NODES, len(node_degree))).index)

    n_before_top = len(df)
    n_nodes_before_top = node_degree.shape[0]
    df = df[df['Source'].isin(top_nodes) & df['Target'].isin(top_nodes)]
    n_removed_top_edges = n_before_top - len(df)
    n_removed_top_nodes = n_nodes_before_top - len(top_nodes)
    if n_removed_top_edges > 0:
        log(f"  Top-{TOP_N_NODES} 核心子网: 保留 {len(top_nodes)} 个最高连接度节点、"
            f"{len(df)} 条边（剔除 {n_removed_top_nodes} 个节点、"
            f"{n_removed_top_edges} 条边）")

    # --- Gephi 格式输出 ---
    gephi_path = out_path.replace('.csv', '_gephi.csv')
    df.to_csv(gephi_path, index=False)
    log(f"已写入: {gephi_path}  ({len(df)} 行 — Gephi 兼容格式)")

    # --- [可选] 完整结果输出 ---
    if OUTPUT_FULL_CSV:
        full_i = triu_i
        full_j = triu_j
        full_df = pd.DataFrame({
            'ASV_A':        [asv_names[i] for i in full_i],
            'ASV_B':        [asv_names[j] for j in full_j],
            'Spearman_rho': np.round(rho_vec, 4),
            'p_value':      np.array([f"{v:.4e}" for v in p_vals]),
            'q_value':      np.array([f"{v:.4e}" for v in q_vals]),
            'significant':  sig_mask,
        })
        if has_loo:
            full_df['LOO_stability'] = np.round(stability_arr, 4)

        full_df = full_df.sort_values(['significant', 'Spearman_rho'],
                                      ascending=[False, False])
        full_df.to_csv(out_path, index=False)
        log(f"已写入: {out_path}  ({len(full_df)} 行 — 完整结果)")

    # --- 汇总统计 ---
    final_rho_arr = df['Spearman_rho'].values  # 度筛选后的 ρ
    n_final = len(df)
    n_pos = int((final_rho_arr > 0).sum())
    n_neg = n_final - n_pos
    asv_in_network = len(set(df['Source']).union(df['Target']))
    avg_degree = 2 * n_final / asv_in_network if asv_in_network > 0 else 0

    print()
    print("=" * 56)
    print("  分 析 汇 总")
    print("=" * 56)
    print(f"  参与分析的 ASV 数:       {p}")
    print(f"  ASV 配对总数:            {len(rho_vec)}")
    print(f"  FDR 显著 (q < {Q_THRESHOLD}):       {n_sig}")
    if has_loo:
        n_stable = int((sig_mask & (stability_arr >= stability_cutoff)).sum())
        print(f"  + LOO 稳定 (≥{stability_cutoff:.0%}):   {n_stable}")
    print(f"  ─────────────────────────────────────")
    print(f"  Gephi 网络 (|ρ|≥{MIN_ABS_RHO}, degree≥{MIN_DEGREE}):")
    print(f"    节点: {asv_in_network}")
    print(f"    边:   {n_final}")
    print(f"    正相关: {n_pos} ({100*n_pos/n_final:.1f}%)")
    print(f"    负相关: {n_neg} ({100*n_neg/n_final:.1f}%)")
    print(f"    平均度: {avg_degree:.2f}")
    print("=" * 56)

    return df


# --------------------------------------------------------------
#  7. Leave-One-Out 稳定性检验（矩阵版）
# --------------------------------------------------------------
def loo_stability(data, triu_i, triu_j, rho_full, rho_vec, p_vals, q_vals,
                  asv_names):
    """
    Leave-One-Out 稳定性检验（矩阵版）。

    对已通过 FDR 的边，逐样本留出后完整重算 Spearman 矩阵 + FDR，
    记录每条边在多少比例的 LOO 轮次中仍保持显著且符号不变。

    与 ASV-ARG 版的区别:
      - 每轮重算完整的 p×p 相关矩阵再提取上三角
      - 利用预计算的 triu_i/triu_j 避免重复索引操作
      - 使用全局 FDR 而非逐组 FDR
    """
    n = data.shape[0]
    p = len(asv_names)
    sig_mask = q_vals < Q_THRESHOLD
    n_sig = sig_mask.sum()

    if n_sig == 0:
        log("无显著边，跳过 LOO。")
        return np.zeros(len(rho_vec), dtype=float)

    log(f"Leave-One-Out 稳定性检验: {n} 轮 × {n_sig} 条显著边 ...")

    # 记录显著边 ASV 名称 × 符号 → 索引
    sig_i_names = [asv_names[i] for i in triu_i[sig_mask]]
    sig_j_names = [asv_names[j] for j in triu_j[sig_mask]]
    sig_signs   = np.sign(rho_vec[sig_mask]).astype(int)

    edge2idx = {}
    for k in range(n_sig):
        edge2idx[(sig_i_names[k], sig_j_names[k], sig_signs[k])] = k

    survival = np.zeros(n_sig, dtype=float)
    all_idx = np.arange(n, dtype=int)
    data_vals = data.values  # 预取 numpy 数组加速

    for loo in range(n):
        # 留出第 loo 个样本
        keep = np.delete(all_idx, loo)

        # 重算 Spearman 矩阵（子集版）
        sub = data_vals[keep]                # (n-1) × p
        sub_n = n - 1

        ranks = np.argsort(np.argsort(sub, axis=0), axis=0) + 1  # (n-1)×p
        ranks = ranks.astype(float)

        centered = ranks - ranks.mean(axis=0)
        cov = centered.T @ centered / (sub_n - 1)
        sd = np.sqrt(np.sum(centered ** 2, axis=0) / (sub_n - 1))
        loo_rho = cov / np.outer(sd, sd)
        loo_rho = np.clip(loo_rho, -1.0, 1.0)

        # 上三角 ρ
        loo_rho_vec = loo_rho[triu_i, triu_j]

        # p 值
        t = loo_rho_vec * np.sqrt(
            (sub_n - 2) / np.maximum(1 - loo_rho_vec ** 2, 1e-300)
        )
        loo_p = 2 * t_dist.sf(np.abs(t), df=sub_n - 2)
        loo_p = np.maximum(loo_p, 1e-300)

        # 全局 FDR
        loo_q = np.zeros_like(loo_p)
        order = np.argsort(loo_p)
        ratios = len(loo_p) / np.arange(1, len(loo_p) + 1, dtype=float)
        qs = loo_p[order] * ratios
        qs = np.minimum.accumulate(qs[::-1])[::-1]
        loo_q[order] = np.clip(qs, 0.0, 1.0)

        # 本轮显著边
        loo_sig = loo_q < Q_THRESHOLD
        loo_sig_i_names = [asv_names[i] for i in triu_i[loo_sig]]
        loo_sig_j_names = [asv_names[j] for j in triu_j[loo_sig]]
        loo_sig_signs   = np.sign(loo_rho_vec[loo_sig]).astype(int)

        # 匹配基线显著边
        for k in range(len(loo_sig_i_names)):
            key = (loo_sig_i_names[k], loo_sig_j_names[k], loo_sig_signs[k])
            idx = edge2idx.get(key)
            if idx is not None:
                survival[idx] += 1.0

        if (loo + 1) % 5 == 0 or loo == n - 1:
            log(f"  LOO: {loo + 1}/{n}")

    # 稳定性向量（对应原始 ρ_vec 中显著边的位置）
    stability = survival / n
    n_stable = int((stability >= 0.80).sum())
    log(f"  稳定 (≥80%): {n_stable}/{n_sig} 条边")

    # 构建完整稳定性数组（对所有上三角配对）
    full_stability = np.zeros(len(rho_vec), dtype=float)
    sig_positions = np.where(sig_mask)[0]
    for k, pos in enumerate(sig_positions):
        full_stability[pos] = stability[k]

    return full_stability


# --------------------------------------------------------------
#  7.  节点分类注释
# --------------------------------------------------------------
def read_taxonomy():
    """读取 taxonomy.txt，返回 {ASV: (Phylum, Genus)} 字典。"""
    log(f"读取 {os.path.basename(TAXONOMY_FILE)} ...")
    tax_df = pd.read_csv(TAXONOMY_FILE, sep='\t', header=0)
    tax_df.columns = tax_df.columns.str.strip()
    result = {}
    for _, row in tax_df.iterrows():
        phylum = row['Phylum'] if pd.notna(row['Phylum']) else 'Unassigned'
        genus  = row['Genus']  if pd.notna(row['Genus'])  else 'Unassigned'
        result[row['OTUID']] = (phylum, genus)
    log(f"  已加载 {len(result)} 个 ASV 的分类信息")
    return result


def export_node_taxonomy(gephi_edges, tax_map, out_path):
    """为 Gephi 边列表中的所有节点添加门/属水平注释，写入 txt。"""
    nodes = set(gephi_edges['Source']).union(set(gephi_edges['Target']))
    rows = []
    for node in sorted(nodes):
        phylum, genus = tax_map.get(node, ('Unassigned', 'Unassigned'))
        rows.append({'Node': node, 'Phylum': phylum, 'Genus': genus})
    result = pd.DataFrame(rows)
    result.to_csv(out_path, sep='\t', index=False)
    log(f"已写入: {out_path}  ({len(result)} 个节点)")
    return result


# --------------------------------------------------------------
#  主函数
# --------------------------------------------------------------
def main():
    print("=" * 56)
    print("  ASV 共现网络分析")
    print("  微生物共现模式挖掘")
    print("=" * 56)

    # Ensure output directory exists
    os.makedirs(os.path.join(PROJECT_ROOT, "results"), exist_ok=True)

    # 1 — 读取特征表
    counts, otu_ids = read_asv_table(INPUT_FILE)
    print(f"  原始 ASV 数: {len(otu_ids)}")
    print(f"  样本数:      {counts.shape[0]}")

    # 2 — 过滤 + CPM 归一化
    cpm, _ = filter_asvs(counts, MIN_OCCURRENCE, MIN_TOTAL_CPM)

    # 3 — CLR 变换（成分数据偏差校正，始终启用）
    analysis_data = clr_transform(cpm, PSEUDOCOUNT)
    log("使用 CLR 变换后的数据进行 Spearman 相关分析。")

    # 4 — Spearman 相关矩阵
    rho_full, triu_i, triu_j, rho_vec, p_vals, asv_names = \
        compute_spearman_matrix(analysis_data)

    # 5 — 全局 FDR
    log("执行全局 BH FDR 校正 ...")
    q_vals = bh_fdr(p_vals)
    n_sig = int((q_vals < Q_THRESHOLD).sum())
    log(f"  显著: {n_sig}/{len(p_vals)} (q<{Q_THRESHOLD})")

    # 6 — LOO 稳定性
    stability = None
    if DO_LOO:
        stability = loo_stability(
            analysis_data,
            triu_i, triu_j,
            rho_full, rho_vec, p_vals, q_vals,
            asv_names
        )

    # 7 — 导出
    gephi_df = export_edges(rho_full, triu_i, triu_j, rho_vec, p_vals, q_vals,
                            asv_names, q_vals < Q_THRESHOLD, stability,
                            OUTPUT_EDGES, LOO_STABILITY)

    # 8 — 节点分类注释
    tax_map = read_taxonomy()
    export_node_taxonomy(gephi_df, tax_map, NODE_TAXONOMY_OUT)

    print("\n完成。将 Gephi CSV 导入 Gephi（文件 → 导入电子表格 → 无向图）即可可视化。")


if __name__ == '__main__':
    main()
