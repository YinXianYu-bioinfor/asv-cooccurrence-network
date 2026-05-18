#!/usr/bin/env python3
"""
ASV Co-occurrence Network Analysis
====================================
Building microbial co-occurrence networks based on Spearman rank correlation.

Purpose:
  Takes only a feature table (OTU/ASV abundance table), computes pairwise
  Spearman correlations between ASVs, applies FDR correction and LOO stability
  testing, and exports Gephi-compatible co-occurrence network edge lists.

Method overview:
  1. ASV filtering — occurrence frequency and total abundance filters
  2. CPM normalization — Counts Per Million to correct for sequencing depth
  3. CLR transformation — log-ratio transform to mitigate compositional bias (always on)
  4. Spearman rank correlation matrix — vectorized all-pairs ASV–ASV computation
     ═ rank-transform raw values, then compute Pearson correlation ≡ Spearman ρ
     ═ single covariance matrix for all pairs, no Python-level loops
  5. Upper-triangle extraction — keep only i < j, skip diagonal self-correlations and symmetry
  6. Global BH FDR correction — single Benjamini-Hochberg on all upper-triangle p-values
  7. LOO stability test — full recomputation with each sample held out
  8. Gephi-compatible edge list export (statistically significant + stable + strong edges only)

Input (place in ../data/):
  - otutab.txt          : OTUID | Sample1 | Sample2 | ...
  - taxonomy.txt        : OTUID | Phylum | Genus | ...

Output (written to ../results/):
  - asv_cooccurrence_network_edges_gephi.csv  : significant + stable + strong edges (for Gephi)
  - [optional] asv_cooccurrence_network_edges.csv : all ASV-ASV pairs with statistics (large, off by default)
"""

import pandas as pd
import numpy as np
from scipy.stats import t as t_dist
import time as _time
import sys, os

# ======================== CONFIGURATION ========================
# File paths (relative to project root via ../)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))

INPUT_FILE        = os.path.join(PROJECT_ROOT, "data", "otutab.txt")
TAXONOMY_FILE     = os.path.join(PROJECT_ROOT, "data", "taxonomy.txt")
OUTPUT_EDGES      = os.path.join(PROJECT_ROOT, "results", "asv_cooccurrence_network_edges.csv")
NODE_TAXONOMY_OUT = os.path.join(PROJECT_ROOT, "results", "asv_cooccurrence_network_node_taxonomy.txt")

MIN_OCCURRENCE    = 14       # ASV must appear in >= N samples
MIN_TOTAL_CPM     = 200      # total CPM across all samples >= N

USE_CLR           = True    # always enable Centered Log-Ratio transform (compositional bias correction)
PSEUDOCOUNT       = 0.5     # pseudocount for zero-value handling in CLR

Q_THRESHOLD       = 0.05    # FDR q-value cutoff

DO_LOO            = True    # run Leave-One-Out stability test
LOO_STABILITY     = 0.90    # edge must survive >= this fraction of LOO iterations

MIN_ABS_RHO       = 0.70    # minimum |Spearman ρ| for Gephi export
                            # 0.70 is the mainstream threshold for microbial co-occurrence networks (Faust & Raes 2012)

MIN_DEGREE        = 3       # minimum node degree; removes dangling nodes (degree=1)
                            # nodes connected to only 1 other ASV contribute little to network structure

TOP_N_NODES       = 30      # keep top-N highest-degree nodes for Gephi visualization
                            # raw networks can have ~3000+ nodes, which are impractical to visualize.
                            # Retaining only the top-N core nodes to show the "core subnetwork"
                            # is standard practice in microbiome co-occurrence network papers.

OUTPUT_FULL_CSV   = False   # whether to output the full results CSV (7.56M+ rows, large file)
# ===============================================================


def log(msg):
    """Timestamped log output."""
    print(f"[{_time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------
#  1. Data reading
# --------------------------------------------------------------
def read_asv_table(filepath):
    """
    Read the feature table (OTU/ASV abundance table).

    Input format: first column is OTUID, subsequent columns are sample counts.
    Returns:   counts — DataFrame, rows=samples, columns=ASVs
               otu_ids — list of original OTU/ASV IDs
    """
    log(f"Reading feature table {filepath} ...")
    otu_df = pd.read_csv(filepath, sep='\t', header=0)
    otu_df.columns = otu_df.columns.str.strip()

    sample_cols = [c for c in otu_df.columns if c != 'OTUID']
    log(f"Samples: {len(sample_cols)}, ASVs: {otu_df.shape[0]}")

    # Transpose to samples × ASVs
    counts = otu_df[sample_cols].T
    counts.columns = otu_df['OTUID'].values
    counts.index.name = 'SampleID'

    return counts, otu_df['OTUID'].tolist()


# --------------------------------------------------------------
#  2. ASV filtering
# --------------------------------------------------------------
def filter_asvs(asv_counts, min_occurrence, min_total_cpm):
    """
    Filter low-quality ASVs based on occurrence frequency and total abundance.

    Same filtering logic as the ASV-ARG script: CPM normalization first,
    then filter by occurrence and abundance.
    """
    n_samp = asv_counts.shape[0]
    present = (asv_counts > 0).sum(axis=0)

    lib_sizes = asv_counts.sum(axis=1)
    cpm = asv_counts.div(lib_sizes, axis=0) * 1e6
    total_cpm = cpm.sum(axis=0)

    keep = (present >= min_occurrence) & (total_cpm >= min_total_cpm)
    n_before = len(keep)
    n_after  = keep.sum()
    log(f"ASV filter: {n_before} → {n_after} retained  "
        f"(occurrence≥{min_occurrence}, total_CPM≥{min_total_cpm})")
    return cpm.loc[:, keep], keep


# --------------------------------------------------------------
#  3. [Optional] CLR transform — mitigate compositional bias
# --------------------------------------------------------------
def clr_transform(cpm_df, pseudocount=0.5):
    """
    Centered Log-Ratio (CLR) transformation.

    Rationale:
      Microbiome sequencing data are inherently compositional — each sample's
      total abundance is constrained by library size, so ASVs are not
      independently variable. This induces a compositional negative bias in
      Spearman correlations: even if two ASVs are truly independent, their
      relative abundances can appear negatively correlated due to the
      "sum-to-constant" constraint.

      CLR maps data from simplex space to real space:
        CLR(x_i) = ln(x_i / g(x))
      where g(x) = (∏_{j=1}^p x_j)^{1/p} is the geometric mean of the sample.

      Transformed data are freed from the sum-to-constant constraint, making
      subsequent correlation analysis more faithful to true covariation.

    Notes:
      - Zeros cannot be log-transformed, so a pseudocount is added first.
      - CLR affects Spearman rank correlation less than Pearson (rank invariance),
        but still helps reduce compositional bias.
      - Most effective on deeply filtered ASV matrices (most samples have values).
    """
    if pseudocount > 0:
        log(f"CLR transform: adding pseudocount {pseudocount} for zero handling ...")
        x = cpm_df + pseudocount
    else:
        x = cpm_df.copy()

    # Per-sample geometric mean: exp(mean(ln(x)))
    geom_mean = np.exp(np.log(x).mean(axis=1))

    # CLR: ln(x / g(x))
    clr = np.log(x.div(geom_mean, axis=0))
    log("CLR transform complete.")
    return clr


# --------------------------------------------------------------
#  4. Spearman correlation matrix (vectorized)
# --------------------------------------------------------------
def compute_spearman_matrix(data):
    """
    Vectorized computation of all pairwise Spearman rank correlations.

    Rationale:
      An ASV-ASV co-occurrence network requires p×(p-1)/2 correlation coefficients.
      Pairwise loops are extremely inefficient; we use matrix operations:
        Let X be the n×p rank matrix of ASVs (n samples, p ASVs).

        1. Rank transform: sort each column independently and replace with ranks (1..n)
        2. Center: X_c = X - mean(X)                       (n×p)
        3. Covariance: Cov = X_cᵀ @ X_c / (n-1)            (p×p)
        4. Std deviation: σ_i = sqrt(ΣX_c² / (n-1))        (p-vector)
        5. Correlation matrix: ρ = D⁻¹ Cov D⁻¹,
           where D = diag(σ₀, σ₁, ..., σ_{p-1})            (p×p)
           Equivalent element-wise: ρ_ij = Cov_ij / (σ_i × σ_j)

        Since ρ is symmetric with unit diagonal (ASV with itself),
        only the upper triangle (i < j) is needed.

    Parameters:
      data: DataFrame, rows=samples, columns=ASVs, values=abundance

    Returns:
      rho_full : np.ndarray (p×p), full correlation matrix
      triu_i   : np.ndarray, upper-triangle row indices
      triu_j   : np.ndarray, upper-triangle column indices
      rho_vec  : np.ndarray, upper-triangle ρ values
      p_vals   : np.ndarray, p-values for each upper-triangle pair
      asv_names: list, ASV names (matching matrix column order)
    """
    n, p = data.shape
    total_pairs = p * (p - 1) // 2
    log(f"Vectorized Spearman: {p} ASVs → {total_pairs} pairs (n={n}) ...")

    # Rank transform: sort per column → replace with ranks
    ranks = data.rank(axis=0).values                    # n×p

    # Center
    centered = ranks - ranks.mean(axis=0)                # n×p

    # Covariance matrix (p×p)
    cov = centered.T @ centered / (n - 1)                # p×p

    # Standard deviation vector (p,)
    sd = np.sqrt(np.sum(centered ** 2, axis=0) / (n - 1))

    # ρ matrix: ρ_ij = Cov_ij / (σ_i × σ_j)
    rho_full = cov / np.outer(sd, sd)                    # p×p
    rho_full = np.clip(rho_full, -1.0, 1.0)              # numerical stability

    # Extract upper-triangle indices (excluding diagonal)
    triu_i, triu_j = np.triu_indices(p, k=1)             # both (total_pairs,) arrays

    # Upper-triangle ρ values
    rho_vec = rho_full[triu_i, triu_j]                  # (total_pairs,)

    # t-statistic → two-tailed p-value
    t_stat = rho_vec * np.sqrt(
        (n - 2) / np.maximum(1 - rho_vec ** 2, 1e-300)
    )
    p_vals = 2 * t_dist.sf(np.abs(t_stat), df=n - 2)
    p_vals = np.maximum(p_vals, 1e-300)

    asv_names = data.columns.tolist()
    log(f"  Full matrix: {p}×{p}, upper triangle: {total_pairs} pairs")
    return rho_full, triu_i, triu_j, rho_vec, p_vals, asv_names


# --------------------------------------------------------------
#  5. Global Benjamini–Hochberg FDR correction
# --------------------------------------------------------------
def bh_fdr(pvals):
    """
    Global Benjamini–Hochberg False Discovery Rate (FDR) correction.

    Unlike the ASV-ARG version: ASV-ASV pairs have no natural grouping
    by gene, so we apply a single global correction across all upper-triangle
    p-values.
    """
    m = len(pvals)
    if m == 0:
        return np.array([], dtype=float)

    order = np.argsort(pvals)
    sorted_p = pvals[order]

    # Compute critical values at each rank: q_i = p_i × m / (i+1)
    ratios = m / np.arange(1, m + 1, dtype=float)
    q_sorted = sorted_p * ratios

    # Monotonicity constraint: cumulative minimum from the tail
    # q_i = min_{k >= i} (p_k × m / (k+1))
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)

    # Restore original order
    qvals = np.zeros(m, dtype=float)
    qvals[order] = q_sorted
    return qvals


# --------------------------------------------------------------
#  6. Export
# --------------------------------------------------------------
def export_edges(rho_full, triu_i, triu_j, rho_vec, p_vals, q_vals,
                 asv_names, sig_mask, stability_arr=None,
                 out_path=OUTPUT_EDGES, stability_cutoff=0.80):
    """
    Export edge list CSV (Gephi-compatible format + optional full version).

    Three-tier filtering for Gephi network:
      (a) FDR-significant (q < Q_THRESHOLD)
      (b) LOO-stable (stability ≥ stability_cutoff)  [if LOO enabled]
      (c) Strong correlation (|ρ| ≥ MIN_ABS_RHO)

    Key difference from the ASV-ARG version:
      - Source and Target are both ASVs (same node type)
      - No ARG-related output
    """
    p = len(asv_names)
    has_loo = stability_arr is not None

    # --- Build final filter mask ---
    if has_loo:
        final_mask = sig_mask & (stability_arr >= stability_cutoff) & (np.abs(rho_vec) >= MIN_ABS_RHO)
    else:
        final_mask = sig_mask & (np.abs(rho_vec) >= MIN_ABS_RHO)

    n_sig = sig_mask.sum()
    n_final = final_mask.sum()

    # --- Degree filter: remove low-connectivity nodes ---
    final_i = triu_i[final_mask]
    final_j = triu_j[final_mask]
    final_rho = rho_vec[final_mask]
    final_q = q_vals[final_mask]

    # Build temporary DataFrame for degree calculation
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

    # Count degree per node (undirected, Source+Target combined)
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
        log(f"  Degree filter (degree≥{MIN_DEGREE}): removed {n_removed_nodes} low-degree nodes, "
            f"{n_removed_edges} associated edges")

    # --- Top-N core subnetwork: keep highest-degree nodes ---
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
        log(f"  Top-{TOP_N_NODES} core subnetwork: retained {len(top_nodes)} highest-degree nodes, "
            f"{len(df)} edges (removed {n_removed_top_nodes} nodes, "
            f"{n_removed_top_edges} edges)")

    # --- Gephi format output ---
    gephi_path = out_path.replace('.csv', '_gephi.csv')
    df.to_csv(gephi_path, index=False)
    log(f"Written: {gephi_path}  ({len(df)} rows — Gephi-compatible)")

    # --- [Optional] Full results output ---
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
        log(f"Written: {out_path}  ({len(full_df)} rows — full results)")

    # --- Summary statistics ---
    final_rho_arr = df['Spearman_rho'].values  # ρ values after degree filter
    n_final = len(df)
    n_pos = int((final_rho_arr > 0).sum())
    n_neg = n_final - n_pos
    asv_in_network = len(set(df['Source']).union(df['Target']))
    avg_degree = 2 * n_final / asv_in_network if asv_in_network > 0 else 0

    print()
    print("=" * 56)
    print("  A N A L Y S I S   S U M M A R Y")
    print("=" * 56)
    print(f"  ASVs tested:              {p}")
    print(f"  Total ASV-ASV pairs:      {len(rho_vec)}")
    print(f"  FDR-significant (q < {Q_THRESHOLD}):       {n_sig}")
    if has_loo:
        n_stable = int((sig_mask & (stability_arr >= stability_cutoff)).sum())
        print(f"  + LOO-stable (≥{stability_cutoff:.0%}):   {n_stable}")
    print(f"  ─────────────────────────────────────")
    print(f"  Gephi network (|ρ|≥{MIN_ABS_RHO}, degree≥{MIN_DEGREE}):")
    print(f"    Nodes: {asv_in_network}")
    print(f"    Edges: {n_final}")
    print(f"    Positive: {n_pos} ({100*n_pos/n_final:.1f}%)")
    print(f"    Negative: {n_neg} ({100*n_neg/n_final:.1f}%)")
    print(f"    Avg degree: {avg_degree:.2f}")
    print("=" * 56)

    return df


# --------------------------------------------------------------
#  7. Leave-One-Out stability test (matrix version)
# --------------------------------------------------------------
def loo_stability(data, triu_i, triu_j, rho_full, rho_vec, p_vals, q_vals,
                  asv_names):
    """
    Leave-One-Out stability test (matrix version).

    For each edge that passed FDR, hold out each sample in turn, fully
    recompute the Spearman matrix + FDR, and record the fraction of LOO
    iterations where the edge remains significant with the same sign.

    Differences from the ASV-ARG version:
      - Each iteration recomputes the full p×p correlation matrix, then extracts
        the upper triangle
      - Uses precomputed triu_i/triu_j to avoid repeated indexing
      - Uses global FDR instead of per-group FDR
    """
    n = data.shape[0]
    p = len(asv_names)
    sig_mask = q_vals < Q_THRESHOLD
    n_sig = sig_mask.sum()

    if n_sig == 0:
        log("No significant edges — skipping LOO.")
        return np.zeros(len(rho_vec), dtype=float)

    log(f"Leave-One-Out stability test: {n} iterations × {n_sig} significant edges ...")

    # Record significant edge: ASV names × sign → index
    sig_i_names = [asv_names[i] for i in triu_i[sig_mask]]
    sig_j_names = [asv_names[j] for j in triu_j[sig_mask]]
    sig_signs   = np.sign(rho_vec[sig_mask]).astype(int)

    edge2idx = {}
    for k in range(n_sig):
        edge2idx[(sig_i_names[k], sig_j_names[k], sig_signs[k])] = k

    survival = np.zeros(n_sig, dtype=float)
    all_idx = np.arange(n, dtype=int)
    data_vals = data.values  # pre-fetch numpy array for speed

    for loo in range(n):
        # Hold out the loo-th sample
        keep = np.delete(all_idx, loo)

        # Recompute Spearman matrix (subset version)
        sub = data_vals[keep]                # (n-1) × p
        sub_n = n - 1

        ranks = np.argsort(np.argsort(sub, axis=0), axis=0) + 1  # (n-1)×p
        ranks = ranks.astype(float)

        centered = ranks - ranks.mean(axis=0)
        cov = centered.T @ centered / (sub_n - 1)
        sd = np.sqrt(np.sum(centered ** 2, axis=0) / (sub_n - 1))
        loo_rho = cov / np.outer(sd, sd)
        loo_rho = np.clip(loo_rho, -1.0, 1.0)

        # Upper-triangle ρ
        loo_rho_vec = loo_rho[triu_i, triu_j]

        # p-values
        t = loo_rho_vec * np.sqrt(
            (sub_n - 2) / np.maximum(1 - loo_rho_vec ** 2, 1e-300)
        )
        loo_p = 2 * t_dist.sf(np.abs(t), df=sub_n - 2)
        loo_p = np.maximum(loo_p, 1e-300)

        # Global FDR
        loo_q = np.zeros_like(loo_p)
        order = np.argsort(loo_p)
        ratios = len(loo_p) / np.arange(1, len(loo_p) + 1, dtype=float)
        qs = loo_p[order] * ratios
        qs = np.minimum.accumulate(qs[::-1])[::-1]
        loo_q[order] = np.clip(qs, 0.0, 1.0)

        # Significant edges in this iteration
        loo_sig = loo_q < Q_THRESHOLD
        loo_sig_i_names = [asv_names[i] for i in triu_i[loo_sig]]
        loo_sig_j_names = [asv_names[j] for j in triu_j[loo_sig]]
        loo_sig_signs   = np.sign(loo_rho_vec[loo_sig]).astype(int)

        # Match against baseline significant edges
        for k in range(len(loo_sig_i_names)):
            key = (loo_sig_i_names[k], loo_sig_j_names[k], loo_sig_signs[k])
            idx = edge2idx.get(key)
            if idx is not None:
                survival[idx] += 1.0

        if (loo + 1) % 5 == 0 or loo == n - 1:
            log(f"  LOO: {loo + 1}/{n}")

    # Stability vector (corresponding to positions of significant edges in ρ_vec)
    stability = survival / n
    n_stable = int((stability >= 0.80).sum())
    log(f"  stable (≥80%): {n_stable}/{n_sig} edges")

    # Build full stability array (for all upper-triangle pairs)
    full_stability = np.zeros(len(rho_vec), dtype=float)
    sig_positions = np.where(sig_mask)[0]
    for k, pos in enumerate(sig_positions):
        full_stability[pos] = stability[k]

    return full_stability


# --------------------------------------------------------------
#  7.  Node taxonomy annotation
# --------------------------------------------------------------
def read_taxonomy():
    """Read taxonomy.txt, return {ASV: (Phylum, Genus)} dict."""
    log(f"Reading {os.path.basename(TAXONOMY_FILE)} ...")
    tax_df = pd.read_csv(TAXONOMY_FILE, sep='\t', header=0)
    tax_df.columns = tax_df.columns.str.strip()
    result = {}
    for _, row in tax_df.iterrows():
        phylum = row['Phylum'] if pd.notna(row['Phylum']) else 'Unassigned'
        genus  = row['Genus']  if pd.notna(row['Genus'])  else 'Unassigned'
        result[row['OTUID']] = (phylum, genus)
    log(f"  Loaded taxonomy for {len(result)} ASVs")
    return result


def export_node_taxonomy(gephi_edges, tax_map, out_path):
    """Annotate every node in the Gephi edge list with Phylum / Genus, write to txt."""
    nodes = set(gephi_edges['Source']).union(set(gephi_edges['Target']))
    rows = []
    for node in sorted(nodes):
        phylum, genus = tax_map.get(node, ('Unassigned', 'Unassigned'))
        rows.append({'Node': node, 'Phylum': phylum, 'Genus': genus})
    result = pd.DataFrame(rows)
    result.to_csv(out_path, sep='\t', index=False)
    log(f"Written: {out_path}  ({len(result)} nodes)")
    return result


# --------------------------------------------------------------
#  main
# --------------------------------------------------------------
def main():
    print("=" * 56)
    print("  ASV Co-occurrence Network Analysis")
    print("  Microbial co-occurrence pattern mining")
    print("=" * 56)

    # Ensure output directory exists
    os.makedirs(os.path.join(PROJECT_ROOT, "results"), exist_ok=True)

    # 1 — Read feature table
    counts, otu_ids = read_asv_table(INPUT_FILE)
    print(f"  Raw ASVs: {len(otu_ids)}")
    print(f"  Samples:  {counts.shape[0]}")

    # 2 — Filter + CPM normalization
    cpm, _ = filter_asvs(counts, MIN_OCCURRENCE, MIN_TOTAL_CPM)

    # 3 — CLR transform (compositional bias correction, always on)
    analysis_data = clr_transform(cpm, PSEUDOCOUNT)
    log("Using CLR-transformed data for Spearman correlation analysis.")

    # 4 — Spearman correlation matrix
    rho_full, triu_i, triu_j, rho_vec, p_vals, asv_names = \
        compute_spearman_matrix(analysis_data)

    # 5 — Global FDR
    log("Applying global BH FDR correction ...")
    q_vals = bh_fdr(p_vals)
    n_sig = int((q_vals < Q_THRESHOLD).sum())
    log(f"  Significant: {n_sig}/{len(p_vals)} (q<{Q_THRESHOLD})")

    # 6 — LOO stability
    stability = None
    if DO_LOO:
        stability = loo_stability(
            analysis_data,
            triu_i, triu_j,
            rho_full, rho_vec, p_vals, q_vals,
            asv_names
        )

    # 7 — Export
    gephi_df = export_edges(rho_full, triu_i, triu_j, rho_vec, p_vals, q_vals,
                            asv_names, q_vals < Q_THRESHOLD, stability,
                            OUTPUT_EDGES, LOO_STABILITY)

    # 8 — Node taxonomy annotation
    tax_map = read_taxonomy()
    export_node_taxonomy(gephi_df, tax_map, NODE_TAXONOMY_OUT)

    print("\nDone. Import the Gephi CSV into Gephi (File → Import Spreadsheet → Undirected) for visualization.")


if __name__ == '__main__':
    main()
