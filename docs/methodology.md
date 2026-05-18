# Methodology: ASV Co-occurrence Network Analysis

## Overview

This pipeline constructs a **co-occurrence network** among microbial ASVs from amplicon sequencing data. Co-occurrence networks are a fundamental tool in microbial ecology for inferring potential ecological interactions (symbiosis, competition, niche overlap) from abundance patterns across samples.

Two ASVs that show statistically significant abundance correlations across samples may share ecological niches, have obligate metabolic dependencies, or respond similarly to environmental gradients.

## Step-by-Step Method

### 1. ASV Filtering

Low-abundance and rare ASVs are removed:

- **Occurrence filter**: an ASV must be present (count > 0) in at least `MIN_OCCURRENCE` samples
- **Abundance filter**: an ASV must have a total CPM ≥ `MIN_TOTAL_CPM` across all samples

### 2. CPM Normalization

Raw read counts are converted to Counts Per Million (CPM):

$$\text{CPM}_{ij} = \frac{\text{count}_{ij}}{\sum_k \text{count}_{kj}} \times 10^6$$

### 3. CLR Transformation

Microbiome data are **compositional** — each sample's feature vector sums to its total read count, not to a biologically meaningful total. This creates a "sum-to-constant" constraint that induces artificial negative correlations between features (compositional bias).

The Centered Log-Ratio (CLR) transformation mitigates this by mapping data from simplex space to real space:

$$\text{CLR}(x_i) = \ln\left(\frac{x_i}{g(x)}\right)$$

where $g(x) = \left(\prod_{j=1}^p x_j\right)^{1/p}$ is the geometric mean of all features in a sample.

A pseudocount (default: 0.5) is added to all values before transformation to handle zeros.

### 4. Spearman Correlation Matrix

Unlike the ASV–ARG pipeline (which computes ASV × ARG pairs), the co-occurrence network must compute ASV × ASV pairs — a p×p matrix of all pairwise correlations. This is computationally intensive, so a fully vectorized approach is used:

1. Rank-transform the n×p data matrix
2. Compute the p×p covariance matrix: `Cov = X_cᵀ @ X_c / (n-1)`
3. Convert to correlation matrix: `ρ_ij = Cov_ij / (σ_i × σ_j)`
4. Extract the upper triangle (i < j) to avoid redundancy and self-correlations

### 5. Global FDR Correction

Because there is no natural grouping for ASV–ASV pairs, we apply a **single global** Benjamini–Hochberg correction across all upper-triangular pairs.

### 6. Leave-One-Out (LOO) Stability Test

Same logic as the ASV-ARG pipeline: for each of N samples, hold it out, recompute the full correlation matrix and FDR, and record which edges remain significant with the same sign. Edges must survive in ≥ `LOO_STABILITY` fraction of iterations.

### 7. Post-Processing Filters

**Degree filter**: nodes with degree < `MIN_DEGREE` are removed. These "dangling" nodes connect to only 1–2 others and contribute little to network structure.

**Top-N core subnetwork**: the top `TOP_N_NODES` highest-degree nodes are retained for visualization, producing a legible "core subnetwork" — a standard practice in microbial network publications.

### 8. Network Export

Edges must satisfy: FDR-significant AND LOO-stable AND strong correlation (|ρ| ≥ MIN_ABS_RHO).

## References

- Aitchison, J. (1982). The statistical analysis of compositional data. *Journal of the Royal Statistical Society: Series B*, 44(2), 139–160.
- Faust, K., & Raes, J. (2012). Microbial interactions: from networks to models. *Nature Reviews Microbiology*, 10(8), 538–550.
- Friedman, J., & Alm, E. J. (2012). Inferring correlation networks from genomic survey data. *PLoS Computational Biology*, 8(9), e1002687.
- Love, M. I., et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.
