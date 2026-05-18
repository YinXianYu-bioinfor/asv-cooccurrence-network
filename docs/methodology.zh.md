# 方法论：ASV 共现网络分析

## 概述

本流程用于从扩增子测序数据中构建微生物 ASV 之间的**共现网络（co-occurrence network）**。共现网络是微生物生态学中的基础分析工具，用于从跨样本的丰度变化模式中推断潜在的生态互作关系（共生、竞争、生态位重叠等）。

如果两个 ASV 在多个样本中表现出统计显著的丰度相关性，它们可能共享相似的生态位、存在代谢依赖关系，或对环境梯度有相似的响应。

## 逐步方法说明

### 1. ASV 过滤

去除低丰度和稀有的 ASV：

- **出现频率过滤**：ASV 必须在至少 `MIN_OCCURRENCE` 个样本中出现（count > 0）
- **丰度过滤**：ASV 在所有样本中的总 CPM 必须 ≥ `MIN_TOTAL_CPM`

### 2. CPM 标准化

将原始读数转换为 Counts Per Million（CPM），以校正样本间测序深度差异：

$$\text{CPM}_{ij} = \frac{\text{count}_{ij}}{\sum_k \text{count}_{kj}} \times 10^6$$

### 3. CLR 变换

微生物组数据具有**成分数据（compositional）**特性——每个样本的特征向量之和为该样本的总读数，而非具有生物学意义的总和。这种"总和为常数"的约束会在特征之间引入人为的负相关性（成分偏差）。

Centered Log-Ratio（CLR）变换通过将数据从单纯形空间映射到实数空间来缓解这一问题：

$$\text{CLR}(x_i) = \ln\left(\frac{x_i}{g(x)}\right)$$

其中 $g(x) = \left(\prod_{j=1}^p x_j\right)^{1/p}$ 是一个样本中所有特征的几何平均数。

变换前会向所有值添加伪计数（pseudocount，默认 0.5）以处理零值。

### 4. Spearman 相关矩阵

与 ASV–ARG 流程（计算 ASV × ARG 对）不同，共现网络需要计算 ASV × ASV 对——即所有两两相关性的 p×p 矩阵。这在计算上较为密集，因此采用完全向量化的方法：

1. 对 n×p 数据矩阵进行秩变换
2. 计算 p×p 协方差矩阵：`Cov = X_cᵀ @ X_c / (n-1)`
3. 转换为相关矩阵：`ρ_ij = Cov_ij / (σ_i × σ_j)`
4. 提取上三角元素（i < j）以避免冗余和自相关

### 5. 全局 FDR 校正

由于 ASV–ASV 对之间没有天然的分组依据，我们在所有上三角对之间应用**单一全局** Benjamini–Hochberg 校正。

### 6. 留一法（LOO）稳定性检验

与 ASV-ARG 流程逻辑相同：对 N 个样本中的每一个，将其剔除后重新计算完整的相关矩阵和 FDR，并记录哪些边仍然显著且符号一致。一条边必须在 ≥ `LOO_STABILITY` 比例的 LOO 迭代中存活。

### 7. 后处理过滤

**度过滤**：移除度 < `MIN_DEGREE` 的节点。这些"悬挂"节点仅连接 1–2 个其他节点，对网络结构的贡献有限。

**Top-N 核心子网络**：保留前 `TOP_N_NODES` 个最高度的节点用于可视化，生成可读的"核心子网络"——这是微生物网络论文中的标准做法。

### 8. 网络导出

边的筛选条件为：FDR 显著且 LOO 稳定且强相关（|ρ| ≥ MIN_ABS_RHO）。

## 参考文献

- Aitchison, J. (1982). The statistical analysis of compositional data. *Journal of the Royal Statistical Society: Series B*, 44(2), 139–160.
- Faust, K., & Raes, J. (2012). Microbial interactions: from networks to models. *Nature Reviews Microbiology*, 10(8), 538–550.
- Friedman, J., & Alm, E. J. (2012). Inferring correlation networks from genomic survey data. *PLoS Computational Biology*, 8(9), e1002687.
- Love, M. I., et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.
