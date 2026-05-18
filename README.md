# ASV Co-occurrence Network Analysis

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

**中文** | [English](#english)

<a name="english"></a>

## English

### Overview

This pipeline builds a **co-occurrence network** between Amplicon Sequence Variants (ASVs) from microbial community amplicon sequencing data. It computes all pairwise Spearman correlations, applies multiple-testing correction and stability testing, and exports edge lists ready for network visualization.

The key methodological steps are:

1. **ASV filtering** — occurrence frequency and total abundance filters
2. **CPM normalization** — Counts Per Million to correct for sequencing depth
3. **CLR transformation** — Centered Log-Ratio transform to mitigate compositional bias
4. **Spearman correlation matrix** — vectorized computation of all ASV–ASV pairs
5. **Global FDR correction** — Benjamini–Hochberg across all pairwise tests
6. **Leave-One-Out (LOO) stability test** — full recomputation with each sample held out
7. **Degree & Top-N filtering** — removes low-degree nodes and keeps the core subnetwork
8. **Gephi export** — edge lists ready for network visualization in [Gephi](https://gephi.org/)

### Project Structure

```
asv-cooccurrence-network/
├── data/                       # Input data (gitignored - add your own)
│   ├── README.md               # Data format specifications
│   ├── otutab.txt              # ASV/OTU abundance table
│   └── taxonomy.txt            # ASV taxonomy annotation
├── docs/
│   ├── methodology.md          # Methodology (English)
│   └── methodology.zh.md       # Methodology (中文)
├── results/                    # Output results (gitignored)
│   └── .gitkeep
├── src/
│   ├── asv_cooccurrence_network_analysis.py    # 中文版
│   └── asv_cooccurrence_network_analysis_en.py # English version
├── .gitignore
├── LICENSE
├── README.md
└── requirements.txt
```

### Input Files

| File | Format | Description |
|------|--------|-------------|
| `otutab.txt` | TSV | OTUID \| Sample1 \| Sample2 \| ... |
| `taxonomy.txt` | TSV | OTUID \| Phylum \| Genus \| ... |

See [data/README.md](data/README.md) for detailed format specifications.

### Output Files

| File | Description |
|------|-------------|
| `results/asv_cooccurrence_network_edges_gephi.csv` | Significant + stable edges in Gephi-compatible format (core subnetwork) |
| `results/asv_cooccurrence_network_node_taxonomy.txt` | Node annotation table with Phylum / Genus for each node |

> Note: The full edge table (all ASV-ASV pairs) can be enabled via `OUTPUT_FULL_CSV = True` in the configuration, but is disabled by default as it can exceed 7M+ rows.

### Quick Start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Prepare input files in the `data/` directory:
#    data/otutab.txt
#    data/taxonomy.txt

# 3. Run the analysis
cd src
python asv_cooccurrence_network_analysis.py

# 4. Find results in results/
#    Import *gephi.csv into Gephi (File → Import Spreadsheet → Undirected)
```

### Configuration

All tunable parameters are at the top of `src/asv_cooccurrence_network_analysis.py`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MIN_OCCURRENCE` | 14 | Minimum samples an ASV must appear in |
| `MIN_TOTAL_CPM` | 200 | Minimum total CPM across all samples |
| `USE_CLR` | True | Enable Centered Log-Ratio transformation |
| `PSEUDOCOUNT` | 0.5 | Pseudocount for zero-value handling in CLR |
| `Q_THRESHOLD` | 0.05 | FDR q-value significance cutoff |
| `DO_LOO` | True | Enable Leave-One-Out stability test |
| `LOO_STABILITY` | 0.90 | Minimum fraction of LOO iterations to retain edge |
| `MIN_ABS_RHO` | 0.70 | Minimum \|Spearman ρ\| for Gephi output |
| `MIN_DEGREE` | 3 | Minimum node degree to retain in network |
| `TOP_N_NODES` | 30 | Keep only top-N highest-degree nodes for Gephi visualization |
| `OUTPUT_FULL_CSV` | False | Export the complete ASV-ASV edge table |

### Visualization in Gephi

1. Open Gephi → File → Import Spreadsheet
2. Select `asv_cooccurrence_network_edges_gephi.csv`
3. Choose **Undirected** graph type
4. Color nodes by Phylum (use the exported taxonomy file)
5. Size nodes by degree

### Methodological Notes

#### CLR Transformation

Microbial amplicon data are compositional — each sample's total count is constrained by sequencing depth, which introduces artificial negative correlations between ASVs. The Centered Log-Ratio transformation maps data from simplex space to real space:

```
CLR(x_i) = ln(x_i / g(x))
```

where `g(x)` is the geometric mean of all features in a sample. This alleviates the compositional bias that plagues correlation-based network inference.

#### Top-N Subnetwork

Raw co-occurrence networks can contain thousands of nodes, which are impractical to visualize meaningfully. Following standard practice in microbial network papers, we extract the top-N most connected nodes to visualize the **core subnetwork** — the most ecologically central microbial players.

### Citation

If you use this pipeline in your research, please cite the repository:

```bibtex
@software{asv_cooccurrence_network,
  author = Xianyu Yin,
  title = {ASV Co-occurrence Network Analysis},
  year = {2026},
  url = {https://github.com/YinXianYu-bioinfor/asv-cooccurrence-network}
}
```

### License

MIT License — see [LICENSE](LICENSE).


---

## 中文文档

### 概述

本流程用于构建微生物群落 Amplicon Sequence Variants（ASV）之间的**共现网络（co-occurrence network）**。通过计算所有 ASV 两两之间的 Spearman 相关性，进行多重检验校正和稳定性检验，导出可直接用于网络可视化的边列表。

共现网络是微生物生态学中的基础分析工具，用于从跨样本的丰度变化模式中推断潜在的生态互作关系（共生、竞争、生态位重叠等）。

主要方法步骤：

1. **ASV 过滤** — 基于出现频率和总丰度筛选
2. **CPM 标准化** — Counts Per Million 校正测序深度
3. **CLR 变换** — Centered Log-Ratio 变换减轻成分数据偏差
4. **Spearman 相关矩阵** — 向量化计算所有 ASV–ASV 对
5. **全局 FDR 校正** — 对所有两两检验执行 Benjamini–Hochberg 校正
6. **留一法（LOO）稳定性检验** — 每次剔除一个样本后完整重算
7. **度过滤 & Top-N 子网络** — 去除低连接节点，保留核心子网络
8. **Gephi 导出** — 生成可直接导入 [Gephi](https://gephi.org/) 的边列表

详细方法原理见 [docs/methodology.zh.md](docs/methodology.zh.md)。

### 项目结构

```
asv-cooccurrence-network/
├── data/                       # 输入数据（gitignored，请自行添加）
│   ├── README.md               # 数据格式说明
│   ├── otutab.txt              # ASV/OTU 丰度表
│   └── taxonomy.txt            # ASV 分类注释
├── docs/
│   ├── methodology.md          # 方法论（英文）
│   └── methodology.zh.md       # 方法论（中文）
├── results/                    # 输出结果（gitignored）
│   └── .gitkeep
├── src/
│   ├── asv_cooccurrence_network_analysis.py    # 中文版
│   └── asv_cooccurrence_network_analysis_en.py # English version
├── .gitignore
├── LICENSE
├── README.md
└── requirements.txt
```

### 输入文件

| 文件 | 格式 | 说明 |
|------|------|------|
| `otutab.txt` | TSV | OTUID \| Sample1 \| Sample2 \| ... |
| `taxonomy.txt` | TSV | OTUID \| Phylum \| Genus \| ... |

详细格式说明见 [data/README.md](data/README.md)。

### 输出文件

| 文件 | 说明 |
|------|------|
| `results/asv_cooccurrence_network_edges_gephi.csv` | 显著 + 稳定的边（核心子网络），Gephi 可直接导入 |
| `results/asv_cooccurrence_network_node_taxonomy.txt` | 节点注释表（Phylum / Genus） |

> 注：完整 ASV–ASV 边表可通过 `OUTPUT_FULL_CSV = True` 开启（默认关闭，可能超过 700 万行）。

### 快速开始

```bash
# 1. 安装依赖
pip install -r requirements.txt

# 2. 将输入文件放入 data/ 目录：
#    data/otutab.txt
#    data/taxonomy.txt

# 3. 运行分析
cd src
python asv_cooccurrence_network_analysis.py

# 4. 结果在 results/ 目录下
#    将 *gephi.csv 导入 Gephi（File → Import Spreadsheet → Undirected）
```

### 参数配置

所有可调参数位于 `src/asv_cooccurrence_network_analysis.py` 顶部：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `MIN_OCCURRENCE` | 14 | ASV 至少出现的样本数 |
| `MIN_TOTAL_CPM` | 200 | ASV 在所有样本中的最小总 CPM |
| `USE_CLR` | True | 是否启用 Centered Log-Ratio 变换 |
| `PSEUDOCOUNT` | 0.5 | CLR 中零值处理的伪计数 |
| `Q_THRESHOLD` | 0.05 | FDR q 值显著性阈值 |
| `DO_LOO` | True | 是否执行留一法稳定性检验 |
| `LOO_STABILITY` | 0.90 | 边在 LOO 迭代中需保留的最小比例 |
| `MIN_ABS_RHO` | 0.70 | Gephi 输出所需的最小 \|Spearman ρ\| |
| `MIN_DEGREE` | 3 | 节点保留的最小度（连接数） |
| `TOP_N_NODES` | 30 | Gephi 可视化仅保留前 N 个高度节点 |
| `OUTPUT_FULL_CSV` | False | 是否导出完整 ASV–ASV 边表 |

### Gephi 可视化

1. 打开 Gephi → File → Import Spreadsheet
2. 选择 `asv_cooccurrence_network_edges_gephi.csv`
3. 图类型选择 **Undirected**
4. 按 Phylum 对节点着色（使用导出的 taxonomy 注释文件）
5. 按度（degree）调整节点大小

### 方法说明

#### CLR 变换

微生物扩增子数据具有成分数据特性——每个样本的总读数受测序深度约束，这会在 ASV 之间引入人为的负相关性。Centered Log-Ratio 变换将数据从单纯形空间映射到实数空间：

```
CLR(x_i) = ln(x_i / g(x))
```

其中 `g(x)` 是一个样本中所有特征的几何平均数。该变换可减轻成分数据偏差对相关网络推断的影响。

#### Top-N 核心子网络

原始共现网络可能包含数千个节点，难以有效可视化。遵循微生物网络论文的标准做法，我们提取前 N 个连接最多的节点来可视化**核心子网络**——即生态学上最核心的微生物类群。

### 引用

如果您在研究中使用了本流程，请引用此仓库：

```bibtex
@software{asv_cooccurrence_network,
  author = Xianyu Yin,
  title = {ASV Co-occurrence Network Analysis},
  year = {2026},
  url = {https://github.com/YinXianYu-bioinfor/asv-cooccurrence-network}
}
```

### 许可证

MIT License — 详见 [LICENSE](LICENSE)。

---
