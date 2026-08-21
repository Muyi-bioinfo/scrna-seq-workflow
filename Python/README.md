# Python/ — scanpy 平行版本（规划中）

> 目录结构借鉴 [hossainlab/sc-workflow](https://github.com/hossainlab/sc-workflow)（R/ + Python/ 双语言分层）。
> 上游 cellranger 定量是语言无关的，同一份 FASTQ/矩阵可同时喂给 R 和 Python 两条管线。

## 规划路线（与 R/ 各阶段一一对应）

| 计划脚本 | 对应 R/ | 内容 | 核心库 |
|---------|---------|------|--------|
| `01_load_qc.py` | `R/01_load_qc.R` | 读 h5ad/矩阵 + 质控（线粒体/基因数/doublet） | scanpy, scrublet |
| `02_preprocess_cluster.py` | `R/02_preprocess_cluster.R` | 标准化/HVG/PCA + 批次校正 + UMAP/Leiden | scanpy, harmonypy/scanorama |
| `03_annotate.py` | `R/03_annotate.R` | marker 鉴定 + 注释 | scanpy rank_genes_groups, celltypist |
| `04_multi_group_plot.py` | `R/04_multi_group_plot.R` | 多组可视化 | scanpy, matplotlib |
| `05_diff_gsea.py` | `R/05_diff_gsea.R` | 差异 + 富集（pseudobulk） | scanpy, decoupler |
| `06_velocity.py` | — | RNA 速率 | scvelo（另有独立项目已跑通） |

## 已有 Python 单细胞经验（交叉引用）

- 环境 `scvelo_env` + 独立拟时序项目：velocyto/scvelo 已真实跑通
- 设计原则与 R/ 保持一致：分阶段脚本 + config 集中参数 + batch 归档

## 状态

- [ ] 脚本未开始编写（R 版本为主线；此目录为路线图，面试可讲规划）
- [ ] 参考 sc-workflow 的 `Python/00_*.py ~ 08_*.py` 对齐实现
