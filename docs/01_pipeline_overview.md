# 整体流程文档：scRNA-seq 标准分析流程

> 从 10x 原始 FASTQ 到细胞注释、差异分析与 GSEA 的完整流程说明。
> 每个分析阶段对应 `R/` 中一个脚本，参数集中在 `config/`，可单独重跑。

## 架构设计（借鉴 GitHub 主流项目）

本项目结构调研了 GitHub 上两类主流 scRNA-seq 项目后设计：

| 设计 | 借鉴来源 | 本项目实现 |
|------|---------|-----------|
| 每步一个独立脚本 | hossainlab/sc-workflow（`R/01_qc.R` ~ `R/08_*.R`） | `R/01_load_qc.R` ~ `09_cell_communication.R` |
| 参数集中配置 | chunyuanzhang/scRNA-seq、Prezi atlas pipeline（`config/config.yaml`） | `config/config.yaml` + `config/config.multi.yaml`（双 profile） |
| 样本表 | chunyuanzhang 的 `sampleandpopulation.tsv` | `config/sample_sheet.csv`（sample, group, fastqs, matrix） |
| 顺序 runner + 断点续跑 | hossainlab 的 bash `&&` 串联（升级版） | `run_pipeline.sh`（.done 标记跳过已跑步骤） |
| 每步独立结果目录 | hossainlab 的 `results/01_qc/` | `output/<batch>/01_load_qc/` ... `output/<batch>/09_cell_communication/` |

> 工业级项目多用工作流引擎（[nf-core/scrnaseq](https://github.com/nf-core/scrnaseq) 用 Nextflow，[epigen/dea_seurat](https://github.com/epigen/dea_seurat) 用 Snakemake）。本流程保持无额外依赖，但每步一脚本 + 产物断点的设计使其可以平滑迁移到 Snakemake/Nextflow，也天然支持 HPC（每步 = 一个投递任务）。

## 流程总览

```
原始数据          上游定量            下游分析（分阶段脚本）
─────────        ──────────         ─────────────────────────────
FASTQ    ──► bash/00_run_cellranger.sh ──► filtered_feature_bc_matrix ──► 01_load_qc
(R1:BC+UMI)   (比对+UMI计数)       (barcodes/features/matrix)   (读入+质控)
(R2:cDNA)                                                          │
                                                                   ▼
                                              02_preprocess_cluster（预处理+聚类）
                                              single: 标准化/HVG/缩放/PCA → UMAP → 分群
                                              multi : 同样预处理 + 去批次(harmony/rpca/mnn) → 分群
                                                                   │
                                                                   ▼
                                                             03_annotate（80%时间）
                                                             （注释 + 可选亚群细分）
                                                                   │
                                             ┌──────────────┬──────┴───────┬──────┐
                                             ▼              ▼              ▼      ▼
                                    04_multi_group_plot  05_diff_gsea  06_pseudobulk_de  07+ 高级
                                    (多组可视化)         (探索性DE+GSEA) (确认性DE,样本级) (基因集打分等)
```

## 脚本粒度设计（参考 hossainlab/sc-workflow）

脚本按**分析阶段**打包，而非每个操作一个脚本（12 个操作 → 7 个脚本）：

| 打包理由 | 合并 |
|---------|------|
| 读入和质控总是连着做，分开无独立复跑价值 | 01_load_qc = 读入 + 质控 |
| 预处理到分群是一条流水线，单样本/多样本只是中间多一步去批次 | 02_preprocess_cluster = 预处理 + 去批次 + 聚类 |
| 亚群细分是注释的延伸（对注释结果的精细调整） | 03_annotate = 注释 + 亚群（config 开关） |

这样每个脚本跑完都有一个**完整可讲的阶段性结果**，同时脚本总量（9 个）与 sc-workflow（8 个）相当，不至于碎片化。

## 模块详解

每个脚本 = 一个阶段 + config 对应段落 + output 一个子目录：

| 阶段 | 脚本 | 配置段 | 输入 → 输出 | 产出图 |
|------|------|--------|------------|--------|
| 0 上游定量 |  `bash/00_run_cellranger.sh` | `config/cellranger.yaml`（工具级独立配置） | FASTQ → `output/<batch>/00_cellranger/<id>/outs/` | web_summary.html（质检报告） |
| 1+2 读入+质控+双细胞 | `01_load_qc.R` | `mode:`、`single:`/`multi:`、`qc:`（含 doublet_enable） | 10x 矩阵（单样本）或样本表矩阵/FASTQ（多样本）→ `seurat_qc.rds` | qc_violin_before |
| 3+4+5 预处理+聚类 | `02_preprocess_cluster.R` | `preprocess:`、`cluster:`、`integrate:` | qc → `seurat_clustered.rds` | variable_features、elbow、dimplot_clusters、（多样本）整合前后对比 |
| 6+7 注释+亚群 | `03_annotate.R` | `annotate:`、`subcluster:` | clustered → `seurat_annotated.rds`（+ 亚群 rds） | featureplot_markers、dimplot_annotated、（亚群图） |
| 8 多组可视化 | `04_multi_group_plot.R` | `multi_group:` | annotated → 仅图 | 比例图、分组 DotPlot、热图 |
| 9 差异+GSEA | `05_diff_gsea.R` | `diff_gsea:` | annotated → 差异 rds + 图 | compareCluster、GSEA dotplot |
| 9+ pseudobulk 确认性 DE（样本级） | `06_pseudobulk_de.R` | `pseudobulk:`（配对 design / min 守卫；sample_id 由 01 写入） | annotated → 逐类型 DESeq2 结果 + `summary_degs.csv` | 逐类型 volcano/MA/PCA、全局 pseudobulk PCA、05 vs 06 对比条形图 |
| A 基因集打分（高级） | `07_gene_set_score.R` | `advanced:` | annotated → `seurat_scored.rds` | 通路打分 FeaturePlot |
| B 拟时序（高级） | `08_trajectory_analysis.R` | `advanced$trajectory`（preset / reduction / root_celltype） | annotated → `trajectory_cds.rds` + `trajectory_resolved.yaml` 参数快照 | 伪时间轨迹图、细胞类型轨迹图、分组轨迹图 |
| C 细胞通讯（高级） | `09_cell_communication.R` | `advanced$cellchat` | annotated → 通讯网络结果 | 通讯网络图、通路气泡图、组间差异比较 |

### 模块 0：cellranger 定量（上游）

- `bash bash/00_run_cellranger.sh [主config] [cellranger配置]`，工具级参数（参考基因组/硬件资源/expect_cells）在独立的 `config/cellranger.yaml`，与下游 R 分析解耦、跨项目可复用；样本级信息（fastqs）在样本表
- 10x 双端测序：R1 = 细胞 barcode + UMI，R2 = cDNA 插入片段
- 产出 `outs/web_summary.html` 质检报告（浏览器打开）与 `filtered_feature_bc_matrix`（下游输入）
- 详细讲解见 [02_cellranger_guide.md](02_cellranger_guide.md)

### 模块 1：数据读入

- `Read10X()` 读取 10x 三件套（或 h5），返回稀疏矩阵（dgCMatrix）——**行是基因，列是细胞**
- `CreateSeuratObject()`：`min.cells`（基因至少在 N 个细胞表达）、`min.features`（细胞至少 N 个基因）——阈值在 `qc:` 段（`min_cells` / `nfeature_min`），与后续质控过滤同源
- 多样本场景每个样本加 `group` 列（05 整合按它去批次）与 `sample_id` 列（pseudobulk 的生物学重复）：默认 = 样本表每行；合并矩阵场景（每行是"多供者池"，如 GSE96583）用 `multi$cell_metadata` 按 (group, barcode) join 逐细胞元数据覆盖（供者在 GEO tsne.df 的 `ind` 列）

### 模块 2：质控（QC）

三个核心指标：**nFeature_RNA**（基因数：过低空液滴，过高双细胞）、**nCount_RNA**（UMI 数）、**percent.mt**（线粒体比例：高 = 凋亡细胞，人 `^MT-` 小鼠 `^mt-`）。阈值在 `qc:` 段调整，改完 `--only 01` 重跑即可。

**双细胞检测（scDblFinder）**：一个液滴包裹两个细胞形成"假中间态"细胞群，必须检测并过滤（`qc$doublet_enable` 开关）。多样本时**逐样本检测**（双细胞率随建库批次变化），再 merge。

### 模块 3：预处理

NormalizeData（UMI 缩放到 10000 取对数）→ FindVariableFeatures（vst 找 2000 高变基因）→ ScaleData（中心化，默认只缩放高变基因）→ RunPCA + ElbowPlot（拐点选维度）。

### 模块 4：聚类

RunUMAP → FindNeighbors（KNN → SNN 图）→ FindClusters（resolution 越大群越多）。分辨率在 config 改完 `--only 02` 重跑，无需重跑前面步骤。

### 模块 5：多样本整合

- Merge 多样本（**v5 坑**：merge 后 split layers，harmony 路径须 `JoinLayers()`）
- 去批次三选一（`integrate$method`）：
  - `harmony`（默认）：对 PCA 坐标线性校正，经典、快
  - `rpca` / `mnn`：Seurat v5 统一接口 `IntegrateLayers()`
- 检验：整合前后 UMAP 按 group 着色对比

### 模块 6：细胞注释（占 80% 分析时间）

1. FindAllMarkers（presto 自动加速）→ top10 表（**pct.1 预筛选**，阈值 `annotate$marker_pct`，v5 logFC 公式变化坑）
2. marker 可视化验证：marker 来自**注释预设库**（`config/annotation_presets/<物种>/`，按细胞类型组织 + tier1/tier2 分级 + 别名 + 出处；加载时自动校验物种、marker 在数据中的存在率、cluster_map 类型名一致性）
3. 三部曲：确认群数 → `RenameIdents`（映射在 `annotate$cluster_map`，人工填写）→ 保存 celltype 列

### 模块 7：亚群细分

`FindSubCluster` 对指定大群（config 的 `subcluster$cluster`）在 SNN 图上重新聚类，找内部异质性（如 T 细胞大群 → Naive CD4/CD8 T）。

### 模块 8：多组比较可视化

细胞比例堆叠图（table + ggplot 百分比堆积）→ 分组 DotPlot（`split.by = group`）→ top5 marker 热图（downsample 抽样）。

### 模块 9：差异分析 + GSEA

- **策略一（群水平）**：`celltype.stim` 组合身份 → FindAllMarkers → `compareCluster` 批量 KEGG 富集
- **策略二（刺激前后）**：同一细胞类型处理 vs 对照 → `FindMarkers(logfc.threshold = 0)` → geneList 三部曲（取 logFC → 命名 → 降序排序）→ `GSEA`（Hallmark/KEGG，`set.seed` 保证可复现）

### 模块 9+：pseudobulk 确认性差异分析（sample-aware）

- **为什么**：模块 9 以细胞为统计单元——同一供者的细胞高度相关，当成独立样本就是**伪重复**（pseudoreplication），p 值虚小、DEG 数虚高。GSE96583 本质是 Kang 2017 的 8 供者配对数据，这是组间差异分析的统计要害
- **怎么做**：原始 counts 按 供者×处理×细胞类型 聚合成 pseudobulk（`Matrix::rowsum`，**统计单元 = sample，不再是 cell**）→ 每个细胞类型独立建 DESeq2 负二项模型。同一供者出现在 ≥2 组时自动用配对 design `~ sample_id + condition`（供者作 blocking 因子吸收个体差异）；非配对数据自动退化为 `~ condition`；batch 列可用且不与 condition 混淆时插到最前
- **供者从哪来**：GSE96583 的合并矩阵不带 sample_id（orig.ident 只有 STIM/CTRL 两大池，直接当样本用每组只有 1 个 pseudobulk、无 replicate 可拟合）。供者编码在 GEO 的逐细胞元数据 tsne.df 的 `ind` 列，**01 读入时**按 (group, barcode) join 写入 `sample_id`（merge 前的干净条码，无需处理数字后缀；条码在两组间大量重复，必须带组匹配）。有原生 sample_id 列的数据集把 `multi$cell_metadata` 置空，样本表每行即一个生物学样本
- **守卫**：`min_cells`（pseudobulk 最小细胞数）、每组最少样本数、残差自由度（配对设计下 2v2 会参数饱和）；不满足的类型跳过并记入 `summary_degs.csv` 的 skip_reason
- **结果判读**：05 是探索层、06 是确认层——同 padj 阈值下 06 的显著 DEG 应**更少更保守**，顶层 `summary_degs.csv` 与 `summary_05_vs_06.pdf` 直接对照。全基因结果（含独立过滤 padj=NA 行的 stat）保留，`run_gsea: true` 时按 Wald stat 排序喂 GSEA（复用 diff_gsea 的 gmt）

## 高级模块（run_pipeline.sh --with-advanced）

- **07 基因集打分**：AddModuleScore（基因集 vs 随机对照集的平均表达差，推断通路活性）；UCell/AUCell 等基于排名的方法见扩展学习脚本（未随本仓库上传）
- **08 拟时序（monocle3）**：学习分化/变化轨迹图 + 伪时间排序，研究连续变化过程（如 T 细胞分化、刺激响应）；适合具有连续谱系结构的数据。参考 [hossainlab/sc-workflow](https://github.com/hossainlab/sc-workflow) 的 `05_trajectory_analysis` 模块
  - **谱系预设**（`config/trajectory_presets/`）：轨迹分析需要两个先验——哪些细胞构成连续过程、谁是起点。预设文件 = root + 细胞子集 + 文献依据，主 config 用 `trajectory$preset` 引用；留空 = 全部细胞（无先验）。写作清单见该目录 README
  - **轨迹空间**（`trajectory$reduction`）：`raw` = monocle3 自算 PCA/UMAP（保留条件信号，适合批次=生物学条件，如 STIM/CTRL）/ `harmony` = 注入校正嵌入（技术批次）/ `umap` = 复用 Seurat 布局（坐标与 DimPlot 完全一致，校正性取决于 Seurat 端来源）
  - 每次运行写 `trajectory_resolved.yaml` 快照，记录实际生效的预设/根/子集/轨迹空间/Inf 细胞数，随时可追溯
- **09 细胞通讯（CellChat）**：基于配体-受体数据库推断细胞群之间的通讯强度与通路。
  - 单样本：整体通讯网络（数量/强度网络图、信号角色热图、通路气泡图）
  - 多样本：两组通讯差异比较（如 STIM vs CTRL 的 IFN 信号差异，CellChat 官方教程的经典用法）
  - 参考 [hossainlab/sc-workflow](https://github.com/hossainlab/sc-workflow) 的 `06_cell_communication` 模块
  - CellChat 官方教程 vignette：https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html

> ✅ 08 拟时序（monocle3）与 09 细胞通讯（CellChat）均已装包；08 已全量验证（run02 批次：T 细胞谱系预设，14,214 细胞，Inf=0）。拟时序另有独立学习项目（monocle3/velocyto/scvelo 已真实跑通，未随本仓库上传）。

> 多种整合方法的对比**不单独设脚本**——主线 `02_preprocess_cluster.R` 的多样本分支已支持 harmony/RPCA/FastMNN，改 config 的 `integrate$method` 即可切换；各方法原理对比见扩展学习脚本（未随本仓库上传）。不做功能重复的脚本，保持流程精简。

## 数据来源

| 数据 | 来源 | 用途 |
|------|------|------|
| PBMC 3k（hg19 三件套） | 10x 官网 | 单样本流程 |
| STIM/CTRL 表达矩阵 | GSE96583（IFN-β 刺激 PBMC） | 多样本整合 + 差异 + GSEA + pseudobulk 差异（注：该 GEO 矩阵在提交时被过滤掉了线粒体基因，percent.mt 指标退化为 0；真实 cellranger 输出含 MT 基因，不受影响）。供者元数据 `GSE96583_batch2.total.tsne.df.tsv.gz`（`ind` 列 = 供者，另含作者注释与 doublet 标记）在同目录，是 06 pseudobulk 的 sample_id 来源 |
| pbmc_1k_v3 FASTQ + GRCh38-2024-A | 10x 官网 | cellranger 真实数据练习 |
| 基因集 gmt | MSigDB（KEGG/Hallmark/c7） | GSEA 与基因集打分 |

## 输出产物（按批次归档）

- **批次（batch）**：config 的 `batch` 字段定义本次运行的归档目录。同批次名重跑 = 断点续跑；换批次名 = 全新运行，旧批次结果完整保留
- 每步独立目录：`output/<batch>/<步骤名>/`（结果 rds + `.done` 断点标记 + figures/）
- **参数快照**：每个批次自动保存 `output/<batch>/config_used.yaml`——这个批次是用什么参数跑出来的，随时可追溯
- 图片统一 `utils/save_fig.R`（cairo_pdf 支持中文），自动输出到各步的 `figures/`
- `output/<batch>/logs/<步骤名>.log`：每个批次每步的运行日志（批次内归档）
