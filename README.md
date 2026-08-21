# scrna-seq-pipeline

> 从上游定量到下游分析的 scRNA-seq 完整分析流程。
> 在没有真实项目数据的情况下，基于公开教程 + 10x 官方文档自学搭建，Seurat v5 适配。
> 架构参考 GitHub 主流 scRNA-seq 项目（每步一脚本 + 集中配置 + 顺序 runner，如 hossainlab/sc-workflow；工业级工作流引擎 Snakemake/Nextflow 的对比见 docs/01）。

## 项目故事（怎么学出来的）

1. **跟教程学**：基于 Seurat v4 公开教程，从单样本基础分析到多样本整合分析，共 10 个教学脚本（原始学习轨迹保留在本机，未随本仓库上传）
2. **适配新版本**：教程基于 Seurat v4，实际环境是 Seurat v5 + Linux，逐行改写并踩坑记录（v4/v5 差异见 [docs/Seurat_v4_vs_v5.md](docs/Seurat_v4_vs_v5.md)）
3. **补齐上游**：教程从表达矩阵开始，自学补充了 cellranger 定量部分（本机已真实跑通 tiny 测试数据，产出位于 `output/run01/00_cellranger/tiny_demo/`）
4. **扩展学习**：在主线之外自学了多种整合方法（RPCA/FastMNN/scVI/Liger）和基因集打分（AddModuleScore/UCell/AUCell）
5. **工程化整合**：调研 GitHub 主流项目结构后，把流程重构为**分阶段脚本 + config 集中配置 + 断点续跑 runner**（脚本粒度参考 hossainlab/sc-workflow：按分析阶段打包，而非每个操作一个脚本）

## 目录结构

```
scrna-seq-pipeline/
├── README.md                    # 项目总览（本文件）
├── environment.yaml             # mamba 环境定义（一键重建 scrna 环境）
├── run_pipeline.sh              # 总控：顺序执行 + 断点续跑 + 每步日志
├── config/
│   ├── config.yaml              # 单样本流程配置（默认）
│   ├── config.multi.yaml        # 多样本流程配置（--mode multi 自动加载）
│   ├── sample_sheet.csv         # 多样本样本表（sample, group, fastqs, matrix）
│   ├── cellranger.yaml          # 上游定量工具级配置（与下游 R 解耦，跨项目复用）
│   ├── trajectory_presets/      # 07 拟时序谱系预设（root + 细胞子集 + 文献依据）
│   └── annotation_presets/      # 03 注释 marker 预设库（物种文件夹 + 别名 + 出处）
├── bash/                        # 上游 shell 脚本（cellranger 是 shell 工具，语言无关）
│   └── 00_run_cellranger.sh     # cellranger count 批量定量（读 config + 样本表）
├── R/                           # R 流程脚本（分阶段，每阶段可单独重跑/投递节点）
│   ├── 01_load_qc.R             # 读入 + 质控（单/多样本按 config 自动分支）
│   ├── 02_preprocess_cluster.R  # 预处理 + 去批次 + UMAP + 聚类
│   ├── 03_annotate.R            # 细胞注释 + 亚群细分（config 开关）
│   ├── 04_multi_group_plot.R    # 多组可视化
│   ├── 05_diff_gsea.R           # 差异分析 + GSEA
│   ├── 06_gene_set_score.R      # 高级：基因集打分
│   ├── 07_trajectory_analysis.R # 高级：拟时序（monocle3）
│   └── 08_cell_communication.R  # 高级：细胞通讯（CellChat）
├── Python/                      # Python 流程（规划中：scanpy 平行版本，见 Python/README.md）
├── docs/                        # 流程文档 / cellranger 指南 / R 包编译踩坑
├── utils/                       # 共享工具：utils.R（配置/归档/质控/预处理/聚类）
│                                #          + save_fig.R（图片保存，cairo_pdf 支持中文）
├── data/                        # 数据分层管理（软链接/大文件，不入 git）：
│   ├── raw/                     #   原始 FASTQ（只读，cellranger 输入）
│   ├── matrices/                #   公开数据集/已定量矩阵（PBMC3k、GSE96583）
│   └── reference/               #   参考数据：
│       ├── refdata/             #     cellranger 参考基因组（tiny 演示参考已就位）
│       └── gene_sets/           #     MSigDB 基因集 gmt
└── output/<batch>/              # 每个批次完全自包含（不入 git）：
                                 #   各步骤结果 rds + figures/ + logs/（每步日志）
                                 #   + config_used.yaml（本批次参数快照，可追溯）
```

## 分析流程总览

```mermaid
flowchart LR
    A[FASTQ<br/>10x 原始测序数据] -->|00 cellranger count| B[表达矩阵<br/>barcodes/features/matrix]
    B -->|01 读入 + 质控<br/>+ 双细胞检测| C[Seurat 对象<br/>过滤后]
    C --> D[02 标准化 + 高变基因<br/>ScaleData + PCA]
    D --> E[02 UMAP + 聚类<br/>FindNeighbors/FindClusters]
    E --> F[03 细胞注释<br/>marker 鉴定 + RenameIdents]
    F --> G[03 亚群细分<br/>FindSubCluster]
    G --> H[04 多组比较可视化<br/>比例/DotPlot/热图]
    H --> I[05 差异分析 + GSEA<br/>FindMarkers/compareCluster]

    subgraph 多样本
        C2[多样本 merge] -->|02 harmony / RPCA / FastMNN 去批次| E
    end

    subgraph 高级模块（--with-advanced）
        F -.-> J[06 基因集打分<br/>AddModuleScore]
        F -.-> K[07 拟时序<br/>monocle3]
        F -.-> L[08 细胞通讯<br/>CellChat]
    end
```

## 快速开始

```bash
# 1. 激活环境（首次使用先 mamba env create -f environment.yaml 重建）
mamba activate scrna        # R 4.5.3 + Seurat 5.5.1 + presto 加速

# 2. 上游定量（从 FASTQ 开始；需要先下载参考基因组并填写样本表 fastqs 列）
bash bash/00_run_cellranger.sh
#    工具级参数: config/cellranger.yaml（参考基因组/硬件资源）
#    参考基因组下载与 FASTQ 组织: 见 docs/02_cellranger_guide.md（参考 ~11GB）
#    ⚠️ 没有 FASTQ 时跳过此步——示例数据自带矩阵（data/matrices/ 公开数据集），直接从第 3 步开始

# 3. 单样本流程（PBMC3k 示例，默认配置 config/config.yaml）
bash run_pipeline.sh

# 4. 多样本流程（STIM/CTRL，自动加载 config/config.multi.yaml）
bash run_pipeline.sh --mode multi

# 5. 追加高级模块
bash run_pipeline.sh --mode multi --with-advanced

# 6. 灵活执行（同批次断点续跑：已有产物的步骤自动跳过）
bash run_pipeline.sh --only 02            # 只重跑预处理+聚类（改 resolution 后常用）
bash run_pipeline.sh --from 02 --to 03    # 跑 02~03
bash run_pipeline.sh --force              # 当前批次全部强制重跑
Rscript R/02_preprocess_cluster.R   # 单阶段直接运行（HPC 节点上单独投递）

# 7. 批次管理（config 里 batch 字段）：每次运行的结果与日志按批次归档
#    同批次名重跑 = 断点续跑；换批次名（如 batch: run02）= 全新运行，旧批次结果保留
#    每个批次自动保存 config_used.yaml 参数快照，便于对比追溯
#    output/run01/（结果 + logs/ + config_used.yaml）← 一个完整批次
```

## 关键踩坑记录

| 坑 | 原因 | 解决 |
|----|------|------|
| merge 后 `FindAllMarkers` 报错 | Seurat v5 merge 产生 split layers | `JoinLayers()` |
| 加载 v4 格式 rds 后 `JoinLayers` 报错 | v4 对象 Assay 类是 "Assay" | `scobj[["RNA"]] <- as(scobj[["RNA"]], "Assay5")` |
| cellranger 9 报缺 `--create-bam` | Cell Ranger 9 起该参数必填 | 显式指定 `--create-bam=false` |
| scale.data 体积爆炸 | 全基因缩放 = 基因数×细胞数稠密矩阵 | 默认只缩放高变基因；保存前清空 |
| top marker 不含经典 marker | v5 logFC 公式变化，低表达基因 logFC 虚高 | `annotate$marker_pct` 预筛选后再排序（默认 0.5） |
| GSEA 结果不可复现 | GSEA 用随机置换算 p 值 | `set.seed(123)` |
| 单样本/多样本注释映射混用 | 不同数据集分群数和细胞类型不同 | 配置 profile 分离：config.yaml / config.multi.yaml |
| monocle3 `preprocess_cds` 段错误 | irlba fastpath + OpenBLAS 多线程竞态（实测 r-irlba 2.3.7 + OpenBLAS 0.3.33） | run_pipeline.sh 对 07 步骤单独 `OPENBLAS_NUM_THREADS=1`（R 内 `Sys.setenv` 无效，必须 shell 层） |
| `order_cells` 报 "root_pr_nodes or root_cells must be provided" | Rscript 非交互模式无法弹窗选根（RStudio 交互才会弹） | 谱系预设显式指定根（config/trajectory_presets/） |
| 全细胞混跑轨迹 30% 细胞伪时间 Inf | 异构细胞类型间不存在连续轨迹，图不连通 | 谱系预设子集（如仅 T 细胞），run02 验证 Inf=0 |
| FeaturePlot 与 DimPlot 布局不一致 | `umap_naive` 抢占默认降维（DefaultDimReduc 按名匹配），FeaturePlot 不显式指定时画在校正前的图上 | 所有 FeaturePlot 显式 `reduction = "umap"` |

## GitHub 上传说明

- 仓库只提交**代码和文档**（`R/`、`bash/`、`Python/`、`utils/`、`config/`、`docs/`、`run_pipeline.sh`、`environment.yaml`、`LICENSE`、README）
- 数据和结果（`data/`、`output/`、软链接）已写入 `.gitignore`
- 数据来源与获取方式见 `docs/01_pipeline_overview.md` 末尾

## 环境

- 完整依赖清单见 [environment.yaml](environment.yaml)，一条命令重建环境：
  `mamba env create -f environment.yaml`
- 核心版本：R 4.5.3、Seurat 5.5.1、SeuratObject 5.4.0、monocle3 1.4.27、harmony、scDblFinder、yaml
- ⚠️ 4 个包需从 GitHub 安装（conda 渠道无，见 yaml 头部注释）：SeuratWrappers、CellChat（必需）+ presto、scCustomize（可选）
- cellranger 9.0.1（独立软件，非 R 包，见 docs/02_cellranger_guide.md）

### 高级模块依赖（跑 07/08 前检查）

```bash
# 07 拟时序：monocle3 已安装并全量验证（run02 批次，T 细胞谱系预设）
#   安装：mamba install -c conda-forge -c bioconda r-monocle3
#   注意：运行需单线程 BLAS（run_pipeline.sh 已自动处理，见踩坑表）

# 08 细胞通讯：CellChat 2.2.0.9001 已安装（GitHub 源码编译；CRAN 已下架）
#   安装方式（备查）：devtools::install_github("jinworks/CellChat")
#   编译要点见 docs/03_r_package_compile.md（conda 环境 + R 4.5 头文件）
```
