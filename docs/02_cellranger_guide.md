# cellranger 上游定量教程

> 10x Chromium 测序数据的标准定量工具。本教程覆盖：安装、参考基因组、count 命令逐参数讲解、产出解读。
> 本机已用自带 tiny 测试数据真实跑通：`output/run01/00_cellranger/tiny_demo/outs/`

## 1. 安装与验证

```bash
# 下载（本机已装 9.0.1）
# https://www.10xgenomics.com/support/software/cell-ranger/downloads
tar -xzvf cellranger-9.0.1.tar.gz -C ~/Software

# 验证
cellranger --version
cellranger testrun --id=tiny   # 自带测试数据，验证安装完整（跑一次 ~5 分钟）
```

## 2. 参考基因组

cellranger count 需要"转录组参考"，两种获取方式：

```bash
# 方式一（推荐）：10x 预构建，下载解压即用（~11GB）
# 与 cellranger 9 配套版本 GRCh38-2024-A：
# https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz

# 方式二：自建（genome fasta + gene annotation gtf，需数小时）
cellranger mkref \
  --genome=GRCh38_custom \
  --fasta=GRCh38.fa \
  --genes=GRCh38.gtf
```

**注意**：cellranger 版本与参考版本要配套（9.x 配 2024-A，7.x/8.x 配 2020-A 等），不配套会直接报错。

**本项目放置位置**：下载解压到 `data/reference/refdata/`（该目录不入 git），然后在 `config/cellranger.yaml` 的 `transcriptome` 字段填入 `data/reference/refdata/refdata-gex-GRCh38-2024-A`。

## 3. count 命令逐参数讲解

真实数据示例（10x 官方 pbmc_1k_v3，FASTQ ~600MB）：

FASTQ 获取与放置（本项目约定）：
1. 下载 pbmc_1k_v3 FASTQ（~600MB）：https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
2. 解压到 `data/raw/pbmc_1k_v3_fastqs/`（文件名需符合 10x 命名规范：`*_S*_L*_R1_001.fastq.gz`，R1 = 细胞 barcode + UMI，R2 = cDNA 插入片段）
3. 在 `config/sample_sheet.csv` 的 `fastqs` 列填该目录路径，运行 `bash bash/00_run_cellranger.sh` 批量定量

手动执行单样本 count 的命令如下：

```bash
cellranger count \
  --id=pbmc_1k_v3 \
  --transcriptome=refdata-gex-GRCh38-2024-A \
  --fastqs=pbmc_1k_v3_fastqs \
  --expect-cells=1000 \
  --create-bam=false \
  --localcores=16 \
  --localmem=8
```

| 参数 | 含义 | 说明 |
|------|------|------|
| `--id` | 输出目录名 | 一般用样本名，产出 `<id>/outs/` |
| `--transcriptome` | 参考基因组路径 | 上面下载/构建的目录 |
| `--fastqs` | FASTQ 目录 | 10x 双端：**R1 = 细胞 barcode + UMI**（26bp），**R2 = cDNA 插入片段**；自动识别目录下所有样本 |
| `--expect-cells` | 预估细胞数 | 帮助算法从"细胞 barcode vs 背景噪声"中划阈值；**宁多勿少**，估少了真细胞会被丢掉 |
| `--create-bam` | 是否输出 BAM | **Cell Ranger 9 起必填**；做 velocyto RNA 速率等下游分析需 `true`，否则 `false` 省大量空间 |
| `--localcores` | 并行核心数 | 16 核内性能近似线性 |
| `--localmem` | 内存上限 GB | 人参考基因组建议 ≥8GB |

其他常用参数：`--sample`（FASTQ 目录有多个样本时指定）、`--chemistry`（自动检测失败时手动指定 10x 化学版本）、`--nosecondary`（跳过分群步骤，纯定量提速）。

## 4. 产出解读（`<id>/outs/`）

| 文件 | 内容 |
|------|------|
| **web_summary.html** | 质检报告：细胞数、每个细胞 reads 数、测序饱和度、Q30、比对率，附聚类预览（浏览器打开，面试可直接展示） |
| **filtered_feature_bc_matrix/** | 过滤后的表达矩阵（barcodes.tsv + features.tsv + matrix.mtx）——**下游 Seurat 的输入** |
| filtered_feature_bc_matrix.h5 | 同上，h5 单文件格式 |
| raw_feature_bc_matrix/ | 原始矩阵（未做 barcode 过滤） |
| cloupe.cloupe | Loupe Browser 交互可视化文件 |
| metrics_summary.csv | 定量指标汇总表 |
| molecule_info.h5 | 分子水平信息（velocyto 等下游需要） |

## 5. 本项目实际运行结果

用 cellranger 自带 tiny 测试数据（2 个基因的微型基因组 + 2 条 lane 的模拟 FASTQ）在本机真实跑通了完整 count 流程：

```bash
# output/run01/00_cellranger/count.log 记录了完整执行过程
cellranger count --id=tiny_demo \
  --transcriptome=/home/yangcl/Software/cellranger-9.0.1/external/cellranger_tiny_ref \
  --fastqs=/home/yangcl/Software/cellranger-9.0.1/external/cellranger_tiny_fastq \
  --create-bam=false --localcores=16 --localmem=2
```

产出见 `output/run01/00_cellranger/tiny_demo/outs/`（web_summary.html 可直接打开给面试官看）。

> tiny 数据仅用于验证流程跑通；真实数据分析请按第 3 节命令 + pbmc_1k_v3 数据执行，预计 20~40 分钟（依赖机器配置）。

## 6. 常见报错

| 报错 | 原因 | 解决 |
|------|------|------|
| `--create-bam` missing | Cell Ranger 9 起必填 | 显式指定 true/false |
| Reference version mismatch | 参考与 cellranger 版本不配套 | 下载配套版本参考 |
| `No input FASTQ was found` | 文件名不符合 10x 命名规范 | 用 `--sample` 指定，或检查 `*_S*_L*_R1_001.fastq.gz` 命名 |
| 内存不足被杀 | 参考大、细胞多 | 提高 `--localmem`；或减少 `--localcores` |
