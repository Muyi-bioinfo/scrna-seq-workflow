# Seurat v4 → v5 主要区别

> 自用参考文档。侧重概念理解（为什么这么改），辅以代码对照（怎么改写）。
> 日常写脚本时的速查表见 [README 踩坑表](../README.md#关键踩坑记录)。

## 目录

- [1. Layers 架构（Assay5）](#1-layers-架构assay5)
- [2. 整合流程重构（IntegrateLayers）](#2-整合流程重构integratelayers)
- [3. 差异分析变化](#3-差异分析变化)
- [4. 对象迁移实操](#4-对象迁移实操)
- [5. SCTransform v2](#5-sctransform-v2)
- [6. 可视化变化](#6-可视化变化)

---

## 1. Layers 架构（Assay5）

### 一句话总结

v5 用 **Layers** 替代了 v4 的 **Slots**，一个 Assay 内可以存储多个同类型的矩阵（如 `counts.sample1`、`counts.sample2`），并且原生支持磁盘矩阵（BPCells），突破了内存上限。

### 为什么这么改

v4 的 Assay 结构是"三个固定格子"：

```
v4 Assay
├── @counts      # 只能存一份原始计数矩阵
├── @data        # 只能存一份归一化矩阵
└── @scale.data  # 只能存一份缩放矩阵
```

当处理多样本数据时，你只能选择：要么把所有样本 merge 成一个大矩阵塞进唯一的 `@counts`，要么把样本拆成多个独立的 Seurat 对象分别管理。前者丢失了样本来源信息，后者让跨样本分析变得繁琐。

v5 的 Assay5 把固定 slot 改成了可命名的 **layer**：

```
v5 Assay5
├── counts.ctrl    # CTRL 样本的原始计数
├── counts.stim    # STIM 样本的原始计数
├── data           # 归一化数据（JoinLayers 后合并）
└── scale.data     # 缩放数据
```

核心收益：
- **数据不丢失来源**：每个样本/条件的 counts 可以独立存储，需要时再合并
- **按需加载**：配合 BPCells 包，矩阵可以存在磁盘上，分析时按需读取子集，处理百万级细胞不再需要几百 GB 内存
- **灵活扩展**：可以自己创建任意命名的 layer（如 `counts.batch1`、`data.sct`）

### 代码对照

```r
# 读取数据（layer 参数替代 slot）
LayerData(obj, assay = "RNA", layer = "counts")
# v4: GetAssayData(obj, assay = "RNA", slot = "counts")

# 访问 layer（$ 访问器）
obj[["RNA"]]$counts      # counts layer
obj[["RNA"]]$data        # data layer
obj[["RNA"]]$scale.data  # scale.data layer
# v4: obj@assays$RNA@counts / obj@assays$RNA@scale.data（V4 习惯用 @ 直接掏 slot）

# 查看所有 layer
Layers(obj[["RNA"]])
# 典型输出: "counts.ctrl" "counts.stim" "data" "scale.data"

# 查看特定 layer 的维度
dim(LayerData(obj, assay = "RNA", layer = "counts.ctrl"))
```

### JoinLayers() — 什么时候调用

当你 `merge()` 多个样本或者 `split` 了 layers 之后，counts 是分散在多个 layer 里的。大多数下游函数（`FindVariableFeatures`、`ScaleData`、`FindMarkers`）期望数据在统一的 `counts` / `data` layer 中。

```r
# merge 多样本后必须做
scobj <- merge(stim_obj, ctrl_obj)
scobj <- JoinLayers(scobj)
# 作用：将 counts.stim + counts.ctrl 合并到统一的 counts layer
```

⚠️ **不要在 split layers 后立即 JoinLayers**：如果你打算用 v5 的 `IntegrateLayers()` 做整合，它需要 layers 按 batch 分开。先整合，再 JoinLayers。

### 踩坑提示

- `GetAssayData(obj, slot = "counts")` 在 v5 中仍然能跑，但会报 deprecation warning。建议一律改用 `layer`。
- `obj@assays$RNA@scale.data` 在 v5 中**不一定能正确访问**，因为 Assay5 的内部 slot 名变了。始终用 `obj[["RNA"]]$scale.data`。
- 如果你用 `as(obj[["RNA"]], "Assay5")` 转换对象后 layers 看起来不对，检查是否需要手动 `JoinLayers()`。

---

## 2. 整合流程重构（IntegrateLayers）

### 一句话总结

v5 用 **一个函数 `IntegrateLayers()`** 替代了 v4 的三步法，整合结果存储为 **dimensional reduction** 而非独立 assay，且 PCA 在整合**之前**做。

> ⚠️ **本项目脚本使用的整合方式**：`RunHarmony()`（见下方说明）。这是 v4 兼容写法，在 v5 中仍然可用，但 Seurat 官方现在推荐 `IntegrateLayers(method = HarmonyIntegration)`。两种方式本章都会介绍。

### 为什么这么改

v4 的整合流程：

```r
# v4 三步法
obj_list <- SplitObject(obj, split.by = "batch")
anchors  <- FindIntegrationAnchors(obj_list)
obj      <- IntegrateData(anchors)
# 产物：一个新的 "integrated" assay（包含 corrected counts）
```

问题：
1. **流程冗长**：三个函数串行，中间产物（anchors）占用大量内存
2. **输出膨胀**：整合后多了一个完整的 "integrated" assay，基因数 × 细胞数的矩阵又来一份
3. **方法不统一**：Harmony、CCA、RPCA 各有自己的函数名和参数，切换成本高
4. **PCA 在整合后**：对 corrected counts 跑 PCA，结果难以和原始表达量建立联系

v5 的设计：

```r
# v5：统一接口（推荐 ⭐）
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)  # 按批次拆分 layers
obj <- IntegrateLayers(
  obj,
  method = CCAIntegration,       # 或 HarmonyIntegration, RPCAIntegration, FastMNNIntegration
  orig.reduction = "pca",
  new.reduction = "integrated.cca"
)
# 产物：一个 reduction（只有 embedding 坐标，没有完整的 corrected counts 矩阵）
```

核心改进：
- **统一入口**：换算法只改 `method` 参数
- **轻量输出**：只存 reduction embedding，节省大量内存
- **PCA 前置**：先用原始数据跑 PCA，整合是对 PC 空间做对齐——分析链路更清晰
- **按 batch 独立预处理**：`split()` 后高变基因筛选在每个 batch 内独立进行，避免"只在某个 batch 中变异的噪音基因"被选入

### IntegrateLayers（推荐） vs RunHarmony（本项目使用）

两种方式在 v5 中都可用，但预处理路径不同，**结果可能有差异**：

| | IntegrateLayers(HarmonyIntegration) ⭐ | RunHarmony() |
|---|---|---|
| **地位** | v5 官方推荐，未来主流 | v4 兼容方式，稳定可靠 |
| **数据准备** | `split(obj[["RNA"]], f = obj$group)` | `merge()` + `JoinLayers()`（无需 split） |
| **高变基因筛选** | **每个 batch 独立筛选**，取并集（避免技术噪音基因） | 全局筛选（merge 后一起做） |
| **输入** | `orig.reduction = "pca"` | `reduction = "pca"` |
| **输出** | `new.reduction = "integrated.harmony"` | `reduction.save = "harmony"` |
| **多批次变量** | 单变量拆分 | `group.by.vars` 可直接指定多列 |
| **本项目** | — | ✅ 使用 |

**为什么本项目脚本用了 RunHarmony**：脚本从 v4 教程迁移而来，保留 v4 兼容的写法以确保教程结果一致。对新手来说，`RunHarmony()` 流程更直观（无需理解 `split()` layers），且 Harmony 官方至今仍然维护和推荐这种方式。

### 代码对照

```r
# ── v4 ──────────────────────────────────────────
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj) %>% ScaleData() %>% RunPCA()

obj_list <- SplitObject(obj, split.by = "batch")
anchors  <- FindIntegrationAnchors(obj_list, dims = 1:30)
obj      <- IntegrateData(anchors, dims = 1:30)

DefaultAssay(obj) <- "integrated"
obj <- ScaleData(obj) %>% RunPCA() %>% RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters()

# ── v5 IntegrateLayers（⭐ 推荐，新项目用这个）───
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj) %>% ScaleData() %>% RunPCA()

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)  # 关键：按 batch 拆分 layers
obj <- IntegrateLayers(
  obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony",
  dims = 1:30
)

obj <- RunUMAP(obj, reduction = "integrated.harmony", dims = 1:30)
obj <- FindNeighbors(obj, reduction = "integrated.harmony", dims = 1:30)
obj <- FindClusters(obj)

# ── v5 RunHarmony（✅ 本项目使用，v4 兼容）───────
obj <- merge(data1, data2)
obj <- JoinLayers(obj)     # merge 后必须
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj) %>% ScaleData(features = rownames(obj)) %>% RunPCA()

library(harmony)
obj <- RunHarmony(obj, reduction = "pca",
                  group.by.vars = "group",
                  reduction.save = "harmony")

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj)
# 项目实际脚本：02_multiple_samples/scripts/02_data_integeration.R
```

### PCA 时机对比

| | v4 | v5（两种方式通用） |
|---|---|---|
| PCA 执行时机 | 整合**后**，对 integrated assay 做 | 整合**前**，对原始 RNA data 做 |
| UMAP 输入 | integrated assay 的 PC | integration reduction 的 embedding |
| 好处/代价 | corrected 矩阵是"干净的"，但丢失与原始表达量的关联 | 链路清晰，统一用 reduction；原始数据质量差时会放大噪声 |

### 踩坑提示

- **FindMarkers 默认 assay 陷阱**：整合后如果你的 active assay 是 "RNA"，`FindMarkers` 会直接用 RNA 的 data layer——这是对的。但如果你在 v4 习惯先 `DefaultAssay(obj) <- "integrated"` 再找 marker，v5 中没有 "integrated" assay 了，必须用 RNA assay + 正确的 group 信息。
- **IntegrateLayers 必须先 split layers**：`split(obj[["RNA"]], f = obj$batch)` 是调用 IntegrateLayers 的前置步骤，忘记 split 会直接报错。
- **两种方式结果可能不同**：IntegrateLayers 按 batch 独立选高变基因（排除单 batch 特异的噪音基因），RunHarmony 用全局高变基因。大多数情况下结果相似，但大规模异质数据中差异明显。本项目用 RunHarmony 产生的 harmony reduction 在全部分析中直接使用即可。
- **RunHarmony 不需要提前 split layers**：项目脚本中 merge → JoinLayers → RunHarmony 是正确的流程。

---

## 3. 差异分析变化

### 一句话总结

v5 的 `FindMarkers` / `FindAllMarkers` 在 logFC 计算方式、默认参数和底层引擎上均有变化，**相同数据跑出来的 marker 列表可能显著不同**。

### logFC 计算：细胞级 → 组级伪计数

这是最隐蔽但影响最大的变化。

```
v4: logFC = log2( mean(cells_in_group1 + pseudocount) / mean(cells_in_group2 + pseudocount) )
     ↑ 对每个细胞加伪计数后再求均值

v5: logFC = log2( (mean(cells_in_group1) + pseudocount) / (mean(cells_in_group2) + pseudocount) )
     ↑ 先求每组的平均表达量，再加伪计数
```

后果：当一个基因只在 group1 表达、在 group2 完全不表达时，v5 的 logFC 会**显著偏高**。v5 开发者自己也承认"这些估计值对低表达基因不稳定"。

**建议**：v5 中筛选 marker 时用更严格的 `logfc.threshold`，不要直接用默认值，例如 `FindAllMarkers(obj, logfc.threshold = 0.25)`（手动改回 v4 的默认值）。

### 默认参数变化

| 参数 | v4 默认 | v5 默认 | 影响 |
|------|---------|---------|------|
| `min.pct` | 0.1 | 0.01 | v5 更宽松，更多基因进入检验 |
| `logfc.threshold` | 0.25 | 0.1 | v5 更宽松，更多基因被判为显著 |
| `pseudocount.use` | 1（细胞级） | 1（组级，行为不同） | 同上，logFC 普遍偏高 |

两个默认参数都更宽松了 + logFC 算法天然偏高 = **v5 返回的 marker 数量可能比 v4 多数倍，且可能包含核糖体蛋白等生物学意义不大的基因**。

### presto 加速

v5 检测到 `presto` 包已安装时，自动使用 presto 的 Wilcoxon 实现，速度大幅提升。你的 `scrna` 环境中已安装 presto，所以 `FindAllMarkers` 会直接用最快的路径。对结果无影响，纯粹是性能优化。

### 伪批量模式（pseudobulk）

v5 推荐在多条件/多样本场景下使用 `AggregateExpression()` + 伪批量 DE，而不是对单细胞做 Wilcoxon：

```r
# 伪批量：先按样本+细胞类型聚合，再做 DE
pseudo <- AggregateExpression(obj, assays = "RNA", group.by = c("celltype", "sample"))
# 然后用 DESeq2 / edgeR / limma 对聚合后的矩阵做差异分析
```

这是 v5 的推荐实践，但属于进阶用法，基础教程中未涉及。

### 代码对照

```r
# 找所有 cluster 的 marker（v5 写法）
DefaultAssay(obj) <- "RNA"  # 确保用的是 RNA assay
all_markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  logfc.threshold = 0.25,   # 建议手动设回 v4 默认
  min.pct = 0.1              # 建议手动设回 v4 默认
)
# v4: 参数名相同，但默认值不同，且当时没有 presto 加速

# 指定 layer（推荐明确写出）
FindMarkers(obj, ident.1 = "B", ident.2 = "T", layer = "data")
# v4: FindMarkers(obj, ident.1 = "B", ident.2 = "T", slot = "data")
```

### 踩坑提示

- 如果你发现 v5 的 marker 列表和教程/预期差很远，**先把 `logfc.threshold` 和 `min.pct` 改回 v4 默认值再试**。
- 整合后找 marker 前，一定确认 `DefaultAssay(obj)` 是 "RNA" 而不是什么奇怪的 assay。
- `FindConservedMarkers` 在 v5 中有兼容性修复，去掉了 v4 输出中的 NA 行。

---

## 4. 对象迁移实操

### 一句话总结

`UpdateSeuratObject()` **只改版本号**，不转换底层数据结构。要把 v4 对象完全变成 v5 对象，需要额外一步 `as(, "Assay5")`。

### UpdateSeuratObject() 做了什么

```r
obj_v5 <- UpdateSeuratObject(obj_v4)
```

做的事：
- ✅ 将对象的 version 标记从 v4.x 改为 v5.x
- ✅ 更新 DimReduc keys、graph 引用、command log 等元数据
- ❌ **不会**把旧的 `Assay` 转成新的 `Assay5`
- ❌ **不会**创建 layers 结构

检查是否真正转成了 Assay5：

```r
class(obj_v5[["RNA"]])
# 如果输出 "Assay"   → 还是 v4 结构，只是版本号改了
# 如果输出 "Assay5"  → 真正升级完成
```

### 正确的完整转换

```r
# 第一步：更新元数据
obj_v5 <- UpdateSeuratObject(obj_v4)

# 第二步：创建新的 Assay5，手动转换
obj_v5[["RNA5"]] <- as(object = obj_v5[["RNA"]], Class = "Assay5")

# 第三步：检查 layers 是否正常
Layers(obj_v5[["RNA5"]])
# 应该看到: "counts" "data" "scale.data"

# 第四步（可选）：删除旧 assay，重命名
obj_v5[["RNA"]] <- NULL
obj_v5[["RNA"]] <- obj_v5[["RNA5"]]
obj_v5[["RNA5"]] <- NULL
```

⚠️ 注意：转换后对象体积会翻倍（新老 assay 同时存在），第四步清理非常必要。

### 版本共存管理

你的项目通过 conda 环境隔离来管理版本：

```bash
# 当前 scrna 环境是 v5
mamba activate scrna

# 如果需要跑 v4 原版教程（保留不动），可以用独立环境
# 但原版教程的 01_ / 02_ 目录保持原样不动，只在 v5_scripts/ 中改写
```

这种"原版不动 + 新版单独改写"的策略比来回切换 R 版本靠谱得多。

### 踩坑提示

- 如果别人的 `.rds` 文件用 v4 保存，你用 v5 的 `readRDS()` 读进来，结构还是 v4 的。必须先 `UpdateSeuratObject()` + `as(,"Assay5")` 转换。
- `Matrix` 包版本和 `SeuratObject` 版本有耦合，不要混用不同 conda 环境的 `.rds` 文件。
- 转换前后跑 `dim(obj)` 确认细胞数/基因数没变。

---

## 5. SCTransform v2

### 一句话总结

v5 中 `SCTransform()` 默认使用 **v2 算法**（正则化负二项回归），速度更快、内存更省。如果需要复现 v4 结果，设置 `vst.flavor = "v1"`。

### v1 vs v2 区别

| | v1 (vst.flavor = "v1") | v2 (vst.flavor = "v2", 默认) |
|---|---|---|
| 回归模型 | 泊松 + 负二项（两步） | 正则化负二项（一步） |
| 速度 | 较慢，大样本时明显 | 更快，差异在大数据集上显著 |
| 输出 slot | `SCT` assay 的 `scale.data` | 同左，但内部计算方式不同 |
| 重现性 | v1 结果稳定 | v2 结果可能与 v1 有细微差异 |

### 代码对照

```r
# 使用 v2（默认，推荐）
obj <- SCTransform(obj, vst.flavor = "v2")

# 使用 v1（兼容 v4 结果）
obj <- SCTransform(obj, vst.flavor = "v1")

# v4: SCTransform(obj)  ← 当时只有 v1 算法存在
```

### 使用建议

- 新项目一律用 `vst.flavor = "v2"`
- 如果教程/文献基于 v1 结果，需要复现时手动切换到 v1
- 原版教程中如果涉及 SCTransform，改写时默认用 v2 并加注释说明

### 踩坑提示

- SCTransform 产出的 assay 名为 "SCT" 而非 "RNA"，后续所有操作（找高变基因、PCA、UMAP）要用 `DefaultAssay(obj) <- "SCT"`。这个 v4 和 v5 是一样的。
- v2 可能会输出去除的基因略有不同，不影响下游结论。

---

## 6. 可视化变化

### 一句话总结

可视化方面的变化是六章中最小的，主要是**参数名统一**（`slot` → `layer`）和一些**新参数**的加入。v4 的绘图代码绝大多数可以直接在 v5 里跑。

### 新参数

| 函数 | 新增参数 | 作用 | 加入版本 |
|------|---------|------|---------|
| `DimPlot` | `stroke.size` | 控制点的描边粗细 | v5.2 |
| `DimPlot` | `label.size.cutoff` | 只标注大于该尺寸的 cluster | v5.3 |
| `FeaturePlot` | `stroke.size` | 同 DimPlot | v5.3 |
| `VlnPlot` | `raster.dpi` | 栅格化渲染的 DPI | v5.3 |

### slot → layer 迁移

大部分可视化函数接受 `slot` / `layer` 参数来指定用哪个数据。v5 中全部统一为 `layer`：

```r
# FeaturePlot / VlnPlot / DotPlot / DoHeatmap
FeaturePlot(obj, features = "CD3D", layer = "data")
# v4: FeaturePlot(obj, features = "CD3D", slot = "data")

VlnPlot(obj, features = "CD3D", layer = "data")
# v4: VlnPlot(obj, features = "CD3D", slot = "data")
```

### FeaturePlot 新增的 keep.scale 参数

```r
# "feature"：每个基因独立配色（默认）
# "all"：所有基因用统一配色范围
# NULL：不做额外缩放
FeaturePlot(obj, features = c("CD3D", "CD14"), keep.scale = "all")
```

这是一个纯增强，v4 中没有对应参数。

### 踩坑提示

- 如果 `DimPlot` 突然不显示 cluster 标签，检查 metadata 列名是否与降维坐标列名冲突（如不要把细胞命名为 "UMAP_1"）。
- 项目中统一使用 `save_fig()` 工具函数（基于 `cairo_pdf`），和 Seurat 版本无关，不用调整。

---

## 附录：速查清单

改写 v4 脚本为 v5 时，按这个清单逐项检查：

- [ ] `obj@assays$RNA@*` → `obj[["RNA"]]$*`
- [ ] `GetAssayData(obj, slot = "xxx")` → `LayerData(obj, layer = "xxx")`
- [ ] `merge()` 后 → `JoinLayers()`
- [ ] 整合：新项目用 `IntegrateLayers(method = ...)`（需先 split layers）；兼容方式用 `RunHarmony()`（本项目使用）
- [ ] 整合后 UMAP/FindNeighbors/FindClusters 基于 reduction 而非 assay
- [ ] `FindAllMarkers` 的 `logfc.threshold` 和 `min.pct` 手动设回合适的值
- [ ] 旧 `.rds` 文件 → `UpdateSeuratObject()` + `as(, "Assay5")`
- [ ] `SCTransform` 需要 v1 结果时显式指定 `vst.flavor = "v1"`
- [ ] 可视化函数 `slot` 参数 → `layer`

---

> 本文基于 Seurat 5.5.x / SeuratObject 5.4.x，对应项目 `scrna` 环境。
> 搜索参考: [Seurat Official Announcements](https://satijalab.org/seurat/articles/announcements), [SeuratObject NEWS](https://github.com/satijalab/seurat-object/blob/main/NEWS.md), [Seurat GitHub Issues](https://github.com/satijalab/seurat/issues)
