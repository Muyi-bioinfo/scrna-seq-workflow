###############################################################################
### 07_trajectory_analysis.R — 高级模块：拟时序分析（monocle3）
###
### 输入：output/<batch>/03_annotate/seurat_annotated.rds（需 celltype 列与 umap 降维）
### 输出：output/<batch>/07_trajectory_analysis/（轨迹图 + 伪时间图 + cds 对象 rds）
###
### ⚠️ 坑：irlba + OpenBLAS 多线程必段错误（实测 r-irlba 2.3.7 + OpenBLAS 0.3.33），
###        需 OPENBLAS_NUM_THREADS=1（run_pipeline.sh 已自动处理；RStudio/手动 Rscript
###        需自己设，R 内 Sys.setenv 无效）
### 谱系预设（root + 细胞子集）与轨迹空间三模式说明见 docs/01 高级模块一节
###############################################################################

source("utils/utils.R")
cfg <- load_config()
step_dir <- setup_step("07_trajectory_analysis")
log_step("高级模块 | 拟时序分析（monocle3）")
library(monocle3)
library(SeuratWrappers)

scobj <- read_prev("03_annotate", "seurat_annotated.rds")
Idents(scobj) <- "celltype"

# ---- 1. 加载轨迹预设（可选）：谱系先验 = root + 细胞子集 ----
# 优先级：主 config 显式参数 > 预设 > 默认值（见 utils.R 的 load_trajectory_preset）
traj <- load_trajectory_preset(cfg$advanced$trajectory)

# ---- 2. 谱系子集（可选）----
# ⚠️ 异构细胞全混跑会因轨迹图不连通产生大量 Inf 伪时间（run02 实测 30%），
#    论文惯例是只对构成连续过程的细胞跑（如仅 T 细胞）
scobj <- subset_by_celltypes(scobj, traj$subset_celltypes)

# ---- 3. Seurat 对象 → monocle3 CDS 对象 ----
# as.cell_data_set 使用 RNA assay 的 data 层（已归一化表达）
cds <- as.cell_data_set(scobj, assay = "RNA")

# ---- 4. 轨迹空间选择（reduction：在哪个空间建图）----
# raw=自算 PCA/UMAP 保全部信号（批次=生物学条件时用）/ harmony=注入校正嵌入（技术批次）/
# umap=复用 Seurat UMAP（坐标与 DimPlot 一致）；完整决策说明见 docs/01
num_dim <- if (!is.null(traj$num_dim)) traj$num_dim else 50
reduction <- if (!is.null(cfg$advanced$trajectory$reduction)) cfg$advanced$trajectory$reduction else "raw"
if (reduction == "raw") {
  cds <- preprocess_cds(cds, num_dim = num_dim)
  cds <- reduce_dimension(cds)
} else if (reduction == "harmony") {
  if (!"harmony" %in% Reductions(scobj)) {
    stop("reduction=harmony 需要 Seurat 对象有 harmony 降维（多样本流程需先跑 02 整合）")
  }
  reducedDims(cds)[["PCA"]] <- Embeddings(scobj, "harmony")   # 校正后的嵌入当 PCA 用
  cds <- reduce_dimension(cds)
} else if (reduction == "umap") {
  if (!"umap" %in% Reductions(scobj)) {
    stop("reduction=umap 需要 Seurat 对象有 umap 降维")
  }
  reducedDims(cds)[["UMAP"]] <- Embeddings(scobj, "umap")     # 直接复用分群图坐标
} else {
  stop("未知 reduction: ", reduction, "（支持 raw / harmony / umap）")
}
message("-- 轨迹空间: ", reduction)

# ---- 5. 聚类 + 学习轨迹图 ----
# ⚠️ utils/utils.R 定义了 Seurat 版 cluster_cells()（需要 reduction 参数），
#   会遮蔽 monocle3 的同名函数，必须用 monocle3:: 显式调用
cds <- monocle3::cluster_cells(cds)
cds <- learn_graph(cds)

# ---- 6. 伪时间分配：根决定轨迹方向（哪个细胞类型是"起点"）----
# ⚠️ Rscript 非交互模式下 monocle3 无法弹出选根窗口，必须显式指定根，
#   否则报错 "either root_pr_nodes or root_cells must be provided"
root_ct <- traj$root_celltype
if (is.null(root_ct) || !nzchar(root_ct)) {
  # 未配置根细胞类型时退化为最大细胞群（方向无生物学意义，仅保证流程可跑通）
  root_ct <- names(sort(table(cds$celltype), decreasing = TRUE))[1]
  message("-- 未配置 root_celltype，自动使用最大细胞群作为根（方向无生物学意义）: ", root_ct)
}
root_cells <- colnames(cds)[cds$celltype == root_ct]
if (length(root_cells) == 0) {
  stop("根细胞类型不存在: ", root_ct,
       "，请检查 config 的 advanced$trajectory$root_celltype 或预设文件")
}
message("-- 根细胞类型: ", root_ct, "（", length(root_cells), " 个细胞）")
cds <- order_cells(cds, root_cells = root_cells)

# ---- 7. 参数快照（可追溯：这个批次用了什么先验、什么轨迹空间）----
resolved <- list(
  preset            = traj$preset_name,
  description       = traj$description,
  reduction         = reduction,                                # 轨迹空间：raw / harmony / umap
  root_celltype     = root_ct,
  subset_celltypes  = as.list(sort(unique(scobj$celltype))),   # 实际参与轨迹的细胞类型
  num_dim           = if (reduction == "raw") num_dim else NA, # 仅 raw 模式使用
  n_cells           = ncol(cds),
  n_inf_pseudotime  = sum(is.infinite(pseudotime(cds)))
)
yaml::write_yaml(resolved, file.path(step_dir, "trajectory_resolved.yaml"))
message("-- 参数快照: ", file.path(step_dir, "trajectory_resolved.yaml"))

# ---- 8. 可视化 ----
# 伪时间着色：颜色从深到浅代表分化/变化方向
save_fig(plot_cells(cds, color_cells_by = "pseudotime",
                    label_groups_by_cluster = FALSE,
                    label_leaves = FALSE, label_branch_points = FALSE),
         "trajectory_pseudotime", type = "dimplot")

# 细胞类型着色 + 轨迹图：验证轨迹方向是否与已知生物学一致
save_fig(plot_cells(cds, color_cells_by = "celltype",
                    label_cell_groups = FALSE,
                    label_leaves = FALSE, label_branch_points = FALSE),
         "trajectory_celltype", type = "dimplot")

# 多样本数据按分组着色，观察两组沿轨迹的分布差异
if (cfg$mode == "multi") {
  group_col <- cfg$multi_group$group_col
  save_fig(plot_cells(cds, color_cells_by = group_col,
                      label_cell_groups = FALSE,
                      label_leaves = FALSE, label_branch_points = FALSE),
           "trajectory_by_group", type = "dimplot")
}

# ---- 9. 保存 CDS 对象（供下游如拟时序差异基因分析）----
saveRDS(cds, file = file.path(step_dir, "trajectory_cds.rds"))

log_step("高级模块 | 拟时序分析完成 → 结果见 output/07_trajectory_analysis/")
