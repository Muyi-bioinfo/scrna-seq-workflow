###############################################################################
### 02_preprocess_cluster.R — 模块 3+4+5：预处理 → 去批次（多样本）→ UMAP → 聚类
###
### multi 模式的去批次方法由 integrate$method 控制（harmony / rpca / mnn）
### 输入：output/<batch>/01_load_qc/seurat_qc.rds
### 输出：output/<batch>/02_preprocess_cluster/seurat_clustered.rds + 各阶段图
###############################################################################

source("utils/utils.R")
cfg <- load_config()
step_dir <- setup_step("02_preprocess_cluster")
log_step(paste("模块 3+4+5 | 预处理 + 聚类 | 模式:", cfg$mode,
               if (cfg$mode == "multi") paste("| 去批次:", cfg$integrate$method) else ""))

scobj <- read_prev("01_load_qc", "seurat_qc.rds")

# ---- 第一步：预处理（标准化/高变基因/缩放/PCA，函数在 utils/utils.R）----
scobj <- preprocess_scobj(scobj, cfg)

# ---- 第二步：多样本去批次 ----
if (cfg$mode == "multi") {
  method <- cfg$integrate$method

  if (method == "harmony") {
    # ⚠️ v5 关键：merge 后产生 split layers，走 harmony 前必须 JoinLayers
    scobj <- JoinLayers(scobj)

    # 校正前的 UMAP 对照图（维度与整合后一致，前后对比才严谨）
    # ⚠️ umap_naive 会抢占默认降维：后续不显式指定 reduction 的 FeaturePlot
    #   会画在校正前的图上（DefaultDimReduc 按名匹配），见 README 踩坑表
    scobj <- RunUMAP(scobj, reduction = "pca",
                     dims = dims_from(cfg$integrate), reduction.name = "umap_naive")
    save_fig(DimPlot(scobj, reduction = "umap_naive", group.by = "group"),
             "umap_before_integration", type = "dimplot")

    library(harmony)
    scobj <- RunHarmony(scobj, reduction = "pca",
                        group.by.vars = "group", reduction.save = "harmony")
    reduction_use <- "harmony"
  } else {
    # Seurat v5 推荐：IntegrateLayers 统一接口，保留 split layers 直接整合
    library(SeuratWrappers)
    method_fun <- switch(method,
                         rpca = RPCAIntegration,   # 内存低速度快，适合相似数据集
                         mnn  = FastMNNIntegration, # 批次效应复杂时更稳健
                         stop("未知的整合方法: ", method, "（支持 harmony / rpca / mnn）"))
    scobj <- IntegrateLayers(scobj, method = method_fun,
                             orig.reduction = "pca",
                             new.reduction = paste0("integrated.", method),
                             verbose = FALSE)
    scobj[["RNA"]] <- JoinLayers(scobj[["RNA"]])   # 整合后合并 layers（官方写法）
    reduction_use <- paste0("integrated.", method)
  }
} else {
  reduction_use <- cfg$cluster$reduction   # 单样本直接用 PCA
}

# ---- 第三步：UMAP + 聚类（函数在 utils/utils.R，维度/分辨率来自 config）----
# 多样本用 integrate 段的维度与分辨率，单样本用 cluster 段
cluster_cfg <- if (cfg$mode == "multi") cfg$integrate else cfg$cluster
scobj <- cluster_cells(scobj,
                       reduction = reduction_use,
                       dims = dims_from(cluster_cfg),
                       resolution = cluster_cfg$resolution,
                       fig_name = "dimplot_clusters")

# 多样本额外画按组着色的 UMAP，检验批次效应是否消除
if (cfg$mode == "multi") {
  save_fig(DimPlot(scobj, reduction = "umap", group.by = "group", label = TRUE),
           paste0("dimplot_by_group_", cfg$integrate$method), type = "dimplot")
}

slim_save(scobj, "seurat_clustered.rds", step_dir)
log_step("模块 3+4+5 完成 → 查看 output/02_preprocess_cluster/figures/dimplot_clusters.pdf 确认分群后注释")
