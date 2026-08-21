###############################################################################
### 03_annotate.R — 模块 6+7：细胞注释 + 亚群细分（占 80% 分析时间）
###
### 注释映射在 config 的 annotate$cluster_map，改 config 无需改代码
### 输入：output/<batch>/02_preprocess_cluster/seurat_clustered.rds
### 输出：output/<batch>/03_annotate/seurat_annotated.rds（+ 开 subcluster 时 subclustered.rds）
###############################################################################

source("utils/utils.R")
cfg <- load_config()
step_dir <- setup_step("03_annotate")
log_step(paste("模块 6+7 | 细胞注释",
               if (isTRUE(cfg$subcluster$enable)) "+ 亚群细分" else ""))

scobj <- read_prev("02_preprocess_cluster", "seurat_clustered.rds")
library(ggplot2)   # theme() 等绘图函数（utils.R 只加载了 Seurat/dplyr，ggplot2 需显式加载）

# 注释预设：细胞类型 → marker 基因（物种文件夹 + 存在性校验，见 utils.R）
preset <- load_annotation_preset(cfg$annotate, scobj)

# ---- 第一步：FindAllMarkers（每个群 vs 其余细胞差异检验，presto 自动加速）----
all_markers <- FindAllMarkers(object = scobj, only.pos = TRUE)
saveRDS(all_markers, file = file.path(step_dir, "all_markers_table.rds"))

# ⚠️ pct.1 预筛选（阈值在 config 的 annotate$marker_pct）：
#   v5 logFC 公式变化导致低表达基因 logFC 虚高，不筛会混入假 marker
top10 <- all_markers %>%
  filter(pct.1 > cfg$annotate$marker_pct) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1:10) %>%
  ungroup()
write.csv(top10, file = file.path(step_dir, "top10_markers.csv"), row.names = FALSE)

# ---- 第二步：marker 可视化（双视图，互为补充）----
# 2a. DotPlot 总览：tier1 marker × cluster 矩阵（行按预设的细胞类型分组排列），
#     判断「每个 cluster 是谁」；点大小=表达比例，颜色=平均表达
#     基因先与数据求交集（数据缺失的 marker 不画，如 GEO 矩阵被过滤掉的基因）
p_dot <- DotPlot(scobj, features = intersect(preset$flat_tier1, rownames(scobj)),
                 group.by = "seurat_clusters") +
  RotatedAxis() +                        # 横轴 cluster 标签 45° 斜排，避免重叠
  theme(axis.text.y = element_text(size = 7))
save_fig(p_dot, "dotplot_markers_by_cluster", width = 14, height = 12)

# 2b. 每细胞类型一张 FeaturePlot（该类型全部 tier1）：验证 marker 的
#     空间一致性（这些基因是否集中在 UMAP 同一片区域）；
#     同样与数据求交集，全部缺失的类型跳过（如 GEO 矩阵被过滤掉的基因）
#     单独子目录存放；ncol=2 面板更大；图高按基因行数自适应
if (!is.null(preset$celltypes)) {
  ct_figdir <- file.path(step_dir, "figures", "featureplot_celltypes")
  dir.create(ct_figdir, recursive = TRUE, showWarnings = FALSE)
  old_figdir <- getOption("fig_outdir")
  options(fig_outdir = ct_figdir)
  for (ct in names(preset$celltypes)) {
    ct_genes <- intersect(preset$celltypes[[ct]]$markers$tier1, rownames(scobj))
    if (length(ct_genes) == 0) next
    fname <- gsub("[^A-Za-z0-9]+", "_", ct)   # 类型名含空格/符号时仍可安全作文件名
    # 尺寸：人眼质检工作图（不进论文），以「fit-width 一屏放下」为准——
    # 动态 ncol：≤4 基因 2 列（面板 6×4.8）/ ≥5 基因 3 列（面板 4.7×4.8），
    # 行数封顶 2 行，最坏 14×10.4in 不滚动
    # ⚠️ 必须显式 reduction="umap"：umap_naive 会抢占默认降维（见 README 踩坑表），
    #    不指定会画在校正前的图上，与 dimplot_annotated 布局不一致
    ncol_ct <- if (length(ct_genes) <= 4) 2 else 3
    n_rows <- ceiling(length(ct_genes) / ncol_ct)
    save_fig(FeaturePlot(scobj, features = ct_genes, order = TRUE, ncol = ncol_ct,
                         reduction = "umap"),
             paste0("featureplot_", fname),
             width = if (ncol_ct == 2) 12 else 14,
             height = 4.8 * n_rows + 0.8)
  }
  options(fig_outdir = old_figdir)   # 恢复：后续图（注释 UMAP 等）仍进主 figures/
}

# ---- 第三步：注释 ----
# 确认群的个数：看 02 的 dimplot_clusters + 上面的 marker 图，再改 config 的 cluster_map
Idents(scobj) <- "seurat_clusters"
cluster_map <- unlist(cfg$annotate$cluster_map)   # yaml 的 map 结构 → 命名向量
# 类型名 vs 预设一致性（别名命中只提示，未知名才 warning）
check_annotation_names(cluster_map, preset)
# 映射证据校验：机器按 marker 打分复查每条编号映射（编号重排后贴错群在此暴露），
# 证据表写 annotation_evidence.csv 供人工对照复核
check_mapping_evidence(scobj, cluster_map, preset, step_dir)
# 双向检查：① config 写了对象中不存在的群编号（调低 resolution 后常见）
if (!all(names(cluster_map) %in% levels(Idents(scobj)))) {
  warning("config 中 cluster_map 的群编号与对象分群不完全一致，请检查")
}
# ② 对象中有未被映射的群——不拦截会静默变成数字 celltype 流向下游
unmapped <- setdiff(levels(Idents(scobj)), names(cluster_map))
if (length(unmapped) > 0) {
  warning("以下群未被 cluster_map 映射，将保持数字编号进入 celltype 列: ",
          paste(unmapped, collapse = ", "), "——请补齐 config 的 annotate$cluster_map")
}
scobj <- RenameIdents(scobj, cluster_map)
scobj[["celltype"]] <- Idents(scobj)              # 注释结果写入 metadata

save_fig(DimPlot(scobj, reduction = "umap", label = TRUE),
         "dimplot_annotated", type = "dimplot")
slim_save(scobj, "seurat_annotated.rds", step_dir)

# ---- 第四步（可选）：亚群细分（config 的 subcluster$enable 控制）----
if (isTRUE(cfg$subcluster$enable)) {
  log_step(paste("亚群细分 | 目标群:", cfg$subcluster$cluster,
                 "| resolution:", cfg$subcluster$resolution))

  # ⚠️ FindSubCluster 要求 Idents 是分群编号：注释后 Idents 是细胞类型名，必须先切回
  Idents(scobj) <- "seurat_clusters"

  # 对指定大群在 SNN 图上重新聚类（v5 中 SNN 图名为 RNA_snn）
  sub_name <- paste0("RNA_snn_res.", cfg$subcluster$resolution,
                     "_c", cfg$subcluster$cluster, "_sub")
  scobj <- FindSubCluster(scobj,
                          cluster = cfg$subcluster$cluster,
                          graph.name = "RNA_snn",
                          subcluster.name = sub_name,
                          resolution = cfg$subcluster$resolution)
  Idents(scobj) <- sub_name

  save_fig(DimPlot(scobj, reduction = "umap", label = TRUE),
           paste0("dimplot_subcluster_c", cfg$subcluster$cluster), type = "dimplot")

  # 细分完成切回细胞类型（metadata 的 celltype 列不受影响）
  Idents(scobj) <- "celltype"
  slim_save(scobj, "seurat_subclustered.rds", step_dir)
}

log_step("模块 6+7 完成 → 下一步: Rscript R/04_multi_group_plot.R")
