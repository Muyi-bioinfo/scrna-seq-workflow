###############################################################################
### 09_cell_communication.R — 高级模块：细胞通讯分析（CellChat）
###
### CellChat 基于配体-受体数据库推断细胞群间通讯强度与通路：
###   single 整体通讯网络 / multi 两组通讯差异比较（如 STIM vs CTRL）
### 输入：output/<batch>/03_annotate/seurat_annotated.rds（需有 celltype 列）
### 输出：output/<batch>/09_cell_communication/（网络图、通路气泡图、结果 rds）
###############################################################################

source("utils/utils.R")
cfg <- load_config()
step_dir <- setup_step("09_cell_communication")
log_step(paste("高级模块 | 细胞通讯分析 | 模式:", cfg$mode))
library(CellChat)

scobj <- read_prev("03_annotate", "seurat_annotated.rds")
Idents(scobj) <- "celltype"

# 配体-受体数据库（物种在 config 的 advanced$cellchat$species）
db_species <- if (!is.null(cfg$advanced$cellchat$species)) cfg$advanced$cellchat$species else "human"

# ---- 单个 CellChat 对象的完整流程（封装为函数，单/多样本共用）----
run_cellchat <- function(scobj_subset, species) {
  cellchat <- createCellChat(object = scobj_subset, group.by = "celltype", assay = "RNA")
  cellchat@DB <- if (species == "human") CellChatDB.human else CellChatDB.mouse

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # 通讯概率（triMean 截断均值法；类型 3 是 CellChat 推荐的默认值）
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)   # 网络"枢纽"细胞群
  cellchat
}

# ---- 单样本模式：整体通讯网络 ----
if (cfg$mode == "single") {
  cellchat <- run_cellchat(scobj, db_species)
  saveRDS(cellchat, file = file.path(step_dir, "cellchat.rds"))

  # ⚠️ netVisual_circle 用 base 图形直接绘制（Rscript 下返回的 recordedplot 为空），
  #   必须用 save_fig_draw 在设备内直接捕获
  save_fig_draw(netVisual_circle(cellchat@net$count,
                                 weight.scale = TRUE, label.edge = FALSE),
                "circle_communication_count")
  save_fig_draw(netVisual_circle(cellchat@net$weight,
                                 weight.scale = TRUE, label.edge = FALSE),
                "circle_communication_strength")

  # 细胞群信号角色热图：出向/入向信号强度
  save_fig(netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing"),
           "heatmap_signaling_role_outgoing", type = "heatmap")

  # 通路贡献气泡图：纵轴为信号通路
  save_fig(netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL,
                            remove.isolate = FALSE),
           "bubble_signaling_pathways", type = "dotplot")

# ---- 多样本模式：两组通讯差异比较 ----
} else {
  group_col <- cfg$multi_group$group_col
  groups <- sort(unique(scobj[[]][[group_col]]))
  if (length(groups) < 2) stop("通讯差异比较需要至少两个分组")

  cellchat_list <- list()
  for (g in groups) {
    message("-- 分析分组: ", g)
    # 按细胞名筛选（避免 NSE 表达式在不同环境下的兼容性问题）
    cells_g <- colnames(scobj)[scobj[[]][[group_col]] == g]
    cellchat_list[[g]] <- run_cellchat(subset(scobj, cells = cells_g), db_species)
  }
  names(cellchat_list) <- groups
  saveRDS(cellchat_list, file = file.path(step_dir, "cellchat_list.rds"))

  # 合并两组对象做差异比较（CellChat 经典用法，如 IFN-β 刺激 PBMC）
  object <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

  # 1. 通讯数量与强度的总体比较
  save_fig(compareInteractions(object, show.legend = FALSE, group = 1:2),
           "compare_interaction_counts", type = "barplot")

  # 2. 两组通讯网络图并列对比（同样 base 图形，需 save_fig_draw 捕获）
  save_fig_draw(netVisual_diffInteraction(object, weight.scale = TRUE),
                "diff_network_all")

  # 3. 通讯通路强度差异热图
  save_fig(rankNet(object, mode = "comparison", stacked = FALSE, do.stat = TRUE),
           "ranknet_pathway_comparison", type = "heatmap")
}

log_step("高级模块 | 细胞通讯分析完成 → 结果见 output/09_cell_communication/")
