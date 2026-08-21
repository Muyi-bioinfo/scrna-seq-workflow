###############################################################################
### 04_multi_group_plot.R — 模块 8：多组比较可视化（比例图 + 分组 DotPlot + 热图）
###
### 输入：output/<batch>/03_annotate/seurat_annotated.rds（需 group 和 celltype 列）
### 输出：output/<batch>/04_multi_group_plot/figures/（只出图，无 rds）
###############################################################################

source("utils/utils.R")
cfg <- load_config()
# 04 是多样本模块：单样本配置没有 multi_group 段，直接跑会得到晦涩报错，先显式拦截
if (is.null(cfg$multi_group) || is.null(cfg$multi_group$group_col) ||
    !nzchar(cfg$multi_group$group_col)) {
  stop("04 需要多样本配置（config.multi.yaml 的 multi_group$group_col），单样本流程请跳过本步")
}
step_dir <- setup_step("04_multi_group_plot")
log_step(paste("模块 8 | 多组比较可视化 | 分组列:", cfg$multi_group$group_col))
library(ggplot2)   # theme() 等绘图函数（utils.R 只加载了 Seurat/dplyr，ggplot2 需显式加载）

scobj <- read_prev("03_annotate", "seurat_annotated.rds")
group_col <- cfg$multi_group$group_col
preset <- load_annotation_preset(cfg$annotate, scobj)   # DotPlot 用预设 tier1 marker

# ---- 1. 各组细胞类型比例堆叠图 ----
data <- as.data.frame(table(scobj[[]][[group_col]], scobj$celltype))
colnames(data) <- c("group", "CellType", "Freq")
df <- data %>%
  group_by(group) %>%
  mutate(Percent = Freq / sum(Freq)) %>%
  ungroup() %>%
  as.data.frame()

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = Percent, fill = CellType)) +
  ggplot2::geom_bar(position = "fill", stat = "identity", color = "white", width = 0.95) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::theme_classic()
save_fig(p1, "barplot_cell_proportion", type = "barplot")

# ---- 2. 分组 DotPlot：split.by 在同一图中比较两组 marker 表达 ----
markers <- preset$flat_tier1
p2 <- DotPlot(scobj, features = markers,
              cols = c("blue", "red"),
              dot.scale = 6,
              split.by = group_col) +
  RotatedAxis() +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 8))
save_fig(p2, "dotplot_grouped", type = "dotplot")

# ---- 3. marker 热图（需要 scale.data；downsample 抽样控制内存）----
# prev_path 只取路径不读入，marker 表缺失时跳过而非报错
markers_file <- prev_path("03_annotate", "all_markers_table.rds")
if (file.exists(markers_file)) {
  top5 <- readRDS(markers_file) %>%
    filter(pct.1 > cfg$annotate$marker_pct) %>%   # 预筛选阈值与 03 同源（annotate$marker_pct）
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice(1:5) %>%
    ungroup() %>%
    pull(gene) %>%
    unique()
  scobj <- ScaleData(scobj, features = top5)
  p3 <- DoHeatmap(subset(scobj, downsample = 50), features = top5, size = 3) +
    theme(axis.text.y = element_text(size = 8))
  save_fig(p3, "heatmap_top5_markers", type = "heatmap")
  scobj[["RNA"]]$scale.data <- matrix()   # 用完即删，控制内存
} else {
  warning("未找到 marker 表（需先运行 03_annotate），跳过热图")
}

log_step("模块 8 完成 → 下一步: Rscript R/05_diff_gsea.R")
