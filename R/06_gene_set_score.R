###############################################################################
### 06_gene_set_score.R — 高级模块：基因集打分（AddModuleScore 通路活性推断）
###
### 输入：output/<batch>/03_annotate/seurat_annotated.rds + config 的 advanced$gene_set_gmt
### 输出：output/<batch>/06_gene_set_score/seurat_scored.rds + 通路打分 FeaturePlot
###############################################################################

source("utils/utils.R")
cfg <- load_config()
step_dir <- setup_step("06_gene_set_score")
log_step(paste("高级模块 B | 基因集打分 | 基因集文件:", cfg$advanced$gene_set_gmt))
library(clusterProfiler)
set.seed(123)   # ⚠️ AddModuleScore 随机抽样对照基因集，设种子保证打分可复现

scobj <- read_prev("03_annotate", "seurat_annotated.rds")

# gmt 文件检查：缺失时给出配置层面的提示，而非文件系统报错
gmt_file <- cfg$advanced$gene_set_gmt
if (is.null(gmt_file) || !file.exists(gmt_file)) {
  stop("基因集文件不存在: ", gmt_file, "，请检查 config 的 advanced$gene_set_gmt")
}

# gmt → 基因集列表（每个通路一个基因向量）→ AddModuleScore 打分
genesets <- read.gmt(gmt_file)
signatures <- split(genesets$gene, genesets$term)
message("-- 基因集数量: ", length(signatures))
scobj <- AddModuleScore(scobj, features = signatures, name = "geneset")

# 把 "geneset1/2/3..." 列名替换为通路名，便于作图（全部走 [[]] API）
old_names <- grep("geneset\\d", colnames(scobj[[]]), value = TRUE)
for (i in seq_along(old_names)) {
  scobj[[names(signatures)[i]]] <- scobj[[old_names[i]]]
  scobj[[old_names[i]]] <- NULL
}

# 展示通路：config 的 advanced$gene_set_show 指定（优先），留空取 gmt 前两条；
# 每个通路单独一张 FeaturePlot（多通路同图在 8×7 里每个面板过窄，标签重叠）
show_pathways <- cfg$advanced$gene_set_show
if (is.null(show_pathways) || length(show_pathways) == 0) {
  show_pathways <- names(signatures)[1:min(2, length(signatures))]
}
show_pathways <- intersect(show_pathways, names(signatures))
if (length(show_pathways) == 0) {
  stop("advanced$gene_set_show 指定的通路均不在 gmt 文件中，请检查通路名")
}
for (pw in show_pathways) {
  fname <- gsub("[^A-Za-z0-9]+", "_", pw)   # 通路名含空格/符号时仍可安全作文件名
  # ⚠️ 显式 reduction="umap"（umap_naive 抢占默认降维的坑，见 README 踩坑表）
  save_fig(FeaturePlot(scobj, features = pw, label = TRUE, repel = TRUE,
                       reduction = "umap"),
           paste0("featureplot_", fname), type = "featureplot")
}

slim_save(scobj, "seurat_scored.rds", step_dir)
log_step("高级模块 B 完成 —— 全部高级模块结束")
