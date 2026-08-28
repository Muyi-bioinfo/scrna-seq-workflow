###############################################################################
### 05_diff_gsea.R — 模块 9：差异分析 + GSEA 富集分析
###
### 两种策略：① 群水平 find marker → compareCluster 批量 KEGG 富集
###           ② 单群刺激前后 FindMarkers → geneList 排序 → GSEA
### 输入：output/<batch>/03_annotate/seurat_annotated.rds；输出：差异结果 rds + 富集图
###############################################################################

source("utils/utils.R")
cfg <- load_config()
# 05 是多样本模块：单样本配置没有 multi_group / diff_gsea 段，
# 直接跑会静默退化为"无分组分析"，产出看似成功的错误结果——先显式拦截
if (is.null(cfg$multi_group) || is.null(cfg$diff_gsea)) {
  stop("05 需要多样本配置（config.multi.yaml 的 multi_group 与 diff_gsea 段），单样本流程请跳过本步")
}
step_dir <- setup_step("05_diff_gsea")
log_step("模块 9 | 差异分析 + GSEA")
set.seed(123)   # ⚠️ GSEA 依赖随机置换计算 p 值，设种子保证可复现
library(ggplot2)   # theme()/facet_grid() 等绘图函数

scobj <- read_prev("03_annotate", "seurat_annotated.rds")
group_col <- cfg$multi_group$group_col
# celltype.stim 的编码分隔符：细胞类型名可能含下划线（自定义注释如 T_reg），
# 用 "_" 编码再按 "_" 解码会错位，故用类型名/组名里不会出现的 "^^"
STIM_SEP <- "^^"
STIM_SEP_RE <- "\\^\\^"   # separate() 的 sep 按正则解析，需转义

# ---- 策略一：群水平富集（compareCluster）----
scobj[["celltype.stim"]] <- paste(scobj$celltype, scobj[[]][[group_col]], sep = STIM_SEP)
Idents(scobj) <- "celltype.stim"

sce_markers <- FindAllMarkers(object = scobj, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25)   # thresh.use 是 v4 废弃参数名
saveRDS(sce_markers, file = file.path(step_dir, "split_markers.rds"))

library(clusterProfiler)
library(enrichplot)
library(tidyr)

# 物种来自 config 的 diff_gsea$organism（默认 human）：决定 OrgDb 与 KEGG 的物种背景
organism <- if (!is.null(cfg$diff_gsea$organism)) cfg$diff_gsea$organism else "human"
orgdb <- if (organism == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
kegg_org <- if (organism == "human") "hsa" else "mmu"
library(orgdb, character.only = TRUE)

markers <- sce_markers %>%
  filter(p_val_adj < 0.001) %>%
  separate(cluster, into = c("celltype", "group"), sep = STIM_SEP_RE, remove = FALSE)
gid <- bitr(unique(markers$gene), "SYMBOL", "ENTREZID", OrgDb = orgdb)
colnames(gid)[1] <- "gene"
markers <- merge(markers, gid, by = "gene")
# ⚠️ bitr 存在 1:many 映射（一个 SYMBOL 对应多个 ENTREZID）：去重，避免同一基因多行重复计数
# ⚠️ distinct() 不带 .keep_all 时只保留列出的列——必须 .keep_all = TRUE，
#    否则 celltype/group 被丢弃，后续 filter 和 compareCluster 直接报错
markers <- distinct(markers, gene, cluster, ENTREZID, .keep_all = TRUE)

# 剔除低价值细胞类型（如红细胞/巨核细胞），列表在 config 的 diff_gsea$exclude_celltypes；
# 注意：列表必须与注释名完全一致（死字符串会静默失效）
exclude_ct <- cfg$diff_gsea$exclude_celltypes
if (!is.null(exclude_ct) && length(exclude_ct) > 0) {
  markers <- markers %>% filter(!celltype %in% exclude_ct)
  message("-- 已从富集剔除细胞类型: ", paste(exclude_ct, collapse = ", "))
}

x <- compareCluster(ENTREZID ~ celltype + group, data = markers,
                    fun = "enrichKEGG", organism = kegg_org)
p <- dotplot(x, label_format = 60, x = "group") +
  facet_grid(~celltype) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p, "compareCluster_dotplot", type = "dotplot")

# ---- 策略二：单群刺激前后 GSEA ----
ident_compare <- cfg$diff_gsea$ident_compare
if (!is.null(ident_compare) && nzchar(ident_compare)) {
  # 比较方向：处理组在前（ident.1），对照在后（ident.2）
  # 优先用 config 的 diff_gsea$ident1_group 指定，否则按字母倒序（STIM 在 CTRL 前）
  groups <- sort(unique(scobj[[]][[group_col]]), decreasing = TRUE)
  if (!is.null(cfg$diff_gsea$ident1_group) &&
      cfg$diff_gsea$ident1_group %in% groups) {
    groups <- c(cfg$diff_gsea$ident1_group, setdiff(groups, cfg$diff_gsea$ident1_group))
  }
  if (length(groups) >= 2) {
    ident1 <- paste0(ident_compare, STIM_SEP, groups[1])
    ident2 <- paste0(ident_compare, STIM_SEP, groups[2])
    # 三组以上设计：策略二只做一次两两比较，多余分组静默丢弃会误导读图——显式提醒
    if (length(groups) > 2) {
      warning("分组多于 2 个（", paste(groups, collapse = ", "), "）：策略二仅比较 ",
              ident1, " vs ", ident2, "，其余分组被忽略；多组两两比较请改用策略一或拆分 config 重跑")
    }
    message("-- 比较: ", ident1, " vs ", ident2)

    # 差异基因：logfc.threshold = 0 保留全部基因用于 GSEA 排序
    diff_res <- FindMarkers(scobj, ident.1 = ident1, ident.2 = ident2,
                            logfc.threshold = 0)
    saveRDS(diff_res, file = file.path(step_dir, "diff_stim_vs_ctrl.rds"))

    # geneList 三部曲：取 logFC → 命名 → 降序排序（排序很重要！）
    geneList <- diff_res$avg_log2FC
    names(geneList) <- rownames(diff_res)
    geneList <- sort(geneList, decreasing = TRUE)

    gmt_files <- list(kegg = cfg$diff_gsea$gmt_kegg, hallmark = cfg$diff_gsea$gmt_hallmark)
    for (nm in names(gmt_files)) {
      gmt <- gmt_files[[nm]]
      if (is.null(gmt) || !file.exists(gmt)) {
        warning("基因集文件不存在，跳过 ", nm, ": ", gmt)
        next
      }
      y <- GSEA(geneList, TERM2GENE = read.gmt(gmt))
      # 零结果保护：gmt 与数据基因名不匹配（如鼠基因名跑人 gmt）时 GSEA 返回空表，
      # 后续 dotplot 会报晦涩错误并中断脚本（主结果已产出但 .done 不落盘，重跑代价大）
      y_df <- as.data.frame(y)
      if (nrow(y_df) == 0) {
        warning(nm, " 基因集 GSEA 无任何富集结果（检查 gmt 物种与格式），跳过作图: ", gmt)
        next
      }
      # ⚠️ 正负两个方向都存在时才按 .sign 拆分 facet（单向结果 facet 会报错中断整个脚本）
      if (length(unique(y_df$.sign)) >= 2) {
        save_fig(dotplot(y, showCategory = 12, split = ".sign") + facet_grid(~.sign),
                 paste0("gsea_", nm, "_dotplot"), type = "gsea")
      } else {
        save_fig(dotplot(y, showCategory = 12),
                 paste0("gsea_", nm, "_dotplot"), type = "gsea")
      }
    }
  } else {
    warning("分组数不足 2，跳过刺激-对照 GSEA")
  }
}

log_step("模块 9 完成 —— 全流程结束！结果见 output/ 各步骤目录")
