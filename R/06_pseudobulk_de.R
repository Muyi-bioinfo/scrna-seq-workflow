###############################################################################
### 06_pseudobulk_de.R — 模块 9+：pseudobulk 确认性差异分析（sample-aware DESeq2）
###
### 05 以细胞为统计单元（伪重复，p 值虚高）；本模块按 样本×细胞类型 聚合原始 counts，
### 逐类型独立建 DESeq2：同一样本出现在 ≥2 组 → 配对 design ~ sample_id + condition，
### 否则 ~ condition；batch 可用且未与 condition 混淆时插到最前。
### ⚠️ 只吃 RNA assay 的 raw counts（DESeq2 负二项模型不吃标准化/SCT/Harmony 产物）。
###
### 输入：output/<batch>/03_annotate/seurat_annotated.rds（sample_id 由 01 写入）
### 输出：output/<batch>/06_pseudobulk_de/ —— by_celltype/<type>/（全基因结果/显著表）
###      + figures/<type>/（volcano/MA/PCA/heatmap）；顶层 summary_degs.csv 含单细胞口径对照
###############################################################################

source("utils/utils.R")
cfg <- load_config()
# 06 是多样本模块：config 没有 pseudobulk 段（单样本无生物学重复）时直接拦截，
# 不静默跑出无 replicate 的假结果——与 05 的 multi_group 拦截同一模式
if (is.null(cfg$pseudobulk)) {
  stop("06 需要 config 的 pseudobulk 段（见 config.multi.yaml）；单样本流程请跳过本步")
}
pb <- cfg$pseudobulk
step_dir <- setup_step("06_pseudobulk_de")
log_step("模块 9+ | pseudobulk 确认性差异分析（DESeq2, sample-aware）")
set.seed(123)
library(ggplot2)   # 火山图/PCA/对比条形图

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("未安装 DESeq2（conda 安装: mamba install -n scrna -c bioconda -c conda-forge ",
       "bioconductor-deseq2 bioconductor-apeglm）")
}

# ---- config 读取（缺省值兜底，与 05 的 organism 处理同一模式）----
get_cfg <- function(x, default) if (is.null(x)) default else x
sample_col     <- get_cfg(pb$sample_col, "sample_id")
condition_col  <- get_cfg(pb$condition_col, "group")
group_col      <- if (!is.null(cfg$multi_group$group_col)) cfg$multi_group$group_col else condition_col
batch_col      <- get_cfg(pb$batch_col, NULL)     # 可选：null/缺失/单一取值时退出 design
celltype_col   <- get_cfg(pb$celltype_col, "celltype")
min_cells      <- get_cfg(pb$min_cells, 20)
min_samp_group <- get_cfg(pb$min_samples_per_group, 2)
padj_cut       <- get_cfg(pb$padj_cutoff, 0.05)
lfc_cut        <- get_cfg(pb$log2fc_cutoff, 1)
test_grp       <- pb$test_group
ref_grp        <- pb$reference_group
lfc_shrink     <- isTRUE(get_cfg(pb$lfc_shrink, TRUE))
run_gsea       <- isTRUE(get_cfg(pb$run_gsea, FALSE))
target_ct      <- get_cfg(pb$target_celltypes, NULL)   # 白名单：NULL/[] = 全部分析
exclude_ct     <- get_cfg(pb$exclude_celltypes, NULL)
if (is.null(test_grp) || is.null(ref_grp)) {
  stop("config 需显式指定 pseudobulk$test_group 与 pseudobulk$reference_group")
}
if (identical(test_grp, ref_grp)) stop("test_group 与 reference_group 不能相同")

scobj <- read_prev("03_annotate", "seurat_annotated.rds")
meta <- scobj[[]]

# ---- 1. 必要列校验（Seurat 5.5 读 metadata 列用 [[]][[列名]]，见 CLAUDE.md）----
missing_cols <- setdiff(c(condition_col, celltype_col), colnames(meta))
if (length(missing_cols) > 0) {
  stop("Seurat metadata 缺少必要列: ", paste(missing_cols, collapse = ", "),
       "（检查 config 的 pseudobulk$condition_col / celltype_col）")
}

# ---- 2. 样本列校验（01 已写入：常规场景 = 样本表每行；合并矩阵场景 =
#        multi$cell_metadata 按 (group, barcode) join 供者，见 config.multi.yaml）----
if (!sample_col %in% colnames(meta)) {
  stop("metadata 缺少样本列 \"", sample_col, "\"——01 的多样本分支会写它：常规场景取样本表每行，",
       "合并矩阵场景用 multi$cell_metadata 按 (group, barcode) join 供者")
}

# ---- 3. 按 sample × condition × celltype 聚合原始 counts ----
cell_df <- data.frame(
  cell      = colnames(scobj),
  sample_id = as.character(scobj[[sample_col]][[1]]),
  condition = as.character(scobj[[]][[condition_col]]),
  celltype  = as.character(scobj[[]][[celltype_col]]),
  batch     = if (!is.null(batch_col) && batch_col %in% colnames(meta)) {
                as.character(scobj[[]][[batch_col]])
              } else NA_character_,
  stringsAsFactors = FALSE)

n_before <- nrow(cell_df)
cell_df <- cell_df[!is.na(cell_df$sample_id), , drop = FALSE]
if (!is.null(exclude_ct) && length(exclude_ct) > 0) {
  drop_ct <- intersect(exclude_ct, unique(cell_df$celltype))
  if (length(drop_ct) > 0) {
    cell_df <- cell_df[!cell_df$celltype %in% drop_ct, ]
    message("-- 已按 config 剔除细胞类型: ", paste(drop_ct, collapse = ", "))
  } else {
    warning("exclude_celltypes 与数据中的注释名无一匹配（死字符串？）: ",
            paste(exclude_ct, collapse = ", "))
  }
}
# 白名单：非空时只保留目标类型（与 exclude 叠加 = 先白后黑；留空 = 全部）
if (!is.null(target_ct) && length(target_ct) > 0) {
  miss <- setdiff(target_ct, unique(cell_df$celltype))
  if (length(miss) > 0) {
    warning("target_celltypes 与数据中的注释名不匹配（死字符串？）: ",
            paste(miss, collapse = ", "))
  }
  cell_df <- cell_df[cell_df$celltype %in% target_ct, , drop = FALSE]
  message("-- 白名单生效，仅分析: ", paste(sort(unique(cell_df$celltype)), collapse = ", "))
}
message("-- 聚合用细胞: ", nrow(cell_df), "/", n_before, "（剔除无供者信息、黑/白名单过滤）")

# 聚合走主流接口（Seurat de_vignette 的 pseudobulk 做法）：AggregateExpression 求和，
# return.seurat = TRUE 时 group.by 列自动进 meta.data——不解析细胞名、无矩阵方向问题
scobj_sub <- subset(scobj, cells = cell_df$cell)   # 无供者信息的细胞先剔除
pb_seu <- AggregateExpression(scobj_sub, assays = "RNA", return.seurat = TRUE,
                              group.by = c(sample_col, condition_col, celltype_col))
counts_pb <- GetAssayData(pb_seu, layer = "counts")   # 基因 × pseudobulk（DESeq2 方向）

# pb_meta 直接取 meta.data 的分组列。⚠️ Seurat 会给数字开头的分组值加 "g" 前缀
#（合法细胞名，与 AverageExpression 的 g0/g1 同款行为）——剥掉；剥后值必须在
# 细胞层集合里才生效（防误伤以 g 开头的真实 id）
strip_g <- function(v, valid) {
  s <- sub("^g", "", v)
  ifelse(s %in% valid, s, v)
}
pb_meta <- data.frame(
  pseudobulk_id = colnames(pb_seu),
  sample_id     = strip_g(as.character(pb_seu[[sample_col]][[1]]), unique(cell_df$sample_id)),
  condition     = strip_g(as.character(pb_seu[[]][[condition_col]]), unique(cell_df$condition)),
  celltype      = strip_g(as.character(pb_seu[[]][[celltype_col]]), unique(cell_df$celltype)),
  stringsAsFactors = FALSE)
# 每格细胞数/批次按三元组计数对齐；键用控制字符 \01 连接（类型名/组名里不会出现）
key_cells <- paste(cell_df$sample_id, cell_df$condition, cell_df$celltype, sep = "\01")
key_pb <- paste(pb_meta$sample_id, pb_meta$condition, pb_meta$celltype, sep = "\01")
# batch 聚合：同一 pseudobulk 跨批次时 batch 因素不可用（记 NA）
batch_pb <- tapply(cell_df$batch, key_cells,
                   function(x) if (all(!is.na(x)) && length(unique(x)) == 1) unique(x)[1] else NA_character_)
pb_meta$batch <- unname(batch_pb[key_pb])
pb_meta$n_cells <- as.integer(table(key_cells)[key_pb])
if (anyNA(pb_meta$n_cells)) {
  stop("pseudobulk 与细胞层三元组对不上（AggregateExpression 分组名改写？），",
       "请检查 sample/condition/celltype 列的取值")
}
rownames(pb_meta) <- pb_meta$pseudobulk_id

# ---- 4. pseudobulk 水平过滤与整理 ----
pb_meta$kept <- pb_meta$n_cells >= min_cells
write.table(pb_meta, file = file.path(step_dir, "pseudobulk_cell_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
n_drop <- sum(!pb_meta$kept)
if (n_drop > 0) {
  message("-- min_cells 过滤: 丢弃 ", n_drop, " 个不足 ", min_cells,
          " 细胞的 pseudobulk（明细见 pseudobulk_cell_counts.tsv）")
}
pb_meta <- pb_meta[pb_meta$kept, , drop = FALSE]
pb_meta$kept <- NULL
if (nrow(pb_meta) == 0) {
  stop("min_cells = ", min_cells, " 过滤后没有任何 pseudobulk 剩下，请放宽阈值")
}

# 未知组别显式拦截（factor 化会把它们变 NA，静默丢数据）
unknown_grp <- setdiff(unique(pb_meta$condition), c(ref_grp, test_grp))
if (length(unknown_grp) > 0) {
  warning("组别 ", paste(unknown_grp, collapse = ", "),
          " 不在 test/reference 中，相关 pseudobulk 已排除")
  pb_meta <- pb_meta[pb_meta$condition %in% c(ref_grp, test_grp), , drop = FALSE]
}
# colData 列名固定（sample_id/condition/batch），design 公式与 config 原列名解耦
pb_meta$condition <- factor(pb_meta$condition, levels = c(ref_grp, test_grp))
pb_meta$sample_id <- factor(pb_meta$sample_id)

# batch 是否进 design：列存在、≥2 个非 NA 取值、且未与 condition 混淆
use_batch <- FALSE
if (!is.null(batch_col) && batch_col %in% colnames(meta)) {
  bt <- table(pb_meta$batch, pb_meta$condition)
  if (length(unique(na.omit(pb_meta$batch))) < 2) {
    message("-- batch 单一取值（或不完整），不进 design（列保留在 metadata 备查）")
  } else if (any(rowSums(bt > 0) == 1)) {
    warning("batch 与 condition 完全混淆（每个 batch 只出现在一个组），",
            "同时进 design 会不可估——batch 退出 design")
  } else {
    use_batch <- TRUE
    message("-- batch 因素进入 design: ~ batch + ...")
  }
}
write.table(pb_meta, file = file.path(step_dir, "pseudobulk_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 图形样式（期刊风）----
# 调色板与 docs/assets/make_pipeline_figure.py 同一体系（经 CVD/对比度校验）：
#   上/下调用发散对（暖红/冷蓝），不显著用中性灰退到背景；
#   条件色（PCA/椭圆）用分类槽 1/2（冷蓝 = 对照，暖橙 = 刺激，与 volcan
# 语义同向：暖 = 刺激方向）。文字一律墨色，不借颜色传义
PAL_UP <- "#e34948"; PAL_DOWN <- "#2a78d6"; PAL_NS <- "#c3c2b7"
PAL_CTRL <- "#2a78d6"; PAL_STIM <- "#eb6834"
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
HAIR <- "#e1e0d9"; BASELINE <- "#c3c2b7"; SURFACE <- "#fcfcfb"

#' 统一主题：浅色面板 + 发丝网格 + 顶部图例，去默认灰底
theme_pb <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = HAIR, linewidth = 0.3),
      axis.line = element_line(color = BASELINE, linewidth = 0.4),
      axis.text = element_text(color = INK),
      axis.title = element_text(color = INK2),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.background = element_rect(fill = SURFACE, color = NA),
      legend.position = "top",
      legend.title = element_blank(),
      legend.key.size = grid::unit(10, "pt")
    )
}

#' 样本层 PCA（期刊风）：条件着色 + 90% 正态椭圆 + 供者标签（ggrepel）。
#' 小分组（<3 点）自动跳过椭圆（stat_ellipse 会在 2 点时报错）
plot_pb_pca <- function(df, pve, title) {
  # 图例不带样本数：单臂供者剔除保证进 PCA 的样本必然配对等数（每个供者双条件齐全）
  cond_levels <- c(ref_grp, test_grp)
  p <- ggplot(df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = setNames(c(PAL_CTRL, PAL_STIM), cond_levels)) +
    labs(title = title,
         x = sprintf("PC1 · %.1f%% variance", 100 * pve[1]),
         y = sprintf("PC2 · %.1f%% variance", 100 * pve[2])) +
    theme_pb()
  if (min(table(df$condition)) >= 3) {
    p <- p + stat_ellipse(type = "norm", level = 0.9,
                          linetype = "dashed", linewidth = 0.5)
  }
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(aes(label = sample_id), size = 2.8, color = INK,
                                      max.overlaps = 20, segment.color = MUTED,
                                      segment.size = 0.3, box.padding = 0.3)
  }
  p
}

#' 取 top 高变 pseudobulk 做 PCA，返回得分与解释度（fallback PCA 用）
pca_scores <- function(mat, ntop = 2000) {
  mat <- log2(as.matrix(mat) + 1)
  v <- apply(mat, 1, var)
  mat <- mat[order(v, decreasing = TRUE), , drop = FALSE]
  mat <- mat[seq_len(min(ntop, nrow(mat))), , drop = FALSE]
  pc <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  list(df = as.data.frame(pc$x), pve = summary(pc)$importance[2, ])
}

# ---- 5. GSEA 预备（可选；默认关）----
if (run_gsea) {
  if (is.null(cfg$diff_gsea)) {
    warning("run_gsea=true 但 config 缺少 diff_gsea 段（gmt 路径），跳过 pseudobulk GSEA")
    run_gsea <- FALSE
  } else if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    warning("clusterProfiler 未安装，跳过 pseudobulk GSEA")
    run_gsea <- FALSE
  } else {
    library(clusterProfiler)   # GSEA()/read.gmt()
    library(enrichplot)        # dotplot()
  }
}

# ---- 6. 逐细胞类型 DESeq2 ----
# 逐类型列表 = 注释全集（含最终被跳过的类型：summary 里给 skip_reason，不留黑洞）
cts_all <- sort(unique(as.character(scobj[[]][[celltype_col]])))

# 单细胞对照基准（同对比口径）：celltype.stim 组合身份与 05 同一编码约定（^^），
# 逐类型 FindMarkers(STIM vs CTRL)，与 DESeq2 取相同基因全集、相同 BH 校正比较显著数
scobj[["celltype.stim"]] <- paste(scobj$celltype, scobj[[]][[group_col]], sep = "^^")
Idents(scobj) <- "celltype.stim"

#' 细胞类型名 → 文件名安全目录名（表内仍保留原名）
sanitize_ct <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
ct_dirs <- sanitize_ct(cts_all)
names(ct_dirs) <- cts_all   # 目录名映射覆盖注释全集（跳过的类型只是不建目录）
if (any(duplicated(ct_dirs))) stop("细胞类型名安全化后重名，无法建目录: ",
                                   paste(unique(ct_dirs[duplicated(ct_dirs)]), collapse = ", "))

summary_rows <- list()
# 实际分析集合：cell_df 已含黑/白名单过滤（用于 summary 的 skip_reason 区分）
analyzed_ct <- unique(cell_df$celltype)
for (ct in cts_all) {
  ct_dir <- file.path(step_dir, "by_celltype", ct_dirs[[ct]])   # 表（数据）
  sub_meta <- droplevels(pb_meta[pb_meta$celltype == ct, , drop = FALSE])
  row <- data.frame(celltype = ct, skip_reason = NA_character_,
                    design = NA_character_,
                    n_pseudobulk_ref = sum(sub_meta$condition == ref_grp),
                    n_pseudobulk_test = sum(sub_meta$condition == test_grp),
                    n_tested_genes = NA_integer_,
                    n_sig_padj_lfc = NA_integer_, n_up = NA_integer_, n_down = NA_integer_,
                    n_sig_up = NA_integer_,
                    n_sc_up_bh = NA_integer_, n_pb_up_bh = NA_integer_)
  # 黑/白名单外的类型与聚合后无样本的类型——明确记入 summary 而非无声消失
  if (!ct %in% analyzed_ct) {
    row$skip_reason <- if (!is.null(exclude_ct) && ct %in% exclude_ct) "config 剔除（exclude_celltypes）"
                       else "非目标类型（config 白名单）"
    summary_rows[[ct]] <- row
    next
  }
  if (!ct %in% pb_meta$celltype) {
    row$skip_reason <- "min_cells 过滤后无 pseudobulk 样本"
    summary_rows[[ct]] <- row
    next
  }
  message("── 细胞类型: ", ct, "（", row$n_pseudobulk_ref, " ", ref_grp, " / ",
          row$n_pseudobulk_test, " ", test_grp, "）")

  # design 构建：配对（同一供者出现在 ≥2 组）→ 供者进 design 吸收个体差异
  sc_tab <- table(sub_meta$sample_id, sub_meta$condition)
  paired <- any(rowSums(sc_tab > 0) > 1)
  if (paired) {
    # 单臂供者剔除：只有一组有合格 pseudobulk 的供者在配对 design 里只添参数、
    # 不添配对信息，白白吃掉残差自由度（小类型常见，剔掉常能救活整型）
    single_arm <- names(which(rowSums(sc_tab > 0) < 2))
    if (length(single_arm) > 0) {
      message("-- 剔除单臂供者 ", length(single_arm), " 个（无配对信息）: ",
              paste(single_arm, collapse = ", "))
      sub_meta <- droplevels(sub_meta[!sub_meta$sample_id %in% single_arm, , drop = FALSE])
      # summary 记录最终建模用的样本数（剔除后必然配对等数，此前记的是剔除前数字）
      row$n_pseudobulk_ref <- sum(sub_meta$condition == ref_grp)
      row$n_pseudobulk_test <- sum(sub_meta$condition == test_grp)
    }
  }

  # 每组最少样本数守卫（单臂剔除后复核）
  if (any(c(sum(sub_meta$condition == ref_grp), sum(sub_meta$condition == test_grp)) < min_samp_group)) {
    warning("类型 ", ct, " 每组 pseudobulk 样本不足 ", min_samp_group,
            "（", ref_grp, ":", row$n_pseudobulk_ref, " / ", test_grp, ":",
            row$n_pseudobulk_test, "），跳过")
    row$skip_reason <- "min_samples_per_group 不足"
    summary_rows[[ct]] <- row
    next
  }
  terms <- c(if (use_batch) "batch", if (paired) "sample_id", "condition")
  fml <- as.formula(paste("~", paste(terms, collapse = " + ")))
  row$design <- paste(deparse(fml), collapse = "")
  # 残差自由度守卫：样本数 - 参数数 < 1 时 DESeq2 拟合退化（结果几乎全 NA）
  mm <- model.matrix(fml, data = sub_meta)
  if (nrow(mm) - ncol(mm) < 1) {
    warning("类型 ", ct, " 自由度不足（", nrow(mm), " 样本 / ", ncol(mm),
            " 参数，design: ", row$design, "），跳过——配对设计建议每组 ≥3 个样本")
    row$skip_reason <- "残差自由度不足"
    summary_rows[[ct]] <- row
    next
  }

  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
  # 图统一进 figures/<type>/（管线约定：所有 PDF 在 figures/ 下），表留在 by_celltype/<type>/
  options(fig_outdir = file.path(step_dir, "figures", ct_dirs[[ct]]))
  cnt <- as.matrix(counts_pb[, sub_meta$pseudobulk_id, drop = FALSE])   # 基因×样本
  storage.mode(cnt) <- "integer"   # AggregateExpression 产出 double；DESeq2 要整数
  # 拟合失败保护：罕见数值问题（如全零/极低计数矩阵）只 warning 跳过本类型，
  # 不中断整步（参考 hossainlab/sc-workflow 的 tryCatch 模式）
  dds <- tryCatch({
    dds_tmp <- DESeq2::DESeqDataSetFromMatrix(countData = cnt, colData = sub_meta, design = fml)
    DESeq2::DESeq(dds_tmp, quiet = TRUE)
  }, error = function(e) {
    warning("类型 ", ct, " 的 DESeq2 拟合失败，跳过: ", conditionMessage(e))
    NULL
  })
  if (is.null(dds)) {
    row$skip_reason <- "DESeq2 拟合失败"
    summary_rows[[ct]] <- row
    next
  }
  res <- DESeq2::results(dds, contrast = c("condition", test_grp, ref_grp), alpha = padj_cut)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)

  # logFC 收缩：只改善低表达基因的效应量展示（volcano/排序讨论），不影响 padj/stat；
  # apeglm 优先（需 coef），否则 ashr（支持 contrast），都没有则不收缩并提示
  lfc_col <- "log2FoldChange"
  if (lfc_shrink) {
    tryCatch({
      coef_nm <- paste0("condition_", test_grp, "_vs_", ref_grp)
      res_s <- NULL
      if (coef_nm %in% DESeq2::resultsNames(dds) && requireNamespace("apeglm", quietly = TRUE)) {
        res_s <- DESeq2::lfcShrink(dds, coef = coef_nm, type = "apeglm", quiet = TRUE)
      } else if (requireNamespace("ashr", quietly = TRUE)) {
        res_s <- DESeq2::lfcShrink(dds, contrast = c("condition", test_grp, ref_grp),
                                   type = "ashr", quiet = TRUE)
      } else {
        message("-- apeglm/ashr 均未安装，log2FC 不收缩（不影响 padj/stat）")
      }
      if (!is.null(res_s)) {
        res_df$log2FC_shrunk <- res_s[match(res_df$gene, rownames(res_s)), "log2FoldChange"]
        lfc_col <- "log2FC_shrunk"
      }
    }, error = function(e) message("-- lfcShrink 失败（不影响主结果）: ", conditionMessage(e)))
  }

  # 显著判定：padj 阈值 + 展示列 log2FC 阈值（收缩列存在时用收缩值，更保守诚实）
  res_df$sig <- !is.na(res_df$padj) & res_df$padj < padj_cut &
    abs(res_df[[lfc_col]]) >= lfc_cut
  # 同 padj-only 口径（不卡 log2FC，STIM 上调）——与 05 对比用
  n_padj_only_up <- sum(!is.na(res_df$padj) & res_df$padj < padj_cut & res_df[[lfc_col]] > 0)

  # ---- 逐类型产出 ----
  write.table(res_df, file = file.path(ct_dir, "deseq2_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # 保留全部基因（独立过滤的 padj=NA 行 stat 仍在）——GSEA 排序要用全量
  sig_df <- res_df[res_df$sig, , drop = FALSE]
  sig_df <- sig_df[order(sig_df$padj), , drop = FALSE]
  write.table(sig_df, file = file.path(ct_dir, "significant_degs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(sub_meta, file = file.path(ct_dir, "metadata.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)   # pseudobulk_id 已是列
  saveRDS(list(counts = cnt, colData = sub_meta, design = row$design),
          file = file.path(ct_dir, "pseudobulk_counts.rds"))

  # 火山图（期刊风）：发散对（暖红上调/冷蓝下调）+ 中性灰不显著点；标签墨色
  vd <- res_df[!is.na(res_df$padj) & is.finite(res_df[[lfc_col]]), , drop = FALSE]
  vd$cat <- factor(ifelse(vd$padj < padj_cut & vd[[lfc_col]] >= lfc_cut, "Up",
                   ifelse(vd$padj < padj_cut & vd[[lfc_col]] <= -lfc_cut, "Down", "NS")),
                   levels = c("Down", "NS", "Up"))
  vd$mlog10p <- pmin(-log10(vd$padj), 350)   # -log10(padj) 截断（极显著基因拉爆坐标轴）
  sig_idx <- which(vd$cat != "NS")
  lab <- head(vd[sig_idx[order(vd$padj[sig_idx])], , drop = FALSE], 10)   # top10 显著基因标签
  # 分层绘制：NS 先画（更小更淡，退到背景），显著点后画（更实，压在最上）——防灰点淹没显著点
  p <- ggplot(vd, aes(x = .data[[lfc_col]], y = mlog10p, color = cat)) +
    geom_point(data = vd[vd$cat == "NS", , drop = FALSE], size = 0.5, alpha = 0.25) +
    geom_point(data = vd[vd$cat != "NS", , drop = FALSE], size = 0.9, alpha = 0.75) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed",
               color = BASELINE, linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed",
               color = BASELINE, linewidth = 0.5) +
    scale_color_manual(values = c(Up = PAL_UP, Down = PAL_DOWN, NS = PAL_NS), drop = FALSE) +
    labs(title = ct,
         subtitle = sprintf("Up %d · Down %d", sum(vd$cat == "Up"), sum(vd$cat == "Down")),
         x = paste0("log2 fold change", if (lfc_col != "log2FoldChange") " (shrunk)" else ""),
         y = expression(-log[10]~adjusted~italic(P)), color = NULL) +
    theme_pb()
  if (nrow(lab) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.8,
                                      color = INK, max.overlaps = 20,
                                      segment.color = MUTED, segment.size = 0.3,
                                      box.padding = 0.3, show.legend = FALSE)
  }
  save_fig(p, "volcano", type = "volcano")

  # MA 图（期刊风 ggplot 版）：x = log10(baseMean) 看计数深度，y = 收缩 log2FC；
  # 不显著点退灰，阈值线为虚线基线色
  ma_df <- res_df[is.finite(res_df$baseMean) & is.finite(res_df[[lfc_col]]), , drop = FALSE]
  ma_df$log10mean <- log10(ma_df$baseMean + 1)
  ma_df$cat <- factor(ifelse(ma_df$sig & ma_df[[lfc_col]] > 0, "Up",
                      ifelse(ma_df$sig & ma_df[[lfc_col]] < 0, "Down", "NS")),
                      levels = c("Down", "NS", "Up"))
  p <- ggplot(ma_df, aes(x = log10mean, y = .data[[lfc_col]], color = cat)) +
    geom_point(size = 0.8, alpha = 0.55) +
    geom_hline(yintercept = 0, color = BASELINE, linewidth = 0.5) +
    geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dashed",
               color = BASELINE, linewidth = 0.5) +
    scale_color_manual(values = c(Up = PAL_UP, Down = PAL_DOWN, NS = PAL_NS), drop = FALSE) +
    labs(title = ct, x = expression(log[10]~mean~normalized~counts),
         y = paste0("log2 fold change", if (lfc_col != "log2FoldChange") " (shrunk)" else ""),
         color = NULL) +
    theme_pb()
  save_fig(p, "MAplot", width = 7, height = 6)

  # vst 变换（PCA 与热图共用）：失败时 PCA 退回 log2 计数、热图跳过
  vsd <- tryCatch(DESeq2::vst(dds, blind = FALSE), error = function(e) {
    message("-- vst 失败，PCA 退回 log2 计数、热图跳过: ", conditionMessage(e))
    NULL
  })
  pc_src <- if (!is.null(vsd)) {
    tryCatch(DESeq2::plotPCA(vsd, intgroup = "condition", returnData = TRUE),
             error = function(e) NULL)
  } else NULL
  if (!is.null(pc_src)) {
    pc_src$sample_id <- as.character(sub_meta$sample_id)   # plotPCA 行序与 colData 一致
    save_fig(plot_pb_pca(pc_src, attr(pc_src, "percentVar"), ct), "PCA")
  } else {
    pcs <- pca_scores(cnt)   # cnt 已是 基因×样本，pca_scores 正期待这个方向
    pcs$df$condition <- sub_meta$condition
    pcs$df$sample_id <- as.character(sub_meta$sample_id)
    save_fig(plot_pb_pca(pcs$df, pcs$pve, ct), "PCA")
  }

  # 热图（期刊标准 pseudobulk 图）：top N 显著 DEG × 本类型样本（vst 值），
  # 列注释 = 条件颜色条；期望样本按条件聚成两簇
  top_deg <- head(sig_df$gene, 50)
  if (!is.null(vsd) && length(top_deg) >= 2 && requireNamespace("pheatmap", quietly = TRUE)) {
    hm <- SummarizedExperiment::assay(vsd)[top_deg, , drop = FALSE]
    hm <- t(scale(t(hm)))   # 行 z-score：基因间可比
    colnames(hm) <- paste0(as.character(sub_meta$sample_id), "_", sub_meta$condition)
    ann <- data.frame(condition = sub_meta$condition, row.names = colnames(hm))
    save_fig_draw(
      pheatmap::pheatmap(
        hm, annotation_col = ann,
        annotation_colors = list(condition = setNames(c(PAL_CTRL, PAL_STIM),
                                                      c(ref_grp, test_grp))),
        main = paste0(ct, " — top ", nrow(hm), " DEGs"),
        color = colorRampPalette(c(PAL_DOWN, "#f5f5f2", PAL_UP))(100),   # 发散对：蓝→白→红
        fontsize = 8, fontsize_row = 6, border_color = NA),
      "heatmap", width = 8, height = 7)
  } else if (length(top_deg) < 2) {
    message("-- ", ct, " 显著 DEG 不足 2 个，跳过热图")
  }

  # 可选 GSEA：全基因按 Wald stat 排序（stat = logFC/SE，低计数噪声基因天然被压后）
  if (run_gsea) {
    geneList <- res_df$stat
    names(geneList) <- res_df$gene
    geneList <- sort(geneList[!is.na(geneList)], decreasing = TRUE)
    gmt_files <- list(kegg = cfg$diff_gsea$gmt_kegg, hallmark = cfg$diff_gsea$gmt_hallmark)
    for (nm in names(gmt_files)) {
      gmt <- gmt_files[[nm]]
      if (is.null(gmt) || !file.exists(gmt)) {
        warning("基因集文件不存在，跳过 ", nm, ": ", gmt)
        next
      }
      set.seed(123)   # ⚠️ GSEA 靠随机置换算 p 值，与 05 同款固定种子
      y <- GSEA(geneList, TERM2GENE = read.gmt(gmt))
      y_df <- as.data.frame(y)
      if (nrow(y_df) == 0) {   # 零结果保护：gmt 与数据基因名不匹配时防晦涩报错（同 05）
        warning(nm, " 基因集 GSEA 无任何富集结果（检查 gmt 物种与格式）: ", gmt)
        next
      }
      write.table(y_df, file = file.path(ct_dir, paste0("gsea_", nm, ".tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      # ⚠️ 正负两方向都存在时才按 .sign 拆 facet（单向结果 facet 会中断脚本，同 05）
      if (length(unique(y_df$.sign)) >= 2) {
        save_fig(dotplot(y, showCategory = 12, split = ".sign") + facet_grid(~.sign),
                 paste0("gsea_", nm, "_dotplot"), type = "gsea")
      } else {
        save_fig(dotplot(y, showCategory = 12),
                 paste0("gsea_", nm, "_dotplot"), type = "gsea")
      }
    }
  }

  # ---- 同对比口径的单细胞对照（伪重复基准）----
  # FindMarkers 以细胞为统计单元：同一供者的细胞被当独立观测，p 值虚小。
  # 与 DESeq2 取相同基因全集 + 相同 BH 校正后比较显著数（STIM 上调），
  # 供 summary_degs.csv 与 summary_05_vs_06 图对照——细胞级口径应 ≥ 样本级口径
  fm <- tryCatch(
    FindMarkers(scobj, ident.1 = paste0(ct, "^^", test_grp),
                ident.2 = paste0(ct, "^^", ref_grp),
                logfc.threshold = 0, min.pct = 0, verbose = FALSE),
    error = function(e) NULL)
  if (!is.null(fm)) {
    common <- intersect(res_df$gene, rownames(fm))
    r_cmp <- res_df[match(common, res_df$gene), ]
    f_cmp <- fm[match(common, rownames(fm)), ]
    r_cmp$padj_bh <- p.adjust(r_cmp$pvalue, "BH")
    f_cmp$padj_bh <- p.adjust(f_cmp$p_val, "BH")
    row$n_pb_up_bh <- sum(!is.na(r_cmp$padj_bh) & r_cmp$padj_bh < padj_cut & r_cmp[[lfc_col]] > 0)
    row$n_sc_up_bh <- sum(!is.na(f_cmp$padj_bh) & f_cmp$padj_bh < padj_cut & f_cmp$avg_log2FC > 0)
  } else {
    warning("类型 ", ct, " 的单细胞 FindMarkers 对照失败（对照列留空）")
  }

  row$n_tested_genes <- nrow(res_df)
  row$n_sig_padj_lfc <- nrow(sig_df)
  row$n_up <- sum(vd$cat == "Up")
  row$n_down <- sum(vd$cat == "Down")
  row$n_sig_up <- as.integer(n_padj_only_up)
  summary_rows[[ct]] <- row
}
options(fig_outdir = file.path(step_dir, "figures"))   # 恢复到本步的 figures/

# ---- 7. 汇总 + 05/06 对比 ----
if (length(summary_rows) == 0) {
  stop("没有任何细胞类型通过 min_samples_per_group 守卫，请检查 config 或数据")
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file = file.path(step_dir, "summary_degs.csv"), row.names = FALSE)

# 头条数字：同对比口径（同基因全集 + 同 BH）下，细胞级（伪重复）口径的显著数
# 应 ≥ 样本级口径——差多少就是伪重复虚高了多少
tot_sc <- sum(summary_df$n_sc_up_bh, na.rm = TRUE)
tot_pb <- sum(summary_df$n_pb_up_bh, na.rm = TRUE)
message("-- 同对比口径显著 DEG（padj<", padj_cut, ", ", test_grp,
        " 上调；相同基因全集 + BH）: 单细胞 FindMarkers ", tot_sc,
        " vs pseudobulk DESeq2 ", tot_pb)

cmp <- summary_df[!is.na(summary_df$n_sc_up_bh) & !is.na(summary_df$n_pb_up_bh), ,
                  drop = FALSE]
if (nrow(cmp) > 0) {
  cmp_long <- rbind(
    data.frame(celltype = cmp$celltype, method = "single-cell FindMarkers (cell-level)",
               n = cmp$n_sc_up_bh),
    data.frame(celltype = cmp$celltype, method = "pseudobulk DESeq2 (sample-level)",
               n = cmp$n_pb_up_bh))
  p <- ggplot(cmp_long, aes(x = celltype, y = n, fill = method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("single-cell FindMarkers (cell-level)" = PAL_CTRL,
                                 "pseudobulk DESeq2 (sample-level)" = PAL_STIM)) +
    labs(y = paste0("显著 DEG 数（padj<", padj_cut, ", ", test_grp, " 上调）"), x = NULL,
         title = "Same-contrast DEG counts: cell-level vs sample-level inference") +
    theme_pb()
  save_fig(p, "summary_05_vs_06", type = "barplot")
}

log_step("模块 9+ 完成 → 下一步: Rscript R/07_gene_set_score.R（--with-advanced）或全流程结束")
