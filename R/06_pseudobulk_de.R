###############################################################################
### 06_pseudobulk_de.R — 模块 9+：pseudobulk 确认性差异分析（sample-aware DESeq2）
###
### 05 以细胞为统计单元（伪重复，p 值虚高）；本模块按 样本×细胞类型 聚合原始 counts，
### 逐类型独立建 DESeq2：同一样本出现在 ≥2 组 → 配对 design ~ sample_id + condition，
### 否则 ~ condition；batch 可用且未与 condition 混淆时插到最前。
### ⚠️ 只吃 RNA assay 的 raw counts（DESeq2 负二项模型不吃标准化/SCT/Harmony 产物）。
###
### 输入：output/<batch>/03_annotate/seurat_annotated.rds（sample_id 由 01 写入）
### 输出：output/<batch>/06_pseudobulk_de/ —— 逐类型 DESeq2 全基因结果/显著表/
###      volcano/MA/PCA；顶层 summary_degs.csv 含与单细胞口径的同对比对照
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
exclude_ct     <- get_cfg(pb$exclude_celltypes, NULL)
if (!identical(get_cfg(pb$method, "DESeq2"), "DESeq2")) {
  stop("pseudobulk$method 目前仅支持 DESeq2，收到: ", pb$method)
}
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
message("-- 聚合用细胞: ", nrow(cell_df), "/", n_before, "（剔除无供者信息与排除类型）")

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

# ---- 5. 全局样本层 PCA（QC）：组分离与供者离群一眼可见 ----
#' 取 top 高变 pseudobulk 做 PCA，返回得分与解释度
pca_scores <- function(mat, ntop = 2000) {
  mat <- log2(as.matrix(mat) + 1)
  v <- apply(mat, 1, var)
  mat <- mat[order(v, decreasing = TRUE), , drop = FALSE]
  mat <- mat[seq_len(min(ntop, nrow(mat))), , drop = FALSE]
  pc <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  list(df = as.data.frame(pc$x), pve = summary(pc)$importance[2, ])
}
if (nrow(pb_meta) >= 6) {
  pcs <- pca_scores(counts_pb[, pb_meta$pseudobulk_id])   # counts_pb 已是 基因×样本
  pcs$df$condition <- pb_meta$condition   # 颜色分组（pca_scores 只返回坐标）
  pcs$df$sample_id <- as.character(pb_meta$sample_id)   # 供者标签（复合键拆出的干净 id）
  p <- ggplot(pcs$df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3, alpha = 0.85) +
    geom_text(aes(label = sample_id), vjust = -0.9, size = 2.6, show.legend = FALSE) +
    labs(x = sprintf("PC1 (%.1f%%)", 100 * pcs$pve[1]),
         y = sprintf("PC2 (%.1f%%)", 100 * pcs$pve[2]),
         title = "Pseudobulk samples (all celltypes, log2 counts PCA)") +
    theme_bw()
  save_fig(p, "all_types_pca")
}

# ---- 6. GSEA 预备（可选；默认关）----
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

# ---- 7. 逐细胞类型 DESeq2 ----
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
for (ct in cts_all) {
  ct_dir <- file.path(step_dir, ct_dirs[[ct]])
  sub_meta <- droplevels(pb_meta[pb_meta$celltype == ct, , drop = FALSE])
  row <- data.frame(celltype = ct, skip_reason = NA_character_,
                    design = NA_character_,
                    n_pseudobulk_ref = sum(sub_meta$condition == ref_grp),
                    n_pseudobulk_test = sum(sub_meta$condition == test_grp),
                    n_tested_genes = NA_integer_,
                    n_sig_padj_lfc = NA_integer_, n_up = NA_integer_, n_down = NA_integer_,
                    n_sig_up = NA_integer_,
                    n_sc_up_bh = NA_integer_, n_pb_up_bh = NA_integer_)
  # 注释里有但聚合后无样本的类型（min_cells 全格不足 / config 剔除）——
  # 明确记入 summary 而非无声消失
  if (!ct %in% pb_meta$celltype) {
    row$skip_reason <- if (!is.null(exclude_ct) && ct %in% exclude_ct) "config 剔除"
                       else "min_cells 过滤后无 pseudobulk 样本"
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
  options(fig_outdir = ct_dir)   # 本类型的图直接落子目录（save_fig 读该 option）
  cnt <- as.matrix(counts_pb[, sub_meta$pseudobulk_id, drop = FALSE])   # 基因×样本
  storage.mode(cnt) <- "integer"   # rowsum 产出 double；DESeq2 要整数（避免静默 round 警告）
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt, colData = sub_meta, design = fml)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
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

  # 火山图（手写 ggplot：不引入 EnhancedVolcano，配色与仓库其余图统一）
  vd <- res_df[!is.na(res_df$padj) & is.finite(res_df[[lfc_col]]), , drop = FALSE]
  vd$cat <- ifelse(vd$padj < padj_cut & vd[[lfc_col]] >= lfc_cut, "Up",
            ifelse(vd$padj < padj_cut & vd[[lfc_col]] <= -lfc_cut, "Down", "NS"))
  vd$cat <- factor(vd$cat, levels = c("Down", "NS", "Up"))
  sig_idx <- which(vd$cat != "NS")
  lab <- head(vd[sig_idx[order(vd$padj[sig_idx])], , drop = FALSE], 10)   # top10 显著基因标签
  p <- ggplot(vd, aes(x = .data[[lfc_col]], y = -log10(padj), color = cat)) +
    geom_point(size = 0.7, alpha = 0.6) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c(Up = "#C0392B", Down = "#2E86C1", NS = "grey70"), drop = FALSE) +
    labs(title = ct, x = paste0("log2 fold change", if (lfc_col != "log2FoldChange") " (shrunk)" else ""),
         y = "-log10(adjusted p)", color = NULL) +
    theme_bw()
  if (nrow(lab) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.6,
                                      max.overlaps = 20, show.legend = FALSE)
  }
  save_fig(p, "volcano", type = "volcano")

  # MA 图：plotMA 直接绘制不返回对象 → save_fig_draw（设备内求值，见 save_fig.R 说明）
  save_fig_draw(DESeq2::plotMA(dds, alpha = padj_cut), "MAplot", width = 7, height = 6)

  # 样本层 PCA：优先 vst + plotPCA(returnData)（官方变换），失败退回 log2 计数 PCA
  pc_src <- tryCatch({
    vsd <- DESeq2::vst(dds, blind = FALSE)
    DESeq2::plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  }, error = function(e) {
    message("-- vst/plotPCA 失败，退回 log2 计数 PCA: ", conditionMessage(e))
    NULL
  })
  if (!is.null(pc_src)) {
    pc_src$sample_id <- as.character(sub_meta$sample_id)   # plotPCA 行序与 colData 一致
    pve <- attr(pc_src, "percentVar")
    p <- ggplot(pc_src, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(size = 3, alpha = 0.85) +
      geom_text(aes(label = sample_id), vjust = -0.9, size = 2.6, show.legend = FALSE) +
      labs(x = sprintf("PC1: %.1f%% variance", 100 * pve[1]),
           y = sprintf("PC2: %.1f%% variance", 100 * pve[2]), title = ct) +
      theme_bw()
    save_fig(p, "PCA")
  } else {
    pcs <- pca_scores(cnt)   # cnt 已是 基因×样本，pca_scores 正期待这个方向
    pcs$df$condition <- sub_meta$condition
    pcs$df$sample_id <- as.character(sub_meta$sample_id)
    p <- ggplot(pcs$df, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(size = 3, alpha = 0.85) +
      geom_text(aes(label = sample_id), vjust = -0.9, size = 2.6, show.legend = FALSE) +
      labs(x = sprintf("PC1 (%.1f%%)", 100 * pcs$pve[1]),
           y = sprintf("PC2 (%.1f%%)", 100 * pcs$pve[2]), title = ct) +
      theme_bw()
    save_fig(p, "PCA")
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

# ---- 8. 汇总 + 05/06 对比 ----
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
    scale_fill_manual(values = c("#95A5A6", "#C0392B")) +
    labs(y = paste0("显著 DEG 数（padj<", padj_cut, ", ", test_grp, " 上调）"), x = NULL,
         title = "Same-contrast DEG counts: cell-level vs sample-level inference") +
    theme_bw()
  save_fig(p, "summary_05_vs_06", type = "barplot")
}

log_step("模块 9+ 完成 → 下一步: Rscript R/07_gene_set_score.R（--with-advanced）或全流程结束")
