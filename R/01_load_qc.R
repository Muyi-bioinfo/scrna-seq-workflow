###############################################################################
### 01_load_qc.R — 模块 1+2：数据读入 + 质控
###
### 按 config 的 mode 分支：single 读 10x 矩阵 / multi 读样本表 merge 后质控
### 输入：config 对应段落；输出：output/<batch>/01_load_qc/seurat_qc.rds + 质控图
###############################################################################

source("utils/utils.R")
cfg <- load_config()
step_dir <- setup_step("01_load_qc")
log_step(paste("模块 1+2 | 数据读入 + 质控 | 模式:", cfg$mode))

# 读入阈值与质控同源（cfg$qc）：min.cells 基因至少 N 细胞表达，min.features 细胞至少 N 基因
min_cells <- if (!is.null(cfg$qc$min_cells)) cfg$qc$min_cells else 3

# ---- 第一步：读入 ----
if (cfg$mode == "single") {
  scdata <- Read10X(data.dir = cfg$single$data_dir)   # 自动识别三件套/h5

  scobj <- CreateSeuratObject(counts = scdata,
                              project = cfg$single$project,
                              min.cells = min_cells,
                              min.features = cfg$qc$nfeature_min)
  # 单样本有多组时在 single$group 填组名
  if (!is.null(cfg$single$group) && nzchar(cfg$single$group)) {
    scobj[["group"]] <- cfg$single$group
  }

  # 双细胞检测（qc$doublet_enable 控制）
  if (isTRUE(cfg$qc$doublet_enable)) {
    scobj <- detect_doublets(scobj)
  }
} else {
  # 多样本：读样本表（sample, group, fastqs, matrix）；matrix 列优先，
  # 为空时自动找本批次 cellranger 输出 output/<batch>/00_cellranger/<sample>/outs/
  sheet <- read.csv(cfg$multi$sample_sheet, stringsAsFactors = FALSE)
  message("-- 样本表: ", cfg$multi$sample_sheet)
  print(sheet)

  sample_list <- lapply(seq_len(nrow(sheet)), function(i) {
    sample_name <- sheet$sample[i]
    matrix_path <- sheet$matrix[i]

    # ⚠️ read.csv 对空单元格返回 NA，且 nzchar(NA) 也是 NA——必须先 is.na 再 nzchar
    if (is.na(matrix_path) || !nzchar(matrix_path)) {
      matrix_path <- file.path(batch_dir(), "00_cellranger", sample_name,
                               "outs", "filtered_feature_bc_matrix.h5")
    }
    if (!file.exists(matrix_path)) {
      stop("找不到样本 ", sample_name, " 的矩阵: ", matrix_path,
           " —— 请检查样本表 matrix 列，或先运行 00_run_cellranger.sh（bash/）")
    }
    message("  -- ", sample_name, " (", sheet$group[i], "): ", matrix_path)

    scdata <- Read10X(data.dir = matrix_path)
    x <- CreateSeuratObject(counts = scdata,
                            project = sample_name,
                            min.cells = min_cells,
                            min.features = cfg$qc$nfeature_min)
    x[["group"]] <- sheet$group[i]   # 分组标签来自样本表（去批次的依据）
    message("     -> ", ncol(x), " cells | ", nrow(x), " genes")
    x
  })

  # 双细胞检测：逐样本检测（双细胞率随建库批次变化），再 merge
  if (isTRUE(cfg$qc$doublet_enable)) {
    sample_list <- lapply(sample_list, detect_doublets)
  }
  # 样本表只有 1 行时无需 merge（merge 传空 list 会报错）
  if (length(sample_list) == 1) {
    scobj <- sample_list[[1]]
  } else {
    scobj <- merge(x = sample_list[[1]], y = sample_list[-1])
  }
  rm(sample_list); gc()
  # ⚠️ merge 后是 split layers：harmony 路径在 02 里 JoinLayers，
  #    rpca/mnn 路径需保留 split layers 供 IntegrateLayers 使用
}

message("-- 读入完成: ", ncol(scobj), " cells | ", nrow(scobj), " genes")

# ---- 过滤双细胞（已检测时 doublet 列存在）----
if ("doublet" %in% colnames(scobj[[]])) {
  scobj <- subset(scobj, doublet == "singlet")
  message("-- 过滤双细胞后: ", ncol(scobj), " cells")
}

# ---- 第二步：质控（阈值全部来自 cfg$qc，函数在 utils/utils.R）----
scobj <- qc_filter(scobj, cfg)

slim_save(scobj, "seurat_qc.rds", step_dir)
log_step("模块 1+2 完成 → 下一步: Rscript R/02_preprocess_cluster.R")
