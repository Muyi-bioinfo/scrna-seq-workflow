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

  # 逐细胞元数据 join（可选）：合并矩阵场景下样本表每行是"多供者池"（如 GSE96583，
  # 供者在 GEO tsne.df 的 ind 列），pseudobulk 需要真实生物学重复。挂法与 demuxlet
  # 教程同款：条码做行名 → AddMetaData（行名对齐条码、与顺序无关，未覆盖细胞自动 NA）。
  # 常规场景（cellranger 每样本一个供者）不配置即可，sample_id = 样本表行名。
  # ⚠️ GEO tsne.df 首列（条码）无表头：read.table 会把它当 rownames——用
  # 「表头字段数 vs 首行数据字段数」探测；连接不手动 close（读取函数会自动开关）
  join_cell_metadata <- function(x, meta_cfg, group_label) {
    open_con <- function() if (grepl("\\.gz$", meta_cfg$cell_metadata)) {
      gzfile(meta_cfg$cell_metadata)
    } else file(meta_cfg$cell_metadata)
    head_lines <- readLines(open_con(), n = 2)
    n_head <- length(strsplit(head_lines[1], "\t", fixed = TRUE)[[1]])
    n_data <- length(strsplit(head_lines[2], "\t", fixed = TRUE)[[1]])
    m <- read.delim(open_con(), stringsAsFactors = FALSE)
    bc_meta <- if (n_data == n_head + 1) rownames(m) else m[[1]]
    if (!meta_cfg$donor_col %in% colnames(m)) {
      stop("cell_metadata 找不到供者列 \"", meta_cfg$donor_col,
           "\"，现有列: ", paste(colnames(m), collapse = ", "))
    }
    grp_col <- meta_cfg$cell_metadata_group_col
    if (is.null(grp_col) || !grp_col %in% colnames(m)) {
      stop("multi$cell_metadata_group_col 需指定 cell_metadata 中的组别列")
    }
    # 本样本的行（组别名大小写归一：STIM ↔ stim），条码做行名交给 AddMetaData——
    # 未覆盖的细胞自动 NA（本函数只负责建列，不能与样本表默认值叠加，见调用处 ⚠️）
    sel <- toupper(as.character(m[[grp_col]])) == toupper(group_label)
    df <- data.frame(sample_id = as.character(m[[meta_cfg$donor_col]])[sel],
                     row.names = bc_meta[sel])
    n_na <- sum(!colnames(x) %in% rownames(df))
    if (n_na > 0) {
      message("  -- ", group_label, " 有 ", n_na, "/", ncol(x),
              " 个细胞无供者信息（06 pseudobulk 聚合时剔除；此处为读入时计数，QC 还会再滤一部分）")
    }
    AddMetaData(x, metadata = df)
  }

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
    # sample_id（生物学重复，pseudobulk 的统计单元）两种来源，二选一：
    # ① 默认 = 样本表每行（cellranger 等每样本一个供者的常规场景）
    # ② 配置 multi$cell_metadata → 按 (group, barcode) join 供者（合并矩阵场景）
    # ⚠️ 两条路不能叠加：AddMetaData 对已存在的列是"合并"（未覆盖细胞保留旧值），
    # 会漏出 "pbmc_stim" 这种池名混进 06 当伪重复——join 场景必须让它建新列
    if (!is.null(cfg$multi$cell_metadata) && nzchar(cfg$multi$cell_metadata)) {
      x <- join_cell_metadata(x, cfg$multi, sheet$group[i])
    } else {
      x[["sample_id"]] <- sheet$sample[i]
    }
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
