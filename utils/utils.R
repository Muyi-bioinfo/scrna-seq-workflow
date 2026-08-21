###############################################################################
### utils.R — 全流程共享工具库
###
### 每个分步脚本（R/01_*.R ~ R/08_*.R）的前三行都是：
###   source("utils/utils.R")
###   cfg <- load_config()                  # 读取 config/config.yaml
###   step_dir <- setup_step("NN_xxx")      # 设置本步的输出目录与图目录
###
### 设计原则：参数全部集中在 config；每步输出到 output/<batch>/<步骤名>/；
### 保存前自动瘦身（清空 scale.data）；只用官方 API（[[]] 等），不用 @ 访问 slot
###############################################################################

# 基础依赖
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# 图片保存工具（cairo_pdf 封装，支持中文；图目录受 options(fig_outdir) 控制）
source("utils/save_fig.R")

# 项目根目录（脚本都从项目根目录运行：Rscript R/01_xxx.R）
PROJECT_ROOT <- normalizePath(".")

# 配置文件：默认 config/config.yaml（单样本），多样本用 config/config.multi.yaml
# 可由环境变量 CONFIG_FILE 覆盖（run_pipeline.sh --mode multi 会自动切换）
CONFIG_FILE  <- Sys.getenv("CONFIG_FILE", "config/config.yaml")

###############################################################################
### 配置与路径
###############################################################################

#' 读取集中配置（嵌套结构，如 cfg$qc$nfeature_min）
load_config <- function(file = CONFIG_FILE) {
  if (!file.exists(file)) stop("找不到配置文件: ", file)
  cfg <- yaml::read_yaml(file)
  message("-- 已加载配置: ", file)
  cfg
}

#' 打印步骤日志
log_step <- function(msg) {
  message("============================================================================")
  message("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message("============================================================================")
}

#' 本批次的结果根目录：output/<batch>/（batch 字段在 config；同批次重跑 = 断点续跑）
batch_dir <- function() {
  if (is.null(cfg$batch) || !nzchar(cfg$batch)) {
    stop("config 缺少 batch 字段（结果归档目录），请先填写")
  }
  file.path(cfg$dirs$base, cfg$batch)
}

#' 设置本步的工作目录（结果子目录 + 图子目录，图路径写入 options(fig_outdir)）
setup_step <- function(step_name) {
  if (!exists("cfg")) stop("请先执行 cfg <- load_config()")
  step_dir <- file.path(batch_dir(), step_name)
  dir.create(file.path(step_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  options(fig_outdir = file.path(step_dir, "figures"))   # 图输出到本步的 figures/
  step_dir
}

#' 上一步输出文件的路径（不读入；用于"缺失则跳过而非报错"的场景）
prev_path <- function(prev_step, file) {
  file.path(batch_dir(), prev_step, file)
}

#' 读取上一步的输出（分步脚本的标准入口，同批次内传递）
read_prev <- function(prev_step, file) {
  fp <- prev_path(prev_step, file)
  if (!file.exists(fp)) stop("找不到上一步输出: ", fp, " —— 请先运行上一步脚本")
  message("-- 读入: ", fp)
  readRDS(fp)
}

#' 保存前瘦身 + saveRDS
#' 清空 scale.data（基因数 × 细胞数的稠密矩阵，体积最大），不影响已有降维和分群
slim_save <- function(scobj, file, step_dir) {
  scobj[["RNA"]]$scale.data <- matrix()
  fp <- file.path(step_dir, file)
  saveRDS(scobj, file = fp)
  message("-- Saved: ", fp, " (", round(file.size(fp) / 1e6, 1), " MB)")
}

###############################################################################
### 可复用分析函数（各分步脚本共用）
###############################################################################

#' 质控：线粒体比例 → 小提琴图 → 细胞过滤（阈值全部来自 cfg$qc）
qc_filter <- function(scobj, cfg) {
  log_step(paste("质控 | 过滤前细胞数:", ncol(scobj)))

  # 计算每个细胞的线粒体基因比例，写入 metadata 新增列 percent.mt
  scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = cfg$qc$mt_pattern)

  save_fig(VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3), "qc_violin_before", type = "vlnplot")

  scobj <- subset(scobj,
                  subset = nFeature_RNA > cfg$qc$nfeature_min &
                    nFeature_RNA < cfg$qc$nfeature_max &
                    percent.mt < cfg$qc$percent_mt_max)
  message("-- 过滤后细胞数: ", ncol(scobj))
  scobj
}

#' 双细胞检测（scDblFinder）
#'
#' 检测结果写入 metadata 的 doublet 列（"singlet"/"doublet"），由调用方过滤。
#' 多样本场景应逐样本检测（双细胞率随建库批次变化），在 merge 前调用。
detect_doublets <- function(scobj) {
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    stop("未安装 scDblFinder（conda 安装: mamba install -c bioconda bioconductor-scdblfinder），",
         "或把 config 的 qc$doublet_enable 设为 false")
  }
  log_step(paste("双细胞检测（scDblFinder）| 细胞数:", ncol(scobj)))
  sce <- as.SingleCellExperiment(scobj)          # Seurat → SCE（用 counts 层）
  sce <- scDblFinder::scDblFinder(sce)           # 模拟双细胞 + 分类
  scobj[["doublet"]] <- sce$scDblFinder.class    # 结果写入 metadata
  message("-- 双细胞: ", sum(scobj$doublet == "doublet"),
          " / 单细胞: ", sum(scobj$doublet == "singlet"))
  scobj
}

#' 预处理：标准化 → 高变基因 → 缩放 → PCA（参数全部来自 cfg$preprocess）
preprocess_scobj <- function(scobj, cfg) {
  log_step(paste("预处理 | nfeatures:", cfg$preprocess$nfeatures,
                 "| npcs:", cfg$preprocess$npcs))

  # 1. 标准化：每个细胞的 UMI 总数缩放到 10000，再取自然对数
  scobj <- NormalizeData(scobj, normalization.method = "LogNormalize",
                         scale.factor = 10000)

  # 2. 高变基因：vst 方法找细胞间方差最大的基因
  scobj <- FindVariableFeatures(scobj, selection.method = "vst",
                                nfeatures = cfg$preprocess$nfeatures)

  top10 <- head(VariableFeatures(scobj), 10)
  save_fig(VariableFeaturePlot(scobj) +
             LabelPoints(plot = VariableFeaturePlot(scobj), points = top10, repel = TRUE),
           "variable_features", type = "default")

  # 3. 缩放：基因表达中心化，PCA 的必备前置
  #    ⚠️ 默认只缩放高变基因；scale_all = true 时缩放全部基因（内存消耗巨大）
  features_scale <- if (isTRUE(cfg$preprocess$scale_all)) rownames(scobj) else VariableFeatures(scobj)
  scobj <- ScaleData(scobj, features = features_scale)

  # 4. PCA：线性降维，结果存入 reduction "pca"
  scobj <- RunPCA(scobj, features = VariableFeatures(scobj),
                  npcs = cfg$preprocess$npcs, reduction.name = "pca")

  # ElbowPlot 帮助选择下游使用的维度数（拐点之后信息量骤降）
  save_fig(ElbowPlot(scobj), "elbow", type = "default")
  scobj
}

#' 聚类：UMAP 非线性降维 + SNN 图 + 分群
#'
#' @param scobj      Seurat 对象（需已有降维）
#' @param reduction  使用的降维（单样本 pca；多样本 harmony/integrated.rpca）
#' @param dims       使用的维度区间
#' @param resolution 分群分辨率
#' @param fig_name   保存的图名
cluster_cells <- function(scobj, reduction, dims, resolution,
                          fig_name = "dimplot_clusters") {
  log_step(paste("聚类 | reduction:", reduction, "| dims:",
                 paste(range(dims), collapse = "-"), "| resolution:", resolution))

  scobj <- RunUMAP(scobj, reduction = reduction, dims = dims, reduction.name = "umap")
  scobj <- FindNeighbors(scobj, reduction = reduction, dims = dims)   # KNN → SNN 图
  scobj <- FindClusters(scobj, resolution = resolution)               # resolution 越大分群越细

  save_fig(DimPlot(scobj, reduction = "umap", label = TRUE), fig_name, type = "dimplot")
  message("-- 分群数: ", length(unique(Idents(scobj))))
  scobj
}

#' 从配置中取维度区间：dims_start:dims_end
dims_from <- function(param) {
  seq(param$dims_start, param$dims_end)
}

###############################################################################
### 轨迹预设（07_trajectory_analysis 专用）
###############################################################################

#' 加载轨迹预设
#'
#' 预设 = 轨迹分析的谱系先验（哪些细胞构成连续过程 + 谁是起点），
#' 文件位于 config/trajectory_presets/<name>.yaml（写作清单见该目录 README）。
#' 优先级：主 config 显式参数 > 预设 > 默认值。
load_trajectory_preset <- function(traj_cfg) {
  preset_name <- traj_cfg$preset
  if (is.null(preset_name) || !nzchar(preset_name)) {
    message("-- 未指定 preset，使用全部细胞（无谱系先验，仅保证流程跑通）")
    return(list(preset_name = "none", description = "全部细胞（无谱系先验）",
                root_celltype = traj_cfg$root_celltype,
                subset_celltypes = traj_cfg$subset_celltypes,
                num_dim = traj_cfg$num_dim))
  }

  preset_file <- file.path("config/trajectory_presets", paste0(preset_name, ".yaml"))
  if (!file.exists(preset_file)) {
    stop("找不到轨迹预设: ", preset_file,
         "，可用预设见 config/trajectory_presets/README.md")
  }
  preset <- yaml::read_yaml(preset_file)
  log_step(paste("轨迹预设:", preset_name, "|", preset$description))

  # 主 config 显式参数优先于预设
  root_ct <- traj_cfg$root_celltype
  if ((is.null(root_ct) || !nzchar(root_ct)) && !is.null(preset$root_celltype)) {
    root_ct <- preset$root_celltype
  }
  subs <- traj_cfg$subset_celltypes
  if (is.null(subs) && !is.null(preset$subset_celltypes)) subs <- preset$subset_celltypes
  num_dim <- traj_cfg$num_dim
  if (is.null(num_dim) && !is.null(preset$num_dim)) num_dim <- preset$num_dim

  list(preset_name = preset_name, description = preset$description,
       root_celltype = root_ct, subset_celltypes = subs, num_dim = num_dim)
}

#' 按预设模式对 Seurat 对象做谱系子集
#'
#' 支持精确细胞类型名与 * / ? 通配模式（如 "CD4*" 匹配所有 CD4 开头的类型）；
#' 全部失配时 stop 并打印数据中可用的 celltype 列表（"自报家底"）。
subset_by_celltypes <- function(scobj, patterns) {
  if (is.null(patterns) || length(patterns) == 0) return(scobj)

  all_ct <- sort(unique(scobj$celltype))
  selected <- character(0)
  for (pat in patterns) {
    hits <- if (grepl("[*?]", pat)) grep(glob2rx(pat), all_ct, value = TRUE)
            else intersect(pat, all_ct)
    if (length(hits) == 0) message("-- 警告: 模式 \"", pat, "\" 未匹配到任何细胞类型")
    selected <- union(selected, hits)
  }

  if (length(selected) == 0) {
    stop("子集未匹配到任何细胞类型！数据中可用的 celltype:\n  ",
         paste(all_ct, collapse = ", "))
  }
  # 按细胞名筛选（与 08 脚本一致，避免 NSE 在不同环境下的兼容性问题）
  scobj <- subset(scobj, cells = colnames(scobj)[scobj$celltype %in% selected])
  message("-- 谱系子集: 保留 ", length(selected), " 个类型, ", ncol(scobj), " 个细胞:\n   ",
          paste(sort(selected), collapse = ", "))
  scobj
}

###############################################################################
### 注释预设（03/04 专用）
###############################################################################

#' 加载注释预设（细胞类型 → marker 基因）
#'
#' 预设文件位于 config/annotation_presets/<species>/<name>.yaml（物种分离文件夹）。
#' 优先级：config 显式 marker_genes（扁平列表）> 预设。
#' 校验：① species 字段与所在文件夹一致 ② marker 在数据中的存在率
#'        （缺失超过 50% 则 warning——人鼠混用/数据集不适用会立刻暴露）
#'
#' @param annotate_cfg cfg$annotate（主 config 的 annotate 段）
#' @param scobj        可选，传入 Seurat 对象时做基因存在性校验
#' @return 列表：preset_name, description, celltypes（原始结构）,
#'         flat_tier1（每类 tier1 拍平）, flat_all（tier1+tier2 拍平）
load_annotation_preset <- function(annotate_cfg, scobj = NULL) {
  # config 显式 marker_genes 优先（扁平列表，与旧用法兼容）
  inline <- annotate_cfg$marker_genes
  if (!is.null(inline) && length(inline) > 0) {
    message("-- 注释 marker 来源: config 显式 marker_genes（覆盖预设）")
    return(list(preset_name = "inline", description = "config 显式 marker_genes",
                celltypes = NULL, flat_tier1 = inline, flat_all = inline))
  }

  preset_name <- annotate_cfg$preset
  if (is.null(preset_name) || !nzchar(preset_name)) {
    stop("config 的 annotate 段需指定 preset（如 human/pbmc）或显式 marker_genes")
  }

  preset_file <- file.path("config/annotation_presets", paste0(preset_name, ".yaml"))
  if (!file.exists(preset_file)) {
    stop("找不到注释预设: ", preset_file, "，可用预设见 config/annotation_presets/README.md")
  }
  preset <- yaml::read_yaml(preset_file)

  # ① 物种校验：species 字段必须与所在文件夹一致（human/pbmc → human）
  folder_species <- basename(dirname(preset_name))
  if (!is.null(preset$species) && nzchar(preset$species) &&
      preset$species != folder_species) {
    stop("预设物种不一致: ", preset_name, " 位于 ", folder_species,
         "/ 但 species 字段为 ", preset$species)
  }

  log_step(paste("注释预设:", preset_name, "|", preset$description))

  # ② 结构统一：每类 tier1/tier2 拍平（去重保序）
  celltypes <- preset$celltypes
  tier_of <- function(ct, tier) {
    g <- ct$markers[[tier]]
    if (is.null(g)) character(0) else g
  }
  flat_tier1 <- unique(unlist(lapply(celltypes, tier_of, tier = "tier1"), use.names = FALSE))
  flat_all   <- unique(unlist(lapply(celltypes, function(ct) c(tier_of(ct, "tier1"),
                                                               tier_of(ct, "tier2"))),
                                    use.names = FALSE))

  # ③ marker 存在性校验（对照数据基因名）
  if (!is.null(scobj)) {
    genes <- rownames(scobj)
    missing <- setdiff(flat_all, genes)
    if (length(missing) > 0) {
      message("-- marker 在数据中缺失 ", length(missing), "/", length(flat_all),
              " 个: ", paste(head(missing, 20), collapse = ", "))
      if (length(missing) / length(flat_all) > 0.5) {
        warning("缺失 marker 超过 50%！请检查物种是否匹配（人鼠混用？）或预设是否适用")
      }
    }
  }

  list(preset_name = preset_name, description = preset$description,
       celltypes = celltypes, flat_tier1 = flat_tier1, flat_all = flat_all)
}

#' 校验 cluster_map 的类型名与注释预设的一致性（不拦截，只提示）
#'
#' 命中预设标准名 → 通过；命中别名 → 提示标准名（数据集惯用名可继续使用）；
#' 完全未知 → warning（可能是拼写错误或预设不适用）。
#'
#' @param cluster_map 命名向量（分群编号 → 细胞类型名），来自 cfg$annotate$cluster_map
#' @param preset      load_annotation_preset() 的返回值
check_annotation_names <- function(cluster_map, preset) {
  celltypes <- preset$celltypes
  if (is.null(celltypes) || length(celltypes) == 0) return(invisible(NULL))

  # 标准名 → 自身；别名 → 标准名
  alias_map <- character(0)
  for (nm in names(celltypes)) {
    aliases <- celltypes[[nm]]$aliases
    if (!is.null(aliases) && length(aliases) > 0) {
      for (al in aliases) alias_map[al] <- nm
    }
  }

  for (v in unique(as.character(cluster_map))) {
    if (v %in% names(celltypes)) next
    if (v %in% names(alias_map)) {
      message("-- cluster_map 的 \"", v, "\" 是预设别名，标准名为 \"", alias_map[[v]], "\"")
    } else {
      warning("cluster_map 的 \"", v, "\" 不在注释预设的细胞类型中，请检查拼写或预设选择")
    }
  }
  invisible(NULL)
}
