###############################################################################
### save_fig.R — 图片保存通用工具（cairo_pdf 封装，支持中文）
###   save_fig(p, "01_dimplot")                  # 默认 7×7
###   save_fig(p, "02_dotplot", type = "dotplot") # 按类型自适应尺寸
###   save_fig(p, "03_custom", width = 10, height = 8)
###   save_fig(p, "04_x", enabled = FALSE)        # 静默模式，不保存
###############################################################################

# 图输出目录：默认 output/figures/；分步脚本中 setup_step() 用 options(fig_outdir) 指向各步目录
fig_dir <- function() {
  getOption("fig_outdir", "output/figures")
}

ensure_fig_dir <- function() {
  dir.create(fig_dir(), recursive = TRUE, showWarnings = FALSE)
}

# 预设尺寸表（可扩展）
# 两类图分开考虑：
#   1. 论文/汇报图（预设管）：期刊双栏 ≈ 7in、单栏 ≈ 3.5in，UMAP 面板惯例 3~4in
#      dimplot 7×7 可直接用作双栏图；featureplot 预设面向单基因面板（如 06 每通路一张）
#   2. 人眼质检工作图（显式尺寸管，不进论文）：以「fit-width 一屏放下」为准——
#      如 03 每类型一张的 FeaturePlot（动态 ncol：≤4 基因 2 列 / ≥5 基因 3 列，
#      行数封顶 2 行，最坏 14×10.4in 不滚动）
#   多面板 grid 一律用显式 width/height 按面板数计算，面板尽量保持近似方形
fig_sizes <- list(
  dimplot     = c(7, 7),
  dotplot     = c(10, 6),
  heatmap     = c(8, 7),
  featureplot = c(6, 6),
  vlnplot     = c(8, 6),
  volcano     = c(7, 6),
  gsea        = c(8, 6),
  barplot     = c(8, 5),
  default     = c(7, 7)
)

#' 保存图片为 PDF（cairo_pdf 设备，支持中文）
#'
#' @param plot    ggplot / Seurat / base R 绘图对象
#' @param name    文件名（不含后缀），如 "01_dimplot"
#' @param type    预设尺寸: dimplot/dotplot/heatmap/featureplot/vlnplot/volcano/gsea/barplot
#' @param width   覆盖宽度（英寸），NULL 则由 type 决定
#' @param height  覆盖高度（英寸），NULL 则由 type 决定
#' @param enabled 是否输出图片，FALSE 时静默跳过
save_fig <- function(plot, name, type = NULL, width = NULL, height = NULL, enabled = TRUE) {
  if (!enabled) return(invisible(NULL))

  ensure_fig_dir()

  # 尺寸优先级：显式 width/height > type 预设 > default
  if (!is.null(type)) {
    size <- fig_sizes[[type]]
    if (is.null(size)) {
      warning("未知的 type '", type, "', 使用 default 尺寸 7×7")
      size <- fig_sizes$default
    }
  } else {
    size <- fig_sizes$default
  }

  if (is.null(width))  width  <- size[1]
  if (is.null(height)) height <- size[2]

  filepath <- file.path(fig_dir(), paste0(name, ".pdf"))

  # ggplot 对象用 ggsave；其他对象用 cairo_pdf + print
  if (inherits(plot, "gg") || inherits(plot, "ggplot")) {
    ggplot2::ggsave(filepath, plot, device = cairo_pdf, width = width, height = height)
  } else {
    cairo_pdf(filepath, width = width, height = height)
    print(plot)
    dev.off()
  }

  message("-- Saved: ", filepath, " (", width, "×", height, "in)")
}

#' 保存"调用时直接绘制"的 base 图形（如 CellChat 的 netVisual_*）
#'
#' ⚠️ 这类函数不返回可重放的绘图对象——Rscript 非交互模式下 recordPlot()
#' 捕获到空 display list，save_fig() 的 print 重放路径失效（图实际画到了
#' 自动打开的 Rplots.pdf 上）。只能在设备打开的状态下执行绘制表达式直接捕获。
#'
#' @param expr   绘制表达式（在 cairo_pdf 设备内求值，惰性求值保证先开设备）
#' @param name   文件名（不含后缀）
#' @param width  宽度（英寸），默认 7
#' @param height 高度（英寸），默认 7
save_fig_draw <- function(expr, name, width = 7, height = 7) {
  ensure_fig_dir()
  filepath <- file.path(fig_dir(), paste0(name, ".pdf"))
  cairo_pdf(filepath, width = width, height = height)
  force(expr)   # 设备打开后才执行绘制
  dev.off()
  message("-- Saved: ", filepath, " (", width, "×", height, "in)")
}
