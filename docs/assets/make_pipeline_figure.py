###############################################################################
### make_pipeline_figure.py — 生成 README 嵌入的管线总览图
###
### 用法：python docs/assets/make_pipeline_figure.py（需 matplotlib，scrna 环境自带）
### 产出：docs/assets/pipeline_overview.png（README 嵌图）+ .svg（矢量版）
###
### 设计规范（dataviz 技能）：分类色只承担「阶段类别」身份——
###   蓝 #2a78d6 = 核心 R 流程 / 橙 #eb6834 = 上游 cellranger / 青绿 #1baf7a = 高级模块
###   数据产物（FASTQ/矩阵）用中性灰（非系列色）；文字一律用墨水色，不用系列色。
###   调色板已过 validate_palette.js（CVD ΔE 9.2 / 普通视觉 27.6，全项 PASS）
###############################################################################

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

# ---- 调色板（dataviz 参考实例，已验证）----
BLUE    = "#2a78d6"   # 核心 R 流程
ORANGE  = "#eb6834"   # 上游 cellranger
AQUA    = "#1baf7a"   # 高级模块
SURFACE = "#fcfcfb"
PILL_BG = "#f0efec"   # 数据产物底色（中性灰）
PILL_BD = "#c3c2b7"
INK     = "#0b0b0b"   # 主墨水
INK2    = "#52514e"   # 次级墨水
MUTED   = "#898781"   # 弱化（箭头/带标签）
HAIR    = "#e1e0d9"   # 发丝线

fig, ax = plt.subplots(figsize=(15.5, 7.75), dpi=100)
fig.patch.set_facecolor(SURFACE)
ax.set_xlim(0, 100); ax.set_ylim(0, 50)
ax.set_aspect("auto"); ax.axis("off")
fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

# ---- 通用绘制函数 ----
def band(x, y, w, h, hue, alpha=0.055):
    """阶段带容器：极浅色洗底 + 发丝边框"""
    ax.add_patch(FancyBboxPatch((x, y), w, h,
        boxstyle="round,pad=0,rounding_size=0.9",
        facecolor=hue, alpha=alpha, edgecolor=HAIR, linewidth=1, zorder=1))

def box(x, y, w, h, hue, title, subs, num=None, dashed=False):
    """阶段框：白底 + 类别色边框 + 骑线编号徽章（避免压住标题）+ 墨水色文字"""
    ax.add_patch(FancyBboxPatch((x, y), w, h,
        boxstyle="round,pad=0,rounding_size=0.7",
        facecolor="white", edgecolor=hue, linewidth=1.8,
        linestyle=(0, (4, 2)) if dashed else "solid", zorder=3))
    cx = x + w / 2
    # 文字纵向布局：标题在上，副文字按行距下排（保底边距，防溢出）
    ty = y + h * 0.60 if len(subs) <= 2 else y + h * 0.68
    ax.text(cx, ty, title, ha="center", va="center",
            fontsize=11.5, fontweight="bold", color=INK, zorder=4)
    y0 = ty - 2.1
    for i, s in enumerate(subs):
        ax.text(cx, y0 - i * 1.5, s, ha="center", va="center",
                fontsize=8.6, color=INK2, zorder=4)
    if num:
        # 徽章骑在框的上边线上（半内半外），不与任何文字重叠
        bw, bh = 5.2, 2.6
        bx, by = x + 1.6, y + h - bh / 2
        ax.add_patch(FancyBboxPatch((bx, by), bw, bh,
            boxstyle="round,pad=0,rounding_size=0.6",
            facecolor=hue, edgecolor=SURFACE, linewidth=1.2, zorder=5))
        ax.text(bx + bw / 2, by + bh / 2, num, ha="center", va="center",
                fontsize=8.4, fontweight="bold", color="white", zorder=6)

def pill(x, y, w, h, title, sub=None):
    """数据产物胶囊：中性灰底 + 灰边框"""
    ax.add_patch(FancyBboxPatch((x, y), w, h,
        boxstyle="round,pad=0,rounding_size=1.4",
        facecolor=PILL_BG, edgecolor=PILL_BD, linewidth=1.4, zorder=3))
    cy = y + h / 2 + (0.9 if sub else 0)
    ax.text(x + w / 2, cy, title, ha="center", va="center",
            fontsize=10.5, fontweight="bold", color=INK, zorder=4)
    if sub:
        ax.text(x + w / 2, y + h / 2 - 1.4, sub, ha="center", va="center",
                fontsize=8.2, color=INK2, zorder=4)

def arrow(p0, p1, dashed=False, color=MUTED):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>",
        mutation_scale=13, linewidth=1.4, color=color,
        linestyle=(0, (4, 2)) if dashed else "solid",
        shrinkA=0, shrinkB=0, zorder=2))

def elbow(pts, dashed=False, color=MUTED):
    """折线箭头：pts 为拐点序列，最后一段带箭头"""
    for a, b in zip(pts[:-2], pts[1:-1]):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, linewidth=1.4,
                linestyle=(0, (4, 2)) if dashed else "solid", zorder=2)
    arrow(pts[-2], pts[-1], dashed=dashed, color=color)

def band_label(x, y, text):
    ax.text(x, y, text, ha="left", va="center",
            fontsize=8.0, fontweight="bold", color=MUTED, zorder=4)

# ---- 左侧：config 面板 ----
ax.add_patch(FancyBboxPatch((1, 6), 15.5, 38,
    boxstyle="round,pad=0,rounding_size=0.9",
    facecolor="white", edgecolor=PILL_BD, linewidth=1.2, zorder=3))
ax.text(8.75, 41.2, "config/", ha="center", va="center",
        fontsize=12, fontweight="bold", color=INK, zorder=4)
ax.text(8.75, 39.2, "centralized parameters\nzero hardcoded in scripts",
        ha="center", va="center", fontsize=7.6, color=INK2, zorder=4)
cfg_items = [
    ("config.yaml", "single-sample profile"),
    ("config.multi.yaml", "multi-sample profile"),
    ("cellranger.yaml", "tool-level parameters"),
    ("annotation_presets/", "marker library (by species)"),
    ("trajectory_presets/", "lineage priors (root + subset)"),
]
for i, (name, desc) in enumerate(cfg_items):
    yy = 35.2 - i * 5.6
    ax.text(3.0, yy, name, ha="left", va="center",
            fontsize=8.4, fontweight="bold", color=INK, zorder=4,
            fontfamily="DejaVu Sans Mono")
    ax.text(3.0, yy - 1.7, desc, ha="left", va="center",
            fontsize=7.2, color=INK2, zorder=4)

# ---- 带 A：上游（cellranger）----
band(19, 37.5, 80, 11.5, ORANGE)
band_label(20.5, 47.6, "UPSTREAM — SHELL")
pill(20.5, 39.2, 8.5, 6, "FASTQ", "10x raw reads")
box(33, 39.2, 13.5, 6, ORANGE, "cellranger count",
    ["alignment · UMI counting"], num="00")
pill(51, 39.2, 10.5, 6, "Expression matrix", "barcodes·features·mtx")

# ---- 带 B：核心 R 流程 ----
band(19, 19, 80, 14.5, BLUE)
band_label(20.5, 31.9, "CORE PIPELINE — R / SEURAT V5")
bw, gap = 14.4, 1.5
bx = [20.5 + i * (bw + gap) for i in range(5)]
by, bh = 20.4, 9.6
box(bx[0], by, bw, bh, BLUE, "Load + QC",
    ["Read10X · merge", "QC filtering", "scDblFinder"], num="01")
box(bx[1], by, bw, bh, BLUE, "Preprocess + Cluster",
    ["Normalize · HVG · PCA", "harmony / RPCA / FastMNN", "UMAP · clustering"], num="02")
box(bx[2], by, bw, bh, BLUE, "Annotate",
    ["FindAllMarkers", "marker presets", "evidence check"], num="03")
box(bx[3], by, bw, bh, BLUE, "Multi-group plot",
    ["cell proportions", "grouped DotPlot", "marker heatmap"], num="04")
box(bx[4], by, bw, bh, BLUE, "DE + GSEA",
    ["FindMarkers", "compareCluster KEGG", "GSEA (fgsea)"], num="05")

# ---- 带 C：高级模块（可选）----
band(19, 3, 80, 12.5, AQUA)
band_label(20.5, 14.2, "ADVANCED MODULES — run_pipeline.sh --with-advanced")
cw = 24.5
cx = [20.5, 20.5 + cw + 3, 20.5 + 2 * (cw + 3)]
cy, ch = 4.2, 8
box(cx[0], cy, cw, ch, AQUA, "Gene-set scoring",
    ["AddModuleScore · pathway activity", "per-signature FeaturePlot"], num="06", dashed=True)
box(cx[1], cy, cw, ch, AQUA, "Trajectory — monocle3",
    ["pseudotime · lineage presets", "root cell type · raw / harmony / umap"], num="07", dashed=True)
box(cx[2], cy, cw, ch, AQUA, "Cell communication",
    ["CellChat · ligand–receptor DB", "STIM vs CTRL comparison"], num="08", dashed=True)

# ---- 箭头 ----
cy_a = 42.2   # 带 A 中线
arrow((29.0, cy_a), (33.0, cy_a))
arrow((46.5, cy_a), (51.0, cy_a))
# 矩阵 → 01（折线走 config 面板与带之间的空隙，从 01 左侧边进入，避免压住带标签）
elbow([(56.25, 39.2), (56.25, 36.2), (18.2, 36.2), (18.2, 25.2), (20.5, 25.2)])
cy_b = 25.2   # 带 B 中线
for i in range(4):
    arrow((bx[i] + bw, cy_b), (bx[i + 1], cy_b))
# 03 → 高级模块三扇（虚线 = 可选；止于带 C 上边线，不伸入带内压住标签）
fan_x = [cx[0] + cw / 2, cx[1] + cw / 2, cx[2] + cw / 2]
for tx in fan_x:
    arrow((bx[2] + bw / 2, by - 0.2), (tx, 15.7), dashed=True)

# ---- config 虚线连接（参数流入各带；B 带锚点避开矩阵折线的纵向段）----
for yy in (43.3, 21.5, 9.3):
    ax.plot([16.5, 19], [yy, yy], color=MUTED, linewidth=1.1,
            linestyle=(0, (3, 2)), zorder=2)

# ---- 底部：批次归档说明条 ----
ax.text(59, 1.2,
        "Every run archives to  output/<batch>/<step>/  —  result .rds + figures/ + logs/ + .done marker  ·  "
        "config_used.yaml parameter snapshot  ·  same batch = resume, new batch = fresh run",
        ha="center", va="center", fontsize=8.2, color=INK2, zorder=4)

# ---- 保存 ----
out_dir = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(out_dir, "pipeline_overview.png"), dpi=200,
            facecolor=SURFACE, bbox_inches="tight", pad_inches=0.05)
fig.savefig(os.path.join(out_dir, "pipeline_overview.svg"),
            facecolor=SURFACE, bbox_inches="tight", pad_inches=0.05)
print("saved: docs/assets/pipeline_overview.{png,svg}")
