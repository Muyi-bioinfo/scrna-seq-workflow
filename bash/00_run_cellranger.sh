#!/usr/bin/env bash
###############################################################################
### 00_run_cellranger.sh — 模块 0：上游定量（cellranger count，批量 shell 版）
###
### 用法：bash bash/00_run_cellranger.sh [主config] [cellranger配置]
###   主config（默认 config/config.multi.yaml）：读 batch 字段与样本表路径
###   cellranger配置（默认 config/cellranger.yaml）：工具级参数（与下游 R 解耦）
### 输入：样本表的 fastqs 列；输出：output/<batch>/00_cellranger/<样本名>/outs/
### 讲解见 docs/02_cellranger_guide.md
###############################################################################

set -euo pipefail

MAIN_CONFIG=${1:-${CONFIG_FILE:-config/config.multi.yaml}}
CR_CONFIG=${2:-config/cellranger.yaml}
[[ -f "$MAIN_CONFIG" ]] || { echo "✘ 找不到主配置: $MAIN_CONFIG"; exit 1; }
[[ -f "$CR_CONFIG" ]]   || { echo "✘ 找不到 cellranger 配置: $CR_CONFIG"; exit 1; }

# ---- 主 config：batch（批次归档）+ 样本表路径（多样本段）----
BATCH=$(awk '$1 == "batch:" {print $2; exit}' "$MAIN_CONFIG")
[[ -n "$BATCH" ]] || { echo "✘ 主配置缺少 batch 字段"; exit 1; }
SAMPLE_SHEET=$(awk '/^multi:/ {in_sec=1; next}
                    in_sec && /^[a-z]/ {exit}
                    in_sec && $1 == "sample_sheet:" {print $2; exit}' "$MAIN_CONFIG")
[[ -n "$SAMPLE_SHEET" && -f "$SAMPLE_SHEET" ]] || { echo "✘ 找不到样本表: $SAMPLE_SHEET"; exit 1; }

# ---- cellranger 配置：全部是顶层键，直接取 ----
cr_get() { awk -v key="$1" '$1 == key":" {print $2; exit}' "$CR_CONFIG"; }
TRANSCRIPTOME=$(cr_get transcriptome)
EXPECT_CELLS=$(cr_get expect_cells)
CREATE_BAM=$(cr_get create_bam | tr '[:upper:]' '[:lower:]')
LOCALCORES=$(cr_get localcores)
LOCALMEM=$(cr_get localmem)
[[ -d "$TRANSCRIPTOME" ]] || {
  echo "✘ 参考基因组不存在: $TRANSCRIPTOME（下载解压到 data/reference/refdata/ 后在 cellranger.yaml 里填）"
  exit 1
}
# ⚠️ 转绝对路径：校验在项目根目录通过，但使用时已 cd 到输出目录，相对路径会失效
TRANSCRIPTOME=$(realpath "$TRANSCRIPTOME")

# ---- 环境自检 ----
command -v cellranger >/dev/null || { echo "✘ 找不到 cellranger，请确认已加入 PATH"; exit 1; }

# ---- 输出目录（批次归档 + 配置快照，可追溯）----
OUTDIR="output/$BATCH/00_cellranger"
mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")   # ⚠️ 同上：logfile 重定向在 cd 之后的子 shell 中，相对路径会重复拼接
cp "$CR_CONFIG" "$OUTDIR/cellranger_config_used.yaml"
echo "批次: $BATCH（定量结果 → $OUTDIR/）"

# ---- 读样本表（sample, group, fastqs, matrix），对 fastqs 非空的样本跑 count ----
# ⚠️ tr -d '\r' 兼容 Windows 编辑的 CSV（CRLF 会让路径带上回车导致目录检查失败）
tail -n +2 "$SAMPLE_SHEET" | tr -d '\r' | while IFS=',' read -r sample group fastqs matrix; do
  [[ -n "${fastqs:-}" ]] || continue
  [[ -d "$fastqs" ]] || { echo "✘ 样本 $sample 的 FASTQ 目录不存在: $fastqs"; exit 1; }
  fastqs=$(realpath "$fastqs")   # 同上：cd 到输出目录后相对路径失效

  echo ""
  echo "▶ ▶ ▶ cellranger count | 样本: $sample（fastqs: $fastqs）"
  logfile="$OUTDIR/cellranger_${sample}.log"

  # cellranger 的产出固定在当前目录的 <id>/ 下，切到输出目录运行
  (
    cd "$OUTDIR"
    cellranger count \
      --id="$sample" \
      --transcriptome="$TRANSCRIPTOME" \
      --fastqs="$fastqs" \
      --expect-cells="$EXPECT_CELLS" \
      --create-bam="$CREATE_BAM" \
      --localcores="$LOCALCORES" \
      --localmem="$LOCALMEM" \
      > "$logfile" 2>&1
  )
  echo "✔ $sample 定量完成（质检报告: $OUTDIR/$sample/outs/web_summary.html）"
done

echo ""
echo "===== 上游定量完成 → 下游 01_load_qc 会自动读取各样本矩阵 ====="
