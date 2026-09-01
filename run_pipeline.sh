#!/usr/bin/env bash
###############################################################################
### run_pipeline.sh — 全流程总控（顺序执行 + 断点续跑 + 每步日志）
###
### 用法：
###   bash run_pipeline.sh                          # 单样本模式全流程（默认）
###   bash run_pipeline.sh --mode multi             # 多样本模式全流程
###   bash run_pipeline.sh --mode multi --with-advanced   # 追加高级模块
###   bash run_pipeline.sh --only 04                # 只跑第 4 步
###   bash run_pipeline.sh --from 03 --to 06        # 跑 03~06
###   bash run_pipeline.sh --force                  # 已有产物的步骤也强制重跑
###   bash run_pipeline.sh --config config/xxx.yaml # 自定义配置文件
###                                                # （默认: single→config.yaml, multi→config.multi.yaml）
###
### 批次与断点续跑：
###   - 每次运行归档到 output/<batch>/：各步骤结果 + logs/ + config_used.yaml 参数快照
###   - 同批次名重跑 = 断点续跑（output/<batch>/<step>/.done 标记）；换批次名 = 全新运行
###############################################################################

set -euo pipefail

RSCRIPT=${RSCRIPT:-Rscript}          # 可指定 Rscript 路径（HPC 节点上可用绝对路径）
MODE="single"
FROM=""
TO=""
ONLY=""
FORCE=false
WITH_ADVANCED=false
CONFIG=""                             # 自定义配置文件（默认按模式自动选择）

# ---- 参数解析 ----
while [[ $# -gt 0 ]]; do
  case $1 in
    --mode)          MODE="$2"; shift 2 ;;
    --from)          FROM="$2"; shift 2 ;;
    --to)            TO="$2"; shift 2 ;;
    --only)          ONLY="$2"; shift 2 ;;
    --force)         FORCE=true; shift ;;
    --with-advanced) WITH_ADVANCED=true; shift ;;
    --config)        CONFIG="$2"; shift 2 ;;
    -h|--help)
      sed -n '3,20p' "$0" | sed 's/^### \{0,1\}//'
      exit 0 ;;
    *) echo "未知参数: $1（--help 查看用法）"; exit 1 ;;
  esac
done

# ---- 配置文件选择（传给 R 脚本的环境变量 CONFIG_FILE）----
if [[ -n "$CONFIG" ]]; then
  export CONFIG_FILE="$CONFIG"
else
  case $MODE in
    single) export CONFIG_FILE="config/config.yaml" ;;
    multi)  export CONFIG_FILE="config/config.multi.yaml" ;;
    *)      echo "未知模式: $MODE（支持 single / multi）"; exit 1 ;;
  esac
fi
[[ -f "$CONFIG_FILE" ]] || { echo "找不到配置文件: $CONFIG_FILE"; exit 1; }

# ---- 步骤清单（步骤名 → 脚本，对应 R/ 下同名脚本、output/ 下同名目录）----
# 脚本粒度：按分析阶段打包（读入+质控 / 预处理+整合+聚类），而非每个操作一个脚本
case $MODE in
  single)
    STEPS=("01_load_qc" "02_preprocess_cluster" "03_annotate") ;;
  multi)
    STEPS=("01_load_qc" "02_preprocess_cluster" "03_annotate" "04_multi_group_plot" "05_diff_gsea" "06_pseudobulk_de") ;;
  *)
    echo "未知模式: $MODE（支持 single / multi）"; exit 1 ;;
esac
if [[ "$WITH_ADVANCED" == true ]]; then
  STEPS+=("07_gene_set_score" "08_trajectory_analysis" "09_cell_communication")
fi

# ---- 批次（每次运行完全自包含归档到 output/<batch>/）----
# ⚠️ set -euo pipefail 下 grep 无匹配会让赋值命令直接失败，脚本会在打印
# 任何信息前静默退出——先 || true 兜底，再显式判断
BATCH=$(grep -E '^batch:' "$CONFIG_FILE" | head -1 | awk '{print $2}' || true)
if [[ -z "$BATCH" ]]; then
  echo "✘ config 缺少 batch 字段，请先在 $CONFIG_FILE 中填写"; exit 1
fi
# 结果根目录与 R 侧 batch_dir() 同源（config 的 dirs$base），不硬编码 output/
BASE=$(grep -E '^[[:space:]]*base:' "$CONFIG_FILE" | head -1 | awk '{print $2}' || true)
[[ -n "$BASE" ]] || BASE="output"
OUTBATCH="$BASE/$BATCH"
LOGDIR="$OUTBATCH/logs"

# ---- 环境自检 ----
command -v "$RSCRIPT" >/dev/null || {
  echo "✘ 找不到 $RSCRIPT —— 请先激活环境（mamba activate scrna）或 export RSCRIPT=/path/to/Rscript"; exit 1; }
echo "使用 Rscript: $($RSCRIPT --version 2>&1 | head -1)"

# ---- 执行 ----
mkdir -p "$LOGDIR" "$OUTBATCH"
cp "$CONFIG_FILE" "$OUTBATCH/config_used.yaml"   # 参数快照：这个批次是用什么参数跑出来的
echo "批次: $BATCH（结果与日志 → $OUTBATCH/）"

# ⚠️ BLAS 单线程固定（全步骤）：irlba（Seurat RunPCA、monocle3 preprocess）在
# OpenBLAS 多线程下段错误（实测 r-irlba 2.3.7 + OpenBLAS 0.3.33，最初在 07 拟时序
# 触发；2026-09-01 装 deseq2/apeglm 刷新环境后 02 的 RunPCA 也触发）。R 内
# Sys.setenv 无效，必须 shell 层设置；宁可损失一点多线程速度也不中途段错误
export OPENBLAS_NUM_THREADS=1

# ---- 步骤过滤参数格式校验：支持数字序号（04）或完整步骤名（04_multi_group_plot）----
for arg_name in FROM TO ONLY; do
  v="${!arg_name}"
  if [[ -n "$v" && ! "$v" =~ ^[0-9]+(_[a-z0-9_]+)?$ ]]; then
    echo "✘ $arg_name 参数格式错误: $v（支持数字序号如 04，或完整步骤名如 04_multi_group_plot）"
    exit 1
  fi
done

FAILED=""
ran_any=false          # 至少匹配到一个步骤才算有效运行（--only 打错字时不能静默假成功）
for step in "${STEPS[@]}"; do
  step_num="${step%%_*}"   # 04_multi_group_plot → 04

  # --from / --to：按数字序号比较
  # ⚠️ 词法比较有坑（"06_gene_set_score" > "06" 会把第 6 步自己滤掉）；
  # 用 10# 前缀避免 08 这类前导零被当八进制
  if [[ -n "$FROM" ]]; then
    from_num="${FROM%%_*}"
    (( 10#$step_num < 10#$from_num )) && continue
  fi
  if [[ -n "$TO" ]]; then
    to_num="${TO%%_*}"
    (( 10#$step_num > 10#$to_num )) && continue
  fi
  # --only：数字序号或完整步骤名均可匹配
  if [[ -n "$ONLY" && "$step" != "$ONLY" && "$step_num" != "$ONLY" ]]; then
    continue
  fi
  ran_any=true

  if [[ -f "$OUTBATCH/$step/.done" ]] && [[ "$FORCE" != true ]]; then
    echo "⏭  跳过 $step（已有产物，--force 强制重跑）"
    continue
  fi

  echo ""
  echo "▶ ▶ ▶ 运行 $step"
  if $RSCRIPT "R/${step}.R" > "$LOGDIR/${step}.log" 2>&1; then
    touch "$OUTBATCH/$step/.done"
    echo "✔ $step 完成（日志: $LOGDIR/${step}.log）"
  else
    echo "✘ $step 失败！查看日志: $LOGDIR/${step}.log（重跑: bash run_pipeline.sh --only $step）"
    FAILED="$step"
    break
  fi
done

# ---- 汇总 ----
echo ""
if [[ -z "$FAILED" ]]; then
  if [[ "$ran_any" != true ]]; then
    echo "===== ✘ 没有匹配到任何步骤（--only/--from/--to 与步骤清单不符），未执行任何分析 ====="
    echo "本模式的步骤清单: ${STEPS[*]}"
    exit 1
  fi
  echo "===== 全部步骤完成！结果见 $OUTBATCH/ ====="
else
  echo "===== 流程在 $FAILED 中断 ====="
  exit 1
fi
