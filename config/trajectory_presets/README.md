# 轨迹预设库（trajectory_presets）

## 什么是预设

07_trajectory_analysis（monocle3 拟时序）需要**谱系先验**：哪些细胞构成生物学上的连续过程、谁是起点（根）。数据本身推不出这个先验，必须由分析者指定。预设就是把这份先验显式化、文件化，跨项目积累复用。

主 config 里用 `advanced$trajectory$preset: <文件名>` 引用；留空 = 全部细胞跑（无谱系先验，仅保证流程跑通）。

## 字段说明

| 字段 | 必填 | 说明 |
|------|------|------|
| `description` | 是 | 一句话描述轨迹，打印到运行日志 |
| `root_celltype` | 是 | 根（起点）：祖细胞 / naive / 未激活状态。留空则自动用最大细胞群（无生物学意义） |
| `subset_celltypes` | 否 | 参与轨迹的细胞类型列表，支持精确名与 `*`/`?` 前缀模式（如 `CD4*`） |
| `num_dim` | 否 | monocle3 预处理维度数，覆盖主 config |

**优先级**：主 config 显式参数 > 预设 > 默认值。

## 匹配与校验规则

- 精确名：与注释（`cluster_map`）**完全一致**才能匹配
- 通配模式：`CD4*` 匹配所有 CD4 开头的类型，日志会打印实际匹配结果
- 全部失配 → 报错并打印数据中全部可用 celltype 列表（照着列表改名字即可）
- 每次运行会在步骤目录写 `trajectory_resolved.yaml` 快照，记录本次实际生效的先验

## 怎么写一个新预设（checklist）

1. **确定连续过程**：看 `03_annotate` 输出的 marker 表（`all_markers_table.rds` / `top10_markers.csv`），确定哪些细胞群构成生物学上的连续分化/激活过程
2. **定起点**：哪个状态是"共同祖先"（祖细胞 / naive / 未激活）——写进 `root_celltype`
3. **查依据**：找一篇代表文献，把引用和一句生物学依据写进文件头注释（`notes` 位置见模板）
4. **对齐命名**：细胞名与你的注释完全一致，或用 `*` 前缀模式；拿不准就故意写错跑一次，脚本会"自报家底"
5. **接入**：主 config 写 `trajectory$preset: <文件名>`，跑通验证后再留在库里

原则：**先验证，后入库**。没跑通过的预设不要放进库里。

## 现有预设

| 文件 | 谱系 | 根 | 细胞 | 依据 |
|------|------|----|------|------|
| `pbmc_T_cell.yaml` | T 细胞 | CD4 Naive T | CD4 Naive/Memory T、T activated、CD8 T | Szabo et al. 2019, Nat Commun |
| `pbmc_myeloid.yaml` | 髓系 | CD14 Mono | CD14 Mono、CD16 Mono | 经典→非经典单核细胞分化轴 |

`_template.yaml` 是写新预设的填空模板。
