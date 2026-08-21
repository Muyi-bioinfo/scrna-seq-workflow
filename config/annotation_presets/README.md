# config/annotation_presets/ — 注释 marker 预设库

> 细胞注释的生物学知识库：**细胞类型 → marker 基因**（跨数据集可复用的部分）。
> cluster 编号 → 细胞类型的映射因数据集/分辨率而异，仍由人工在 config 的
> `annotate$cluster_map` 填写——预设只提供 marker 与类型名一致性校验。

## 目录结构

```
annotation_presets/
├── README.md          # 本文件
├── _template.yaml     # 空白模板（新建预设照此填）
├── human/             # 人源预设（基因符号全大写）
│   └── pbmc.yaml      # 外周血 PBMC
└── mouse/             # 小鼠预设（预留，见 mouse/README.md）
```

## 用法

config 的 annotate 段引用：

```yaml
annotate:
  preset: human/pbmc            # 文件夹/文件名（不含 .yaml）
  # marker_genes: [MS4A1, ...]  # 显式扁平列表（优先于 preset，一般不用）
```

03 的可视化为双视图（互为补充）：
- `dotplot_markers_by_cluster.pdf`：tier1 marker × cluster 矩阵，判断「每个 cluster 是谁」
- `figures/featureplot_celltypes/featureplot_<类型名>.pdf`：每个类型一张
  （单独子目录，ncol=2），验证 marker 的空间一致性

加载时自动校验（见 utils.R 的 `load_annotation_preset()`）：
1. **物种一致性**：preset 的 species 字段必须与所在文件夹一致
2. **基因存在性**：marker 在数据中缺失超过 50% 时 warning（人鼠混用/数据集不适用会立刻暴露）
3. **类型名一致性**：cluster_map 的类型名命中 preset 别名（aliases）时提示标准名，完全未知时 warning

## 写作清单（新建预设）

1. 复制 `_template.yaml` 到对应物种子文件夹，以组织名命名（如 `human/liver.yaml`）
2. 填 description / tissue / species / source / version——source 必须可追溯
3. marker 分级：tier1 只放公认经典（判断依据）；文献特异差异基因不要进 tier1
4. 基因符号与物种匹配（人全大写 / 鼠首字母大写），**禁止人鼠混用**（坑见 mouse/README.md）
5. 有公认阴性 marker 时填 negative（如浆细胞的 MS4A1）
6. 常见同义名填 aliases（如 "CD16 Mono" → "FCGR3A+ Mono"）
7. 自测：`Rscript` 跑一次 03 看校验无 warning

## 现有预设

| 预设 | 物种 | 组织 | 覆盖 |
|------|------|------|------|
| human/pbmc | human | 外周血 PBMC | 14 类（含活化 T/B、pDC、红细胞等 IFN 刺激数据亚群） |
| human/liver | human | 肝脏 | 7 类（肝细胞/胆管/Kupffer/星状等） |
| human/heart | human | 心脏 | 7 类（心肌/SMC/成纤维/免疫等） |
| human/brain | human | 脑 | 9 类（神经元/胶质/周细胞等） |
| human/gastric_tumor | human | 胃癌 | 8 类（上皮/免疫/基质） |
| human/hcc | human | 肝细胞癌 | 7 类（肿瘤/免疫/基质） |

## 待办

- mouse 髓系/组织预设（素材见 mouse/README.md）
- ⚠️ pbmc 以外预设的 marker 均来自汇总 PDF 人工核查，**尚未在真实数据上验证**，
  使用前建议先在目标数据跑 03 看 marker 存在率与 FeaturePlot
