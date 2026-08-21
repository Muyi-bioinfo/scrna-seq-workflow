# mouse/ — 小鼠注释预设

> 预留目录：mouse 预设的写法与 human/ 完全一致（见根目录 `_template.yaml`），
> 只是 species 字段为 `mouse`、基因符号按小鼠习惯（首字母大写，如 `Cd3d`）。

## 待建预设素材（来自《单细胞分析中常见marker基因汇总》，未验证，建库时逐条核对）

### 小鼠特有基因（人源无对应物，注意不要混入 human/ 预设）

| 基因 | 含义 | 用途 |
|------|------|------|
| Ly6c2 | Ly6C 单核标志 | 小鼠单核亚群 |
| Siglech | Siglec-H | 小鼠 pDC |
| F4/80（Adgre1） | 巨噬标志 | 小鼠巨噬/Kupffer |

### 小鼠髓系（原表小写基因，转小鼠规范后）

- Macrophages: Adgre1, Cd14, Fcgr3, Lyz2, Cd68, Cd163
- cDCs: Xcr1, Flt3, Ccr7, Cd1e
- pDCs: Siglech, Clec10a, Clec12a
- Monocytes: Ly6c2, Spn, Cd300e
- Neutrophils: Csf3r, S100a8, Cxcl3

## 已知坑（建库时注意）

1. 原 PDF 神经系统的 Pericyte 列表被成纤维 marker 污染（DCN/LUM/GSN/FGF7/MME 是成纤维的）
2. 原 PDF Endothelial 列表中的 EBF1 是 B 细胞转录因子，非内皮 marker，疑似抄错
3. 原 PDF Inhibitory neurons 的 PCDH15 是内耳毛细胞基因，存疑
