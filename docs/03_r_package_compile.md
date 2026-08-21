# conda 环境 R 包源码编译踩坑

> 从 CRAN / GitHub 源码安装 R 包时的两个高频坑：编译器代际不匹配、R 4.5 头文件变动。
> 均实测踩坑（2026-08，环境：mamba 环境 `scrna`，R 4.5.3 + conda-forge gcc 15.1.0）。

## 1. 编译器代际匹配

- conda 环境的 R 必须用 conda 自带编译器套件（`envs/scrna/bin/` 下的 `x86_64-conda-linux-gnu-gcc/g++`），软链接 `x86_64-conda-linux-gnu-cc` / `-cxx` 也必须指向该套件
- 实测教训：软链接指向系统 `/usr/bin/gcc` 11 时，R 4.5.3 头文件的 C23 `enum :int` 语法全部编译失败——R 4.5 是用 gcc 13+ 构建的，`Rconfig.h` 无条件开启 `HAVE_ENUM_BASE_TYPE`，gcc 11 解析不了
- 若报 `-std=gnu23` 相关错误，把 `Makeconf` 中 `-std=gnu23` 改为 `-std=gnu17`（gcc 15 同样接受）

## 2. Rboolean 未定义（老包编译失败）

- 症状：老包（NMF / foreign / urca 等）直接 `#include <Rinternals.h>` 报 `Rboolean` 未定义——R 4.5 把它移到了 `R.h`
- conda 的 R 不认环境变量 `CFLAGS`，需临时写 `~/.R/Makevars`：

```make
PKG_CFLAGS += -include R.h
PKG_CXXFLAGS += -include R.h
```

- 装完删除 `~/.R/Makevars`，避免影响后续正常安装
- 实测：此方案 + gcc 15 成功安装 CellChat 2.2.0.9001（GitHub，CRAN 已下架）+ NMF 0.28
