#!/usr/bin/env Rscript
# =============================================================================
# 00_config.R  —— 全局配置（路径 / 随机种子 / 库加载）
# 所有下游脚本 source 本文件
# =============================================================================
ROOT  <- Sys.getenv("TET_ROOT", normalizePath(".", mustWork = FALSE))  # 解压目录；可用环境变量 TET_ROOT 覆盖
DIR_RAW     <- file.path(ROOT, "data/raw")
DIR_PROC    <- file.path(ROOT, "data/processed")
DIR_CODE    <- file.path(ROOT, "code")
DIR_RES     <- file.path(ROOT, "results")
DIR_FIG     <- file.path(ROOT, "figures")
for (d in c(DIR_PROC, DIR_FIG)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

SEED <- 20261022
set.seed(SEED)

# 期刊友好型 300 dpi 出版级默认
FIG_DPI  <- 300
FIG_W    <- 8.27   # A4 宽 inch
FIG_H    <- 5.83

# 库加载助手（失败即报错，避免静默缺失）
load_libs <- function(pkgs) {
  for (p in pkgs) {
    if (!suppressWarnings(requireNamespace(p, quietly = TRUE)))
      stop(sprintf("缺少必要 R 包: %s —— 请先安装", p))
    library(p, character.only = TRUE)
  }
}
cat("[config] ROOT =", ROOT, "\n")
