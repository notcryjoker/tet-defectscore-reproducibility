#!/usr/bin/env Rscript
# =============================================================================
# 09_fig1_design.R —— Fig 1 研究设计流程图（不依赖分析结果，可随时生成）
#   版式: 三列式流程（数据→方法→分析），符合期刊 Figure 1 惯例
# =============================================================================
.ofile <- sys.frame(1)$ofile
if (is.null(.ofile) || !nzchar(.ofile)) .ofile <- "."
.cfg_cands <- c(
  file.path(Sys.getenv("TET_ROOT", "."), "code", "00_config.R"),
  file.path(dirname(.ofile), "00_config.R"),
  "code/00_config.R", "00_config.R"
)
.cfg <- .cfg_cands[sapply(.cfg_cands, file.exists)]
if (length(.cfg) == 0L) stop("00_config.R not found; set TET_ROOT or run from deposit root") else source(.cfg[1])
load_libs(c("ggplot2", "grid"))

out <- file.path(DIR_FIG, "fig1_study_design.png")
png(out, width = 3400, height = 2000, res = 300, bg = "white")

x_box <- c(0.06, 0.40, 0.74)   # 三列中心
w     <- 0.27                  # 盒宽
h     <- 0.135                 # 盒高

box <- function(x, y, w, h, fill, label, sub = NULL, fs = 11, subfs = 9) {
  grid::grid.roundrect(
    x = x, y = y, width = w, height = h,
    r = grid::unit(2, "mm"),
    gp = grid::gpar(fill = fill, col = "white", lwd = 1.2, alpha = 0.92)
  )
  if (!is.null(sub)) {
    grid::grid.text(label, x = x, y = y + 0.018, gp = gpar(fontsize = fs, fontface = "bold", col = "white"))
    grid::grid.text(sub,   x = x, y = y - 0.020, gp = gpar(fontsize = subfs, col = "white"))
  } else {
    grid::grid.text(label, x = x, y = y, gp = gpar(fontsize = fs, fontface = "bold", col = "white"))
  }
}

arrow <- function(x1, y1, x2, y2) {
  grid::grid.segments(
    x0 = x1, y0 = y1, x1 = x2, y1 = y2,
    gp = gpar(lwd = 2.4, col = "#54606C", arrow = grid::arrow(length = unit(3, "mm"), type = "closed"))
  )
}

grid.newpage()
pushViewport(viewport(width = 0.98, height = 0.94, x = 0.5, y = 0.53))

# ---- 列标题 -------------------------------------------------------------------
cols <- c("① 数据来源", "② 预锁定方法学", "③ 分析与验证")
col_cols <- c("#1F6FB2", "#0E8074", "#B2561F")
for (i in 1:3)
  grid::grid.text(cols[i], x = x_box[i], y = 0.955,
                  gp = gpar(fontsize = 15, fontface = "bold", col = col_cols[i]))

# ---- ① 数据来源 ----------------------------------------------------------------
box(x_box[1], 0.815, w, h, "#2E86C1", "TCGA-THYM（发现队列）", "bulk RNA-seq, n = 120")
box(x_box[1], 0.625, w, h, "#2E86C1", "GSE29695（外部验证）", "DASL 芯片, n = 37")
box(x_box[1], 0.435, w, h, "#5DADE2", "单细胞参考图谱", "GSE147520 / GSE284479")
box(x_box[1], 0.245, w, h, "#5DADE2", "GSE181815（ICB 探索）", "thymic carcinoma, n = 10")

# ---- ② 方法 --------------------------------------------------------------------
box(x_box[2], 0.815, w, h, "#148F77", "胸腺发育 9 状态基因集", "ETP/DN/DP/SP/Treg/mTEC/cTEC/TRA")
box(x_box[2], 0.625, w, h, "#148F77", "复合评分定义", "CORTICAL（成熟阻滞）vs MEDULLARY（阴性选择完整性）")
box(x_box[2], 0.435, w, h, "#1ABC9C", "DefectScore", "z(CORTICAL) − z(MEDULLARY)，连续可解释")
box(x_box[2], 0.245, w, h, "#1ABC9C", "候选基因预锁定", "关联检验前固定，防数据挖掘")

# ---- ③ 分析与验证 ----------------------------------------------------------------
box(x_box[3], 0.815, w, h, "#CA6F1E", "临床关联", "WHO 谱系 · MG（34/86）· Masaoka")
box(x_box[3], 0.625, w, h, "#CA6F1E", "生存分析", "RFS（复发 13 事件；OS 事件过少不建模）")
box(x_box[3], 0.435, w, h, "#E67E22", "外部验证", "GSE29695 复现 WHO 梯度与 MG 关联")
box(x_box[3], 0.245, w, h, "#E67E22", "探索层", "ICB 反应关联（假设生成）")

# ---- 箭头：数据 → 方法 ----------------------------------------------------------
arrow(x_box[1] + w/2, 0.815, x_box[2] - w/2, 0.815)   # TCGA -> gene sets
arrow(x_box[1] + w/2, 0.625, x_box[2] - w/2, 0.625)   # GSE29695 -> validation plan
arrow(x_box[1] + w/2, 0.435, x_box[2] - w/2, 0.435)   # scRNA refs -> gene sets
arrow(x_box[1] + w/2, 0.245, x_box[3] - w/2, 0.245)   # ICB data -> exploratory

# ---- 箭头：方法内部 -------------------------------------------------------------
arrow(x_box[2], 0.815 - h/2, x_box[2], 0.625 + h/2)
arrow(x_box[2], 0.625 - h/2, x_box[2], 0.435 + h/2)
arrow(x_box[2], 0.435 - h/2, x_box[2], 0.245 + h/2)

# ---- 箭头：方法 → 分析 ----------------------------------------------------------
arrow(x_box[2] + w/2, 0.435, x_box[3] - w/2, 0.815)   # score -> associations
arrow(x_box[2] + w/2, 0.435, x_box[3] - w/2, 0.625)   # score -> survival

# ---- 底部方法学徽标 -------------------------------------------------------------
grid::grid.text(
  "方法学保障：候选基因预锁定 · 外部验证 · 可解释评分（非黑箱 ML）· 关联而非因果",
  x = 0.5, y = 0.075,
  gp = gpar(fontsize = 10, fontface = "italic", col = "#566573")
)
grid::popViewport()
dev.off()
cat("[fig1] saved:", out, "\n")
