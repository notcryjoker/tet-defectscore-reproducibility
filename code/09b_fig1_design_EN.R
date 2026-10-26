#!/usr/bin/env Rscript
# =============================================================================
# 09b_fig1_design_EN.R —— Fig 1 研究设计流程图英文版（投稿主稿用）
#   版式与中文版一致：三列式流程（数据→方法→分析）
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

out <- file.path(DIR_FIG, "fig1_study_design_EN.png")
png(out, width = 3400, height = 2000, res = 300, bg = "white")

x_box <- c(0.17, 0.47, 0.77)   # 三列中心（右移避免左侧裁切）
w     <- 0.28                  # 盒宽
h     <- 0.135                 # 盒高

box <- function(x, y, w, h, fill, label, sub = NULL, fs = 10, subfs = 8.5) {
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

# ---- Column headers -------------------------------------------------------------
cols <- c("① Data sources", "② Pre-locked methodology", "③ Analysis & validation")
col_cols <- c("#1F6FB2", "#0E8074", "#B2561F")
for (i in 1:3)
  grid::grid.text(cols[i], x = x_box[i], y = 0.955,
                  gp = gpar(fontsize = 15, fontface = "bold", col = col_cols[i]))

# ---- ① Data sources -------------------------------------------------------------
box(x_box[1], 0.815, w, h, "#2E86C1", "TCGA-THYM (discovery)", "bulk RNA-seq, n = 120")
box(x_box[1], 0.625, w, h, "#2E86C1", "GSE29695 (external)", "DASL array, n = 37")
box(x_box[1], 0.435, w, h, "#5DADE2", "Single-cell reference atlases", "GSE147520 / GSE284479")
box(x_box[1], 0.245, w, h, "#5DADE2", "GSE181815 (ICB exploration)", "thymic carcinoma, n = 10")

# ---- ② Methods ------------------------------------------------------------------
box(x_box[2], 0.815, w, h, "#148F77", "Nine thymopoiesis cell-state gene sets", "ETP/DN/DP/SP/Treg/mTEC/cTEC/TRA")
box(x_box[2], 0.625, w, h, "#148F77", "Composite score definition", "CORTICAL (arrest) vs MEDULLARY (negative selection)")
box(x_box[2], 0.435, w, h, "#1ABC9C", "DefectScore", "z(CORTICAL) − z(MEDULLARY); interpretable")
box(x_box[2], 0.245, w, h, "#1ABC9C", "Candidate pre-specification", "Fixed before any association test")

# ---- ③ Analysis & validation ----------------------------------------------------
box(x_box[3], 0.815, w, h, "#CA6F1E", "Clinical associations", "WHO histology · MG (34/86) · Masaoka")
box(x_box[3], 0.625, w, h, "#CA6F1E", "Survival analysis", "RFS (13 relapse events; OS not modeled)")
box(x_box[3], 0.435, w, h, "#E67E22", "External validation", "GSE29695: WHO gradient and MG link")
box(x_box[3], 0.245, w, h, "#E67E22", "Exploratory layer", "ICB response (hypothesis-generating)")

# ---- Arrows: data -> methods ----------------------------------------------------
arrow(x_box[1] + w/2, 0.815, x_box[2] - w/2, 0.815)   # TCGA -> gene sets
arrow(x_box[1] + w/2, 0.625, x_box[2] - w/2, 0.625)   # GSE29695 -> validation plan
arrow(x_box[1] + w/2, 0.435, x_box[2] - w/2, 0.435)   # scRNA refs -> gene sets
arrow(x_box[1] + w/2, 0.245, x_box[3] - w/2, 0.245)   # ICB data -> exploratory

# ---- Arrows: within methods -----------------------------------------------------
arrow(x_box[2], 0.815 - h/2, x_box[2], 0.625 + h/2)
arrow(x_box[2], 0.625 - h/2, x_box[2], 0.435 + h/2)
arrow(x_box[2], 0.435 - h/2, x_box[2], 0.245 + h/2)

# ---- Arrows: methods -> analysis ------------------------------------------------
arrow(x_box[2] + w/2, 0.435, x_box[3] - w/2, 0.815)   # score -> associations
arrow(x_box[2] + w/2, 0.435, x_box[3] - w/2, 0.625)   # score -> survival

# ---- Footer ---------------------------------------------------------------------
grid::grid.text(
  "Methodological safeguards: candidate pre-specification · external validation · interpretable scores (no black-box ML) · association, not causation",
  x = 0.5, y = 0.075,
  gp = gpar(fontsize = 10, fontface = "italic", col = "#566573")
)
grid::popViewport()
dev.off()
cat("[fig1-EN] saved:", out, "\n")
