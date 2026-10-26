#!/usr/bin/env Rscript
# =============================================================================
# 03b_scRNA.R —— 单细胞参考验证预锁定基因集（GSE147520 人胸腺 stromal scRNA）
#   输入: data/raw/GSE147520/scRNA/{barcodes,genes,matrix}.mtx.gz + thymo_signatures.rds
#   方法: 轻量 scater 流程（过滤→CPM→PCA→UMAP）→ 逐细胞 signature 评分
#   输出: fig2_scRNA_umap.png（UMAP 着色签名活性）+ fig3_scRNA_sig_panel.png
#   定位: 证明预锁定基因集在真实单细胞数据中非冗余、可分离（验证而非重新定义）
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
load_libs(c("Matrix", "scater", "data.table", "ggplot2", "uwot", "gridExtra", "viridis"))

sc_dir <- file.path(DIR_RAW, "GSE147520", "scRNA")
sig <- readRDS(file.path(DIR_CODE, "thymo_signatures.rds"))

read_mtx <- function(prefix) {
  mtx_f <- file.path(sc_dir, paste0(prefix, "_matrix.mtx.gz"))
  bc_f  <- file.path(sc_dir, paste0(prefix, "_barcodes.tsv.gz"))
  gf    <- file.path(sc_dir, paste0(prefix, "_features.tsv.gz"))
  if (!file.exists(gf)) gf <- file.path(sc_dir, paste0(prefix, "_genes.tsv.gz"))
  if (!all(file.exists(mtx_f, bc_f, gf))) {
    cat(sprintf("[scRNA] 缺少 %s 的关键文件，跳过\n", prefix)); return(NULL)
  }
  m   <- Matrix::readMM(mtx_f)
  bc  <- read.delim(bc_f, header = FALSE, stringsAsFactors = FALSE)[[1]]
  gdf <- read.delim(gf, header = FALSE, stringsAsFactors = FALSE)
  # 10x features/genes.tsv 第二列为 gene symbol，第一列 Ensembl
  sym <- as.character(gdf[[2]])
  na_idx <- is.na(sym) | sym == "" | sym == "NA"
  sym[na_idx] <- as.character(gdf[[1]])[na_idx]
  sym <- sub("\\..*$", "", sym)                       # 去 Ensembl 版本号
  dup <- duplicated(sym)
  if (any(dup)) sym[dup] <- paste0(sym[dup], ".dup", seq_len(sum(dup)))
  rownames(m) <- sym
  colnames(m) <- as.character(bc)
  if (exists("N_SUB", mode = "numeric") && ncol(m) > N_SUB)
    m <- m[, sample(seq_len(ncol(m)), N_SUB)]
  m
}

samples <- c("GSM4466780_F_19wks", "GSM4466783_P_6d")
N_SUB <- 4000
sce_list <- list()
for (s in samples) {
  m <- read_mtx(s)
  if (is.null(m)) next
  # 基础过滤：基因数 200-8000 的细胞（去双细胞与空液滴）
  lib <- Matrix::colSums(m); ngene <- Matrix::colSums(m > 0)
  keep <- ngene >= 200 & ngene <= 8000 & lib >= 1000
  if (sum(keep) < 50) { cat(sprintf("[scRNA] %s: 通过过滤细胞过少(%d)，跳过\n", s, sum(keep))); next }
  m <- m[, keep]
  if (ncol(m) > N_SUB) m <- m[, sample(seq_len(ncol(m)), N_SUB)]
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = m)
  )
  colData(sce)$sample <- s
  sce <- scater::logNormCounts(sce)
  cat(sprintf("[scRNA] %s: %d 细胞通过过滤\n", s, ncol(sce)))
  sce_list[[s]] <- sce
}

if (length(sce_list) == 0) stop("[scRNA] 无可用样本——终止")

# ---- 逐细胞签名评分（rank-based，AUCell 思路的轻量版）---------------------------
score_cells <- function(sce, sig) {
  lc <- logcounts(sce)                           # 稀疏 genes x cells
  sg <- unique(unlist(sig)); sg <- sg[sg %in% rownames(lc)]
  sub <- as.matrix(lc[sg, , drop = FALSE])       # 仅签名基因 x cells，极小
  rk <- apply(sub, 2, rank) / nrow(sub)          # 细胞内签名基因排名分位
  sapply(names(sig), function(nm) {
    g <- intersect(sig[[nm]], sg)
    if (length(g) < 5) return(rep(NA_real_, ncol(sub)))
    colMeans(rk[g, , drop = FALSE])
  })
}

all_scores <- list()
N_SUB <- 4000
for (s in names(sce_list)) {
  sce <- sce_list[[s]]
  if (ncol(sce) > N_SUB) sce <- sce[, sample(seq_len(ncol(sce)), N_SUB)]
  sc <- score_cells(sce, sig)
  all_scores[[s]] <- sc
  # PCA/UMAP
  lc <- logcounts(sce)
  gv <- Matrix::rowMeans(lc^2) - (Matrix::rowMeans(lc))^2
  hv <- names(sort(gv, decreasing = TRUE))[1:min(1000, length(gv))]
  sce <- scater::runPCA(sce, subset_row = hv, ncomponents = 30)
  sce <- scater::runUMAP(sce, dimred = "PCA")
  # 把评分挂到 colData
  scd <- as.data.frame(colData(sce))
  for (nm in colnames(sc)) scd[[nm]] <- sc[, nm]
  scd$DEFECT <- sc[, "CORTICAL"] - sc[, "MEDULLARY"]
  colData(sce) <- DataFrame(scd)
  assign(paste0("sce_", s), sce)
  saveRDS(sce, file.path(DIR_PROC, paste0("sce_", s, ".rds")))
}

# ---- Fig 2: UMAP 着色（CORTICAL / MEDULLARY / DEFECT）----------------------------
plot_umap_panel <- function(sce, s) {
  cd <- as.data.frame(colData(sce))
  red <- reducedDim(sce, "UMAP")
  df <- data.frame(UMAP1 = red[, 1], UMAP2 = red[, 2], cd)
  panels <- lapply(c("CORTICAL", "MEDULLARY", "DEFECT"), function(v) {
    ggplot(df, aes(.data[["UMAP1"]], .data[["UMAP2"]], color = .data[[v]])) +
      geom_point(size = 0.35, alpha = 0.8) +
      scale_color_viridis_c(name = "score") +
      labs(title = v, x = NULL, y = NULL) +
      theme_classic(base_size = 11) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
            legend.key.height = unit(0.4, "cm"))
  })
  gridExtra::grid.arrange(grobs = panels, ncol = 3,
                          top = grid::textGrob(paste0("GSE147520 — ", s),
                                               gp = grid::gpar(fontsize = 12, fontface = "bold")))
}
for (s in names(sce_list)) {
  sce <- get(paste0("sce_", s))
  p <- plot_umap_panel(sce, s)
  ggsave(file.path(DIR_FIG, paste0("fig2_scRNA_umap_", s, ".png")), p,
         width = 13, height = 4.6, dpi = 300)
}
cat("[scRNA] fig2 UMAP panels saved\n")

# ---- Fig 3: 签名相关性矩阵（证明非冗余，单样本）-------------------------------
sc_main <- all_scores[[samples[1]]]
keep <- complete.cases(sc_main)
cormat <- cor(sc_main[keep, , drop = FALSE])
sig_cols <- c("ETP", "DN", "DP", "SP_CD4", "SP_CD8", "TREG", "MTEC", "CTEC",
              "TRA", "MEDULLARY", "CORTICAL")
sig_cols <- intersect(sig_cols, colnames(cormat))
cm <- cormat[sig_cols, sig_cols]
long <- data.frame(
  x = factor(rep(sig_cols, each = length(sig_cols)), levels = sig_cols),
  y = factor(rep(sig_cols, times = length(sig_cols)), levels = rev(sig_cols)),
  r = as.vector(cm)
)
p3 <- ggplot(long, aes(x, y, fill = r)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.4, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-1, 1), name = "Pearson r") +
  labs(title = "Per-cell signature correlation (GSE147520 thymic cells)",
       subtitle = "Pre-locked thymopoiesis signatures are largely non-redundant",
       x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 12))
ggsave(file.path(DIR_FIG, "fig3_scRNA_sig_corr.png"), p3, width = 8, height = 7, dpi = FIG_DPI)
cat("[scRNA] fig3 signature correlation saved\n")
cat("[scRNA] 03b 完成\n")
