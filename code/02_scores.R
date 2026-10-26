#!/usr/bin/env Rscript
# =============================================================================
# 02_scores.R —— 胸腺发育缺陷评分（DefectScore）构造
#   输入: expr_counts.rds (STAR counts, Ensembl 行名) + thymo_signatures.rds
#   方法: ① GSVA-ssGSEA 样本级富集（主）  ② 基因集均值 z-score（交叉验证）
#   输出: scores_gsva.rds, scores_zscore.rds, defect_score.rds
#   候选基因集在关联检验前已预锁定（见 signature_definition.R）—— 禁止在此调整
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
load_libs(c("GSVA", "edgeR", "org.Hs.eg.db", "data.table"))

sig   <- readRDS(file.path(DIR_CODE, "thymo_signatures.rds"))
expr  <- readRDS(file.path(DIR_PROC, "expr_counts.rds"))   # Ensembl x sample
cat(sprintf("[score] 表达矩阵 %d x %d\n", nrow(expr), ncol(expr)))

# ---- Ensembl -> Symbol 映射 --------------------------------------------------
ens <- rownames(expr)
sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL")
sym <- as.character(sym)
names(sym) <- ens
expr_sym <- expr[!is.na(sym), ]
rownames(expr_sym) <- sym[!is.na(sym)]
dup <- duplicated(rownames(expr_sym))
expr_sym <- expr_sym[!dup, ]
cat(sprintf("[score] 映射为 symbol 矩阵 %d x %d\n", nrow(expr_sym), ncol(expr_sym)))

# ---- logCPM ------------------------------------------------------------------
lcpm <- cpm(expr_sym, log = TRUE, prior.count = 1)
lcpm[is.infinite(lcpm)] <- NA

# ---- 基因集过滤（仅保留在矩阵中存在的 symbol）-------------------------------
sig_f <- lapply(sig, function(g) intersect(g, rownames(lcpm)))
sizes <- lengths(sig_f)
cat("[score] 各基因集可用基因数:\n"); print(sizes)
small <- sizes < 5
if (any(small)) cat("[score] 警告: 以下基因集 <5 基因被剔除:", names(sig_f)[small], "\n")
sig_f <- sig_f[!small]

# ---- ① GSVA ssGSEA（主方法）------------------------------------------------
# GSVA 2.x API: 方法专用构造器 ssgseaParam（alpha=0.25 为 ssGSEA 默认指数）
gsva_param <- ssgseaParam(lcpm, sig_f, minSize = 5, maxSize = 500, alpha = 0.25)
gsva_res <- gsva(gsva_param)
colnames(gsva_res) <- colnames(lcpm)
saveRDS(gsva_res, file.path(DIR_PROC, "scores_gsva.rds"))
cat("[score] GSVA-ssGSEA 完成，写出 scores_gsva.rds\n")

# ---- ② 基因集均值 z-score（交叉验证方法）-----------------------------------
zscore_set <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 3) return(rep(NA, ncol(mat)))
  m <- colMeans(mat[g, , drop = FALSE], na.rm = TRUE)
  (m - mean(m, na.rm = TRUE)) / sd(m, na.rm = TRUE)
}
zs <- sapply(sig_f, function(g) zscore_set(lcpm, g))
colnames(zs) <- names(sig_f)
saveRDS(zs, file.path(DIR_PROC, "scores_zscore.rds"))

# ---- DefectScore = z(CORTICAL) - z(MEDULLARY) -------------------------------
# 注意：zs 是 样本 x 签名 矩阵，签名在列 —— 必须用列索引 zs[, "MEDULLARY"]
stopifnot(all(c("MEDULLARY", "CORTICAL") %in% colnames(zs)))
defect <- data.frame(
  sample = rownames(zs),
  MEDULLARY = zs[, "MEDULLARY"],
  CORTICAL  = zs[, "CORTICAL"],
  DEFECTSCORE = zs[, "CORTICAL"] - zs[, "MEDULLARY"]
)
saveRDS(defect, file.path(DIR_PROC, "defect_score.rds"))
cat("[score] DefectScore 构造完成（CORTICAL - MEDULLARY, 已 z-score）\n")
cat(sprintf("[score] DefectScore 范围: %.2f ~ %.2f\n",
            min(defect$DEFECTSCORE, na.rm = TRUE), max(defect$DEFECTSCORE, na.rm = TRUE)))
