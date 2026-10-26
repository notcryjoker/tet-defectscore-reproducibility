#!/usr/bin/env Rscript
# =============================================================================
# 05_external.R —— GSE29695 外部验证（Illumina DASL, GPL8432, Badve 2012）
#   目标：在独立队列复现 ① WHO 谱系 DefectScore 梯度 ② Masaoka 分期关联
#   ③（探索）MG 关联与复发关联
#   实现要点:
#     - 本地解析 series matrix（表达表 + characteristics），不依赖 GEOquery 网络往返
#     - characteristics 行: "!Sample_characteristics_ch1" + tab 分隔的 "字段: 值"
#     - 排除 4 个细胞系样本（tissue == "Cell Line"）
#     - 探针->符号: 本地 GPL8432 注释（suppl annot），多探针按均值聚合
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
load_libs(c("data.table", "edgeR", "rstatix"))

gse <- "GSE29695"
mat_file <- file.path(DIR_RAW, paste0(gse, "_series_matrix.txt.gz"))
stopifnot(file.exists(mat_file))

# ---- 1. 本地解析 series matrix -------------------------------------------------
lines <- readLines(gzfile(mat_file))
i0 <- grep("^!series_matrix_table_begin", lines)
i1 <- grep("^!series_matrix_table_end", lines)
stopifnot(length(i0) == 1, length(i1) == 1, i1 > i0)
tbl_lines <- lines[(i0 + 1):(i1 - 1)]
hdr <- gsub('"', '', strsplit(tbl_lines[1], "\t")[[1]])
stopifnot(hdr[1] == "ID_REF")
gsm_ids <- hdr[-1]
# 数据行写入临时文件后用 fread 快速读入
tmp <- tempfile(fileext = ".txt")
writeLines(tbl_lines[-1], tmp)
expr_tbl <- fread(tmp, header = FALSE, sep = "\t", data.table = FALSE)
unlink(tmp)
colnames(expr_tbl) <- hdr
expr_mat <- as.matrix(expr_tbl[, -1, drop = FALSE])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- expr_tbl[, 1]
colnames(expr_mat) <- gsm_ids
cat(sprintf("[external] 表达矩阵 %d 探针 x %d 样本\n", nrow(expr_mat), ncol(expr_mat)))

# ---- 1b. characteristics 解析（行=字段, 列=样本 -> 转置为 样本 x 字段）----------
char_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
parse_char <- function(line) {
  body <- sub("^!Sample_characteristics_ch1\\t?", "", line)
  parts <- gsub('"', '', trimws(strsplit(body, "\t")[[1]]))
  named <- grepl(": ", parts)
  fname <- if (any(named)) sub(":.*$", "", parts[which(named)[1]]) else NA_character_
  vals <- ifelse(named, sub("^.*?: ", "", parts), parts)
  list(name = fname, vals = vals)
}
parsed <- lapply(char_lines, parse_char)
fv <- do.call(rbind, lapply(parsed, `[[`, "vals"))       # n_fields x n_samples
rownames(fv) <- vapply(parsed, function(x) x[["name"]], character(1))
fields <- as.data.frame(t(fv), stringsAsFactors = FALSE)  # 样本 x 字段
colnames(fields) <- rownames(fv)
rownames(fields) <- gsm_ids                                # 列序与表达矩阵一致
cat(sprintf("[external] characteristics 字段: %s\n", paste(colnames(fields), collapse = ", ")))

# ---- 1c. 排除细胞系 -------------------------------------------------------------
is_tumor <- fields$tissue == "Fresh Frozen Human Tumors"
clin <- fields[is_tumor, , drop = FALSE]
expr_mat <- expr_mat[, rownames(clin), drop = FALSE]
cat(sprintf("[external] 排除细胞系后保留 %d 例肿瘤样本\n", ncol(expr_mat)))

# ---- 1d. 临床变量归一化 ---------------------------------------------------------
# WHO: "Mixed AB"/"A/B"/"B1/B2" 为混合组织学；主分析剔除，敏感性取最严重成分
who_raw <- tolower(trimws(clin$type))
clin$WHO_mixed <- who_raw %in% c("mixed ab", "a/b", "b1/b2")
clin$WHO_worst <- ifelse(who_raw == "a", "a",
                  ifelse(who_raw %in% c("ab", "mixed ab", "a/b"), "ab",
                  ifelse(who_raw == "b1", "b1",
                  ifelse(who_raw == "b2", "b2",
                  ifelse(who_raw == "b3", "b3",
                  ifelse(who_raw == "b1/b2", "b2", NA_character_))))))
clin$WHO <- ifelse(clin$WHO_mixed, NA_character_, clin$WHO_worst)
clin$WHO <- factor(clin$WHO, levels = c("a", "ab", "b1", "b2", "b3"), ordered = TRUE)
# Masaoka: "I_II" / "III_IV" / "NA"（下划线分隔，非斜杠）
stage_tok <- tolower(trimws(clin$stage))
has34 <- grepl("(^|[^a-z])(iii|iv)($|[^a-z])", stage_tok)
has12 <- grepl("(^|[^a-z])(i|ii)($|[^a-z])", stage_tok)
clin$masaoka <- ifelse(has34, "III/IV", ifelse(has12, "I/II", NA_character_))
# 复发（大写 NO/YES）
rl <- tolower(trimws(clin$relapse))
clin$relapse <- ifelse(rl == "yes", 1, ifelse(rl == "no", 0, NA))
cat("[external] WHO 分布:\n"); print(table(clin$WHO, useNA = "ifany"))
cat("[external] Masaoka 分布:\n"); print(table(clin$masaoka, useNA = "ifany"))
cat("[external] 复发分布:\n"); print(table(clin$relapse, useNA = "ifany"))

# ---- 2. 探针->符号桥接（本地 Illumina manifest 注释，多探针取均值）--------------
annot_file <- file.path(DIR_RAW, "GPL8432_annot.txt.gz")
stopifnot(file.exists(annot_file))
al <- readLines(gzfile(annot_file))
p0 <- grep("^\\[Probes\\]", al)
stopifnot(length(p0) == 1)
nxt <- grep("^\\[", al); nxt <- nxt[nxt > p0]
p1 <- if (length(nxt)) min(nxt) - 1 else length(al)
a_hdr <- strsplit(al[p0 + 1], "\t")[[1]]
tmpa <- tempfile(fileext = ".txt")
writeLines(al[(p0 + 2):p1], tmpa)
a <- fread(tmpa, header = FALSE, sep = "\t", data.table = FALSE,
           quote = "", na.strings = NULL)
unlink(tmpa)
colnames(a) <- a_hdr
sym_col <- grep("Symbol", colnames(a), ignore.case = TRUE, value = TRUE)[1]
id_col  <- grep("Probe_Id", colnames(a), ignore.case = TRUE, value = TRUE)[1]
cat(sprintf("[external] manifest 注释 %d 行: ID=%s Symbol=%s\n", nrow(a), id_col, sym_col))
map_sym <- setNames(as.character(a[[sym_col]]), as.character(a[[id_col]]))
keep <- rownames(expr_mat) %in% names(map_sym)
expr_mat <- expr_mat[keep, , drop = FALSE]
sym_vec <- map_sym[rownames(expr_mat)]
ok <- !is.na(sym_vec) & nzchar(sym_vec) & sym_vec != "NA"
expr_mat <- expr_mat[ok, , drop = FALSE]
sym_vec <- sym_vec[ok]
# 多探针 -> 同一符号按均值聚合
agg <- rowsum(expr_mat, group = sym_vec, reorder = TRUE)
expr_sym <- agg / as.numeric(table(sym_vec)[rownames(agg)])
rownames(expr_sym) <- rownames(agg)
cat(sprintf("[external] 桥接后 %d 基因符号 x %d 样本\n", nrow(expr_sym), ncol(expr_sym)))

# ---- 3. DefectScore（与 TCGA 同法: cpm-log + z-score）--------------------------
lcpm <- cpm(expr_sym, log = TRUE, prior.count = 1)
sig <- readRDS(file.path(DIR_CODE, "thymo_signatures.rds"))
zscore_set <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 3) return(rep(NA, ncol(mat)))
  m <- colMeans(mat[g, , drop = FALSE], na.rm = TRUE)
  (m - mean(m, na.rm = TRUE)) / sd(m, na.rm = TRUE)
}
MED <- zscore_set(lcpm, sig$MEDULLARY)
COR <- zscore_set(lcpm, sig$CORTICAL)
defect <- COR - MED
actidx <- (COR + MED) / 2

# ---- 4. 复现分析 ----------------------------------------------------------------
df <- data.frame(sample = colnames(expr_sym), DEFECT = defect, ACTIDX = actidx,
                 stringsAsFactors = FALSE)
rownames(df) <- df$sample
df <- cbind(df, clin[, c("WHO", "WHO_worst", "masaoka", "relapse")])
df <- df[!is.na(df$DEFECT), ]

res <- list(n = nrow(df), who_levels = table(df$WHO, useNA = "ifany"))
# ① WHO 梯度（主分析: 纯组织学有序趋势）
df_who <- df[!is.na(df$WHO), ]
kw <- kruskal.test(DEFECT ~ WHO, data = df_who)
trend <- suppressWarnings(cor.test(as.numeric(df_who$WHO), df_who$DEFECT,
                                   method = "spearman"))
res$who_kruskal_p <- kw$p.value
res$who_trend_rho <- unname(trend$estimate)
res$who_trend_p <- trend$p.value
res$who_group_median <- tapply(df_who$DEFECT, df_who$WHO, median, na.rm = TRUE)
res$who_n_pure <- nrow(df_who)
cat(sprintf("[external] WHO 主分析: Kruskal p=%.3g | rho=%.2f p=%.3g (n=%d)\n",
            kw$p.value, unname(trend$estimate), trend$p.value, nrow(df_who)))
# 敏感性: worst-component（含混合型）
df_w <- df[!is.na(df$WHO_worst), ]
trend_w <- suppressWarnings(cor.test(as.numeric(factor(df_w$WHO_worst,
                       levels = c("a","ab","b1","b2","b3"), ordered = TRUE)),
                       df_w$DEFECT, method = "spearman"))
res$who_trend_worst_rho <- unname(trend_w$estimate)
res$who_trend_worst_p <- trend_w$p.value
cat(sprintf("[external] WHO 敏感性(worst-component): rho=%.2f p=%.3g (n=%d)\n",
            unname(trend_w$estimate), trend_w$p.value, nrow(df_w)))

# ② Masaoka 分期（I/II vs III/IV）
df_ms <- df[!is.na(df$masaoka), ]
if (nlevels(factor(df_ms$masaoka)) == 2) {
  w_ms <- wilcox.test(DEFECT ~ masaoka, data = df_ms)
  w_ms_a <- wilcox.test(ACTIDX ~ masaoka, data = df_ms)
  res$masaoka_wilcoxon_p <- w_ms$p.value
  res$masaoka_actidx_p <- w_ms_a$p.value
  res$masaoka_median <- tapply(df_ms$DEFECT, df_ms$masaoka, median)
  res$masaoka_actidx_median <- tapply(df_ms$ACTIDX, df_ms$masaoka, median)
  res$masaoka_n <- nrow(df_ms)
  cat(sprintf("[external] Masaoka I/II vs III/IV: DefectScore p=%.3g | ActIdx p=%.3g (n=%d)\n",
              w_ms$p.value, w_ms_a$p.value, nrow(df_ms)))
}

# ③ MG（GSE29695 无 MG 注释 -> 跳过）；复发关联（探索）
df_rl <- df[!is.na(df$relapse), ]
if (length(unique(df_rl$relapse)) == 2) {
  w_rl <- wilcox.test(DEFECT ~ factor(relapse), data = df_rl)
  res$relapse_wilcoxon_p <- w_rl$p.value
  res$relapse_n <- nrow(df_rl)
  cat(sprintf("[external] 复发 vs 无复发: Wilcoxon p=%.3g (n=%d)\n",
              w_rl$p.value, nrow(df_rl)))
}

saveRDS(res, file.path(DIR_RES, "external_gse29695.rds"))
saveRDS(df,  file.path(DIR_PROC, "external_defect.rds"))
cat("[external] 写出 external_gse29695.rds + external_defect.rds\n")
