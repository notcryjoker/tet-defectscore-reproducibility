#!/usr/bin/env Rscript
# =============================================================================
# 03_associations.R —— DefectScore 与临床/病理特征的关联（核心结果）
#   ① WHO 谱系梯度（有序，先验假设: C 最缺陷、A 最完整）
#   ② MG 分层（核心临床图）
#   ③ Masaoka 分期
#   ④ 敏感性：z-score 法 vs GSVA 法一致性
#   方法学：候选基因已预锁定；多重比较 FDR；克制建模
# =============================================================================
source("/root/.codebuddy/artifact/2f8a7d30-9c21-4f0b-1d53-0f0283e72000/code/00_config.R")
load_libs(c("data.table", "ggplot2", "ggpubr", "rstatix"))

defect <- readRDS(file.path(DIR_PROC, "defect_score.rds"))
clin   <- readRDS(file.path(DIR_PROC, "clinical_master.rds"))
gsva   <- readRDS(file.path(DIR_PROC, "scores_gsva.rds"))

# sample 匹配：defect$sample 是 submitter_id（含 -01A）；clin$submitter_id 同
m <- merge(defect, clin, by.x = "sample", by.y = "submitter_id", all.x = TRUE)
cat(sprintf("[assoc] 合并 %d 样本\n", nrow(m)))

# ---- WHO 有序化（鲁棒正则解析）------------------------------------------------
# clinical_master 实际值形如:
#   "thymoma; type a" / "thymoma; type a|thymoma; type ab" / "thymoma; type c"
# 策略: ① 从高级别向低级别填充 → 混合组织学取"最严重成分"（worst-component）
#       ② 主分析仅用纯组织学（剔除 | 混合型，n=114）；worst-component 作敏感性（n=120）
normalize_who <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  for (g in c("C", "B3", "B2", "B1", "AB", "A")) {
    pat <- switch(g,
      "C"  = "type c|thymic carcinoma",
      "B3" = "type b3",
      "B2" = "type b2",
      "B1" = "type b1",
      "AB" = "type ab",
      "A"  = "type a(\\b|$|\\|)")
    idx <- is.na(out) & grepl(pat, x)
    out[idx] <- g
  }
  out
}
m$WHO_worst <- normalize_who(m$histologic_diagnosis)
m$is_mixed   <- grepl("\\|", tolower(as.character(m$histologic_diagnosis)))
m$WHO        <- ifelse(m$is_mixed, NA_character_, m$WHO_worst)   # 主分析: 纯组织学
m$WHO_ord    <- factor(m$WHO_worst, levels = c("A","AB","B1","B2","B3","C"),
                       ordered = TRUE)
m <- m[!is.na(m$DEFECTSCORE), ]
cat(sprintf("[assoc] WHO 解析: 纯组织学 n=%d（A/AB/B1/B2/B3/C = %s）；混合型 n=%d\n",
            sum(!is.na(m$WHO)),
            paste(table(m$WHO), collapse="/"), sum(m$is_mixed)))

res <- list()

# ① WHO 梯度（主分析: 纯组织学有序趋势）
m_pure <- m[!is.na(m$WHO), ]
kw <- kruskal.test(DEFECTSCORE ~ WHO_ord, data = m_pure)
# 趋势检验（Spearman 等级相关，有序 WHO）
trend <- cor.test(as.numeric(m_pure$WHO_ord), m_pure$DEFECTSCORE, method = "spearman")
res$WHO_kruskal_p <- kw$p.value
res$WHO_trend_rho <- trend$estimate; res$WHO_trend_p <- trend$p.value
res$WHO_n_pure <- nrow(m_pure)
res$WHO_group_median <- tapply(m_pure$DEFECTSCORE, m_pure$WHO, median, na.rm = TRUE)
cat(sprintf("[assoc] WHO Kruskal p=%.3g | Spearman rho=%.2f p=%.3g (纯组织学 n=%d)\n",
            kw$p.value, trend$estimate, trend$p.value, nrow(m_pure)))
# 敏感性: worst-component 赋值（含混合型, n=120）
trend_w <- suppressWarnings(cor.test(as.numeric(m$WHO_ord), m$DEFECTSCORE,
                                     method = "spearman"))
res$WHO_trend_worst_rho <- unname(trend_w$estimate)
res$WHO_trend_worst_p   <- trend_w$p.value
cat(sprintf("[assoc] WHO 敏感性(worst-component): rho=%.2f p=%.3g (n=%d)\n",
            unname(trend_w$estimate), trend_w$p.value, nrow(m)))

# ② MG 分层（核心）
mg_use <- m[!is.na(m$MG), ]
t_mg <- wilcox.test(DEFECTSCORE ~ MG, data = mg_use)
res$MG_wilcox_p <- t_mg$p.value
res$MG_median <- tapply(mg_use$DEFECTSCORE, mg_use$MG, median, na.rm = TRUE)
cat(sprintf("[assoc] MG Wilcoxon p=%.3g | MG+ median=%.2f MG- median=%.2f\n",
            t_mg$p.value, res$MG_median["MG+"], res$MG_median["MG-"]))

# ③ Masaoka 分期（去字母法: "iia"→"ii"、"ivb"→"iv"；先判 III/IV 再判 I/II）
ms_raw <- gsub("[^iv]", "", tolower(trimws(as.character(m$masaoka_stage))))
m$Masaoka_grp <- ifelse(grepl("^(iii|iv)", ms_raw), "III-IV",
                        ifelse(grepl("^(i|ii)", ms_raw), "I-II", NA))
mk <- m[!is.na(m$Masaoka_grp), ]
t_mk <- wilcox.test(DEFECTSCORE ~ Masaoka_grp, data = mk)
res$Masaoka_wilcox_p <- t_mk$p.value
res$Masaoka_median <- tapply(mk$DEFECTSCORE, mk$Masaoka_grp, median, na.rm = TRUE)
cat(sprintf("[assoc] Masaoka Wilcoxon p=%.3g | I-II n=%d III-IV n=%d\n",
            t_mk$p.value, sum(mk$Masaoka_grp=="I-II"), sum(mk$Masaoka_grp=="III-IV")))

# ④ 敏感性：GSVA-MEDULLARY/CORTICAL 与 z-score 法一致性
gsva_long <- as.data.frame(t(gsva))
gsva_long$sample <- rownames(gsva_long)
mg <- merge(m, gsva_long, by = "sample")   # 同名列自动加 .x(zscore) / .y(gsva)
g_cols <- intersect(c("MEDULLARY","CORTICAL"), colnames(gsva_long))
cor_check <- sapply(g_cols, function(cn) {
  x <- mg[[paste0(cn, ".x")]]   # z-score 法
  y <- mg[[paste0(cn, ".y")]]   # GSVA 法
  if (sum(!is.na(x) & !is.na(y)) < 5) return(NA)
  cor(x, y, use = "pairwise.complete")
})
res$gsva_vs_zscore_corr <- cor_check

# ⑤ MG 多因素 logistic 回归（Methods §2.4 承诺；≤3 协变量防过拟合）
mg_fit <- m[!is.na(m$MG), ]
mg_fit$MG_bin <- ifelse(mg_fit$MG == "MG+", 1, 0)
mg_fit$age <- suppressWarnings(as.numeric(mg_fit$age_at_diagnosis))
mg_fit$sex <- ifelse(tolower(as.character(mg_fit$gender)) == "male", 1, 0)
# 协变量：DefectScore + WHO_ord(worst-component, n=120 完整) + age；缺失行剔除
fit_dat <- mg_fit[!is.na(mg_fit$MG_bin) & !is.na(mg_fit$DEFECTSCORE) & !is.na(mg_fit$age), ]
if (nrow(fit_dat) >= 30 && sum(fit_dat$MG_bin) >= 10) {
  lr <- glm(MG_bin ~ DEFECTSCORE + as.numeric(WHO_ord) + age,
            data = fit_dat, family = binomial)
  s <- summary(lr)$coefficients
  res$MG_logistic <- data.frame(term = rownames(s), beta = s[, 1], se = s[, 2],
                                p = s[, 4], stringsAsFactors = FALSE)
  cat("[assoc] MG logistic 回归:\n"); print(res$MG_logistic)
} else {
  res$MG_logistic <- NULL
  cat("[assoc] MG logistic 样本不足，跳过\n")
}

# ⑥ 探索性（事后分析，稿件中须明确标注为 post-hoc）:
#    胸腺生成总活性 ActIdx = (z(CORTICAL)+z(MEDULLARY))/2
#    动机: DefectScore 仅度量皮质/髓质"平衡"，不度量绝对胸腺生成活性；
#          A 型与 C 型均呈低活性（耗竭），但平衡方向相反
m$ActIdx <- (m$CORTICAL + m$MEDULLARY) / 2
kw_a <- kruskal.test(ActIdx ~ WHO_ord, data = m[!is.na(m$WHO), ])
res$ActIdx_WHO_kruskal_p <- kw_a$p.value
res$ActIdx_WHO_median <- tapply(m$ActIdx[!is.na(m$WHO)], m$WHO[!is.na(m$WHO)],
                                median, na.rm = TRUE)
t_amg <- wilcox.test(ActIdx ~ MG, data = m[!is.na(m$MG), ])
res$ActIdx_MG_wilcox_p <- t_amg$p.value
res$ActIdx_MG_median <- tapply(m$ActIdx, m$MG, median, na.rm = TRUE)
t_amk <- wilcox.test(ActIdx ~ Masaoka_grp, data = m[!is.na(m$Masaoka_grp), ])
res$ActIdx_Masaoka_wilcox_p <- t_amk$p.value
res$ActIdx_Masaoka_median <- tapply(m$ActIdx, m$Masaoka_grp, median, na.rm = TRUE)
cat(sprintf("[assoc] 探索 ActIdx: WHO KW p=%.3g | MG p=%.3g | Masaoka p=%.3g\n",
            kw_a$p.value, t_amg$p.value, t_amk$p.value))

# ⑦ 主检验 FDR（BH，跨 WHO/MG/Masaoka 三个主假设；探索性检验不入 FDR 家族）
pvals <- c(res$WHO_trend_p, res$MG_wilcox_p, res$Masaoka_wilcox_p)
names(pvals) <- c("WHO_trend", "MG", "Masaoka")
res$fdr <- p.adjust(pvals, method = "BH")
cat(sprintf("[assoc] BH-FDR: WHO=%.3g MG=%.3g Masaoka=%.3g\n",
            res$fdr["WHO_trend"], res$fdr["MG"], res$fdr["Masaoka"]))

saveRDS(res, file.path(DIR_RES, "association_results.rds"))
saveRDS(m,  file.path(DIR_PROC, "merged_defect_clinical.rds"))
cat("[assoc] 写出 association_results.rds + merged_defect_clinical.rds\n")
