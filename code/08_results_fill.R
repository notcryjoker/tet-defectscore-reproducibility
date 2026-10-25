#!/usr/bin/env Rscript
# =============================================================================
# 08_results_fill.R —— 读取各步 .rds 结果，编译 Results 文本 + 更新稿件
#   输出: results/results_summary.md（供人工填入 manuscript_draft.md §3）
# =============================================================================
source("/root/.codebuddy/artifact/2f8a7d30-9c21-4f0b-1d53-0f0283e72000/code/00_config.R")

out <- c()
cat_line <- function(...) { line <- paste(...); out <<- c(out, line); cat(line, "\n") }

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# ---- 队列与关联结果 ----
if (file.exists(file.path(DIR_RES, "association_results.rds"))) {
  r <- readRDS(file.path(DIR_RES, "association_results.rds"))
  cat_line("## 3.1 队列与 DefectScore 分布")
  cat_line(sprintf("- TCGA-THYM 可分析 n=120；纯组织学 n=%s（A/AB/B1/B2/B3/C）", r$WHO_n_pure))
  cat_line(sprintf("- DefectScore 取值连续；WHO 各组中位: %s",
                   paste(names(r$WHO_group_median), fmt(r$WHO_group_median), sep="=", collapse=", ")))
  cat_line("")
  cat_line("## 3.2 WHO 谱系梯度（主假设）")
  cat_line(sprintf("- Kruskal-Wallis p = %.3g（纯组织学 n=%s）", r$WHO_kruskal_p, r$WHO_n_pure))
  cat_line(sprintf("- Spearman 趋势 rho = %.2f (p = %.3g)：方向与预锁定单调假设相反",
                   r$WHO_trend_rho, r$WHO_trend_p))
  cat_line("  - 形态: 非单调（倒 U，B1 峰值）；预锁定 A→C 单调递增假设未获支持")
  cat_line(sprintf("- 敏感性（worst-component 含混合型 n=120）: rho = %.2f (p = %.3g)",
                   r$WHO_trend_worst_rho, r$WHO_trend_worst_p))
  cat_line("")
  cat_line("## 3.3 MG 分层（核心临床图）")
  cat_line(sprintf("- Wilcoxon p = %.3g（MG+ vs MG−）", r$MG_wilcox_p))
  cat_line(sprintf("- MG+ 中位 DefectScore = %s；MG− = %s",
                   fmt(r$MG_median["MG+"]), fmt(r$MG_median["MG-"])))
  if (!is.null(r$MG_logistic))
    cat_line(sprintf("- 多因素 logistic（+WHO+年龄）DefectScore OR = %.2f, p = %.3g",
                     exp(r$MG_logistic$beta[2]), r$MG_logistic$p[2]))
  cat_line("")
  cat_line("## 3.4 Masaoka 分期与探索性活性指数")
  cat_line(sprintf("- DefectScore Masaoka I-II vs III-IV: Wilcoxon p = %.3g", r$Masaoka_wilcox_p))
  cat_line(sprintf("- 探索 ActIdx: WHO Kruskal p = %.3g | MG p = %.3g | Masaoka p = %.3g",
                   r$ActIdx_WHO_kruskal_p, r$ActIdx_MG_wilcox_p, r$ActIdx_Masaoka_wilcox_p))
  cat_line(sprintf("- ActIdx Masaoka 中位: I-II = %s；III-IV = %s",
                   fmt(r$ActIdx_Masaoka_median["I-II"]), fmt(r$ActIdx_Masaoka_median["III-IV"])))
  cat_line("")
  cat_line("## 3.5 多重校正（BH-FDR，三个主假设）")
  cat_line(sprintf("- WHO 趋势 FDR = %.3g；MG FDR = %.3g；Masaoka FDR = %.3g",
                   r$fdr["WHO_trend"], r$fdr["MG"], r$fdr["Masaoka"]))
  cat_line("- 注: ActIdx 为事后探索性指标，不纳入主假设 FDR 家族")
  if (!is.null(r$gsva_vs_zscore_corr))
    cat_line(sprintf("- 稳健性: GSVA 法与 z-score 法相关性 %s",
                     paste(names(r$gsva_vs_zscore_corr), fmt(r$gsva_vs_zscore_corr), sep="=", collapse=", ")))
}

# ---- 生存 ----
if (file.exists(file.path(DIR_RES, "survival_rfs.rds"))) {
  s <- readRDS(file.path(DIR_RES, "survival_rfs.rds"))
  cat_line("")
  cat_line("## 3.6 无复发生存（RFS）")
  cat_line(sprintf("- 可分析 n=%d，复发事件 %d（事件稀疏，仅单变量 Cox）", s$n, s$events))
  cat_line(sprintf("- KM log-rank p = %.3g", s$km_logrank_p))
  cat_line(sprintf("- 单变量 Cox HR = %.2f (95%%CI %.2f–%.2f), p = %.3g",
                   s$cox_hr, s$cox_ci[1], s$cox_ci[2], s$cox_p))
  if (!is.null(s$actidx_hr))
    cat_line(sprintf("- 探索 ActIdx: HR = %.2f (95%%CI %.2f–%.2f), p = %.3g",
                     s$actidx_hr, s$actidx_ci[1], s$actidx_ci[2], s$actidx_p))
}

# ---- 外部验证 ----
if (file.exists(file.path(DIR_RES, "external_gse29695.rds"))) {
  e <- readRDS(file.path(DIR_RES, "external_gse29695.rds"))
  cat_line("")
  cat_line("## 3.7 外部验证 GSE29695（独立队列，Illumina DASL）")
  cat_line(sprintf("- 肿瘤样本 n=%d（排除 4 个细胞系）；纯组织学 n=%s", e$n, e$who_n_pure))
  cat_line(sprintf("- WHO Kruskal p = %.3g；趋势 rho = %.2f (p = %.3g)",
                   e$who_kruskal_p, e$who_trend_rho, e$who_trend_p))
  cat_line("  - 方向与 TCGA-THYM 一致（均为负向），构成独立复现")
  cat_line(sprintf("- 敏感性（worst-component n=37）: rho = %.2f (p = %.3g)",
                   e$who_trend_worst_rho, e$who_trend_worst_p))
  if (!is.null(e$who_group_median))
    cat_line(sprintf("- WHO 各组中位 DefectScore: %s",
                     paste(names(e$who_group_median), fmt(e$who_group_median), sep="=", collapse=", ")))
  if (!is.null(e$masaoka_wilcoxon_p))
    cat_line(sprintf("- Masaoka I/II vs III/IV: DefectScore p = %.3g；ActIdx p = %.3g (n=%d)",
                     e$masaoka_wilcoxon_p, e$masaoka_actidx_p, e$masaoka_n))
  if (!is.null(e$relapse_wilcoxon_p))
    cat_line(sprintf("- 复发 vs 无复发（探索）: p = %.3g (n=%d)",
                     e$relapse_wilcoxon_p, e$relapse_n))
}

# ---- ICB ----
if (file.exists(file.path(DIR_RES, "icb_gse181815.rds"))) {
  i <- readRDS(file.path(DIR_RES, "icb_gse181815.rds"))
  cat_line("")
  cat_line("## 3.8 ICB 探索层")
  cat_line(sprintf("- feasible = %s；note = %s", i$feasible, i$note))
}

writeLines(out, file.path(DIR_RES, "results_summary.md"))
cat("[08] 结果摘要已写出 results/results_summary.md\n")
