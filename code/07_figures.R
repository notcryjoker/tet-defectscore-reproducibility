#!/usr/bin/env Rscript
# =============================================================================
# 07_figures.R —— 主图生成（Fig 4/5/6/7/9/10, 300 dpi 出版级）
#   依赖: merged_defect_clinical.rds, external_defect.rds, association_results.rds
#   说明: WHO 图仅用纯组织学（混合型已置 NA，绘图时自动剔除并在正文注明 n=6）
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
load_libs(c("ggplot2","ggpubr","data.table","patchwork","RColorBrewer"))

m <- readRDS(file.path(DIR_PROC, "merged_defect_clinical.rds"))
res <- readRDS(file.path(DIR_RES, "association_results.rds"))
if (!"ActIdx" %in% colnames(m) && file.exists(file.path(DIR_PROC, "defect_score.rds"))) {
  d <- readRDS(file.path(DIR_PROC, "defect_score.rds"))
  m <- merge(m, d[, c("sample","ActIdx")], by = "sample")
}
m_pure <- m[!is.na(m$WHO), ]
m_pure$WHO <- factor(m_pure$WHO, levels = c("A","AB","B1","B2","B3","C"), ordered = TRUE)

who_pal  <- c("A"="#2166ac","AB"="#67a9cf","B1"="#f7f7f7","B2"="#ef8a62",
              "B3"="#d6604d","C"="#b2182b")
tw <- theme_pubclean(base_size = 11)
title_small <- theme(plot.title = element_text(size = 10, hjust = 0))

# ---- Fig 4: WHO 谱系（DefectScore 主 + ActIdx 探索）---------------------------
p4a <- ggboxplot(m_pure, x = "WHO", y = "DEFECTSCORE", fill = "WHO",
                 palette = who_pal, add = "jitter", width = 0.7,
                 title = "A  Thymopoiesis Defect Score by WHO histotype",
                 xlab = "WHO type", ylab = "DefectScore (Cortical − Medullary)") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y.npc = 0.95) +
  theme_pubclean(base_size = 11)
p4b <- ggboxplot(m_pure, x = "WHO", y = "ActIdx", fill = "WHO",
                 palette = who_pal, add = "jitter", width = 0.7,
                 title = "B  Thymopoietic activity index (exploratory)",
                 xlab = "WHO type", ylab = "Activity index (mean z)") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y.npc = 0.95) +
  theme_pubclean(base_size = 11)
p4 <- p4a + p4b + plot_layout(ncol = 2, guides = "collect") & title_small
ggsave(file.path(DIR_FIG, "fig4_who_gradient.png"), p4,
       width = FIG_W * 1.5, height = FIG_H * 0.85, dpi = FIG_DPI)

# ---- Fig 5: MG 分层（核心临床图）----------------------------------------------
mg <- m[!is.na(m$MG), ]
mg$MG <- factor(mg$MG, levels = c("MG-","MG+"))
p5a <- ggboxplot(mg, x = "MG", y = "DEFECTSCORE", fill = "MG",
                 palette = c("MG-"="#2166ac","MG+"="#b2182b"), add = "jitter", width = 0.6,
                 title = "A  DefectScore by myasthenia gravis status",
                 xlab = "", ylab = "DefectScore") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.95) +
  theme_pubclean(base_size = 11)
p5b <- ggboxplot(mg, x = "MG", y = "ActIdx", fill = "MG",
                 palette = c("MG-"="#2166ac","MG+"="#b2182b"), add = "jitter", width = 0.6,
                 title = "B  Activity index by MG status",
                 xlab = "", ylab = "Activity index") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.95) +
  theme_pubclean(base_size = 11)
p5 <- p5a + p5b + plot_layout(ncol = 2, guides = "collect") & title_small
ggsave(file.path(DIR_FIG, "fig5_mg.png"), p5,
       width = FIG_W * 1.3, height = FIG_H * 0.85, dpi = FIG_DPI)

# ---- Fig 6: Masaoka 分期 + 散点景观 --------------------------------------------
mk <- m[!is.na(m$Masaoka_grp), ]
mk$Masaoka_grp <- factor(mk$Masaoka_grp, levels = c("I-II","III-IV"))
p6a <- ggboxplot(mk, x = "Masaoka_grp", y = "DEFECTSCORE", fill = "Masaoka_grp",
                 palette = c("I-II"="#2166ac","III-IV"="#b2182b"), add = "jitter", width = 0.6,
                 title = "A  DefectScore by Masaoka stage", xlab = "", ylab = "DefectScore") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.95) +
  theme_pubclean(base_size = 11)
p6b <- ggboxplot(mk, x = "Masaoka_grp", y = "ActIdx", fill = "Masaoka_grp",
                 palette = c("I-II"="#2166ac","III-IV"="#b2182b"), add = "jitter", width = 0.6,
                 title = "B  Activity index by Masaoka stage", xlab = "", ylab = "Activity index") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.95) +
  theme_pubclean(base_size = 11)
p6 <- p6a + p6b + plot_layout(ncol = 2, guides = "collect") & title_small
ggsave(file.path(DIR_FIG, "fig6_masaoka.png"), p6,
       width = FIG_W * 1.3, height = FIG_H * 0.85, dpi = FIG_DPI)

# ---- Fig 9: 皮质-髓质平衡散点（按 WHO 着色, MG 标形）---------------------------
scatter_pal <- c("A"="#2166ac","AB"="#67a9cf","B1"="#4d4d4d","B2"="#ef8a62",
                 "B3"="#d6604d","C"="#b2182b")
p9 <- ggplot(m, aes(x = MEDULLARY, y = CORTICAL, colour = WHO_worst, shape = MG)) +
  geom_point(alpha = 0.8, size = 2.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  scale_colour_manual(values = scatter_pal, name = "WHO") +
  scale_shape_manual(values = c("MG-" = 16, "MG+" = 17), na.value = 1, name = "MG") +
  labs(title = "Cortical vs medullary programme balance in TET",
       x = "Medullary programme (z)", y = "Cortical programme (z)") +
  theme_pubclean(base_size = 11)
ggsave(file.path(DIR_FIG, "fig9_balance_scatter.png"), p9,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI)

# ---- Fig 7: 外部验证（GSE29695）------------------------------------------------
if (file.exists(file.path(DIR_PROC, "external_defect.rds"))) {
  ext <- readRDS(file.path(DIR_PROC, "external_defect.rds"))
  epal <- c("a"="#2166ac","ab"="#67a9cf","b1"="#f7f7f7","b2"="#ef8a62","b3"="#b2182b")
  p7a <- ggboxplot(ext[!is.na(ext$WHO), ], x = "WHO", y = "DEFECT", fill = "WHO",
                   palette = epal, add = "jitter", width = 0.7,
                   title = "A  GSE29695: DefectScore by WHO",
                   xlab = "WHO type", ylab = "DefectScore") +
    stat_compare_means(method = "kruskal.test", label = "p.format", label.y.npc = 0.95) +
    theme_pubclean(base_size = 11)
  em <- ext[!is.na(ext$masaoka), ]
  p7b <- ggboxplot(em, x = "masaoka", y = "DEFECT", fill = "masaoka",
                   palette = c("I/II"="#2166ac","III/IV"="#b2182b"), add = "jitter", width = 0.6,
                   title = "B  GSE29695: DefectScore by Masaoka",
                   xlab = "", ylab = "DefectScore") +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.95) +
    theme_pubclean(base_size = 11)
  p7 <- p7a + p7b + plot_layout(ncol = 2, guides = "collect") & title_small
  ggsave(file.path(DIR_FIG, "fig7_external.png"), p7,
         width = FIG_W * 1.3, height = FIG_H * 0.85, dpi = FIG_DPI)
  cat("[fig] Fig7 外部验证图已生成\n")
}

# ---- Fig 10: 机制示意图 --------------------------------------------------------
mech <- data.frame(
  stage = factor(c("Cortical programme\n(ETP / DN / DP)",
                   "Medullary programme\n(mTEC / TRA / Treg / SP)",
                   "Clinical correlates\n(histotype / stage / MG)"),
                 levels = c("Cortical programme\n(ETP / DN / DP)",
                            "Medullary programme\n(mTEC / TRA / Treg / SP)",
                            "Clinical correlates\n(histotype / stage / MG)")),
  x = 1:3, y = c(3, 2, 1))
p10 <- ggplot(mech, aes(x, y, label = stage)) +
  geom_label(size = 4.5, fill = "#f7f7f7", label.size = 0.4) +
  geom_segment(aes(x = 1, xend = 2, y = 3, yend = 2),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 1),
               arrow = arrow(length = unit(0.25, "cm"))) +
  xlim(0.6, 3.4) + ylim(0.6, 3.4) + theme_void() +
  labs(title = "Thymopoietic programme balance in thymic epithelial tumours")
ggsave(file.path(DIR_FIG, "fig10_mechanism.png"), p10,
       width = FIG_W, height = FIG_H * 0.5, dpi = FIG_DPI)

cat("[fig] Fig 4/5/6/7/9/10 生成完毕\n")
cat("输出目录:", DIR_FIG, "\n")
