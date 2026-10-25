#!/usr/bin/env Rscript
# =============================================================================
# 10_supplementary.R —— 投稿补充材料生成（基于管线真实输出）
#   Table S1: 预锁定胸腺生成细胞状态基因集（9 签名 + 2 复合评分）
#   Table S2: 样本级 DefectScore / ActIdx / ssGSEA 交叉验证评分 + 临床注释
#             （TCGA-THYM n=120；GSE29695 n=37）
#   Figure S1: 方法一致性散点图（ssGSEA vs 基因集均值 z-score）
# 输出目录: /workspace/supplementary/
# =============================================================================
source("/root/.codebuddy/artifact/2f8a7d30-9c21-4f0b-1d53-0f0283e72000/code/00_config.R")
load_libs(c("data.table", "ggplot2"))

OUT <- "/workspace/supplementary"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Table S1: 基因集 ----------------------------------------------------------
sig_env <- new.env()
sys.source(file.path(ROOT, "code/signature_definition.R"), envir = sig_env)
sig_keys <- c("ETP", "DN", "DP", "SP_CD4", "SP_CD8", "TREG", "MTEC", "TRA", "CTEC")
labels <- c(
  ETP     = "Early thymic progenitor (CD34+ KIT+ CD1a-)",
  DN      = "Double-negative (CD4- CD8-)",
  DP      = "Double-positive (CD4+ CD8+; cortical bulk)",
  SP_CD4  = "Single-positive CD4 (medullary output)",
  SP_CD8  = "Single-positive CD8 (medullary output)",
  TREG    = "Regulatory T cell (FoxP3+; central tolerance)",
  MTEC    = "Medullary TEC (AIRE+; negative-selection machinery)",
  TRA     = "AIRE-dependent tissue-restricted antigen panel",
  CTEC    = "Cortical TEC (DLL4+; positive-selection machinery)"
)
s1 <- rbindlist(lapply(sig_keys, function(k) {
  genes <- unique(get(paste0("GENE_", k), envir = sig_env))
  data.table(Signature = gsub("_", "-", k), `Cell state / programme` = labels[k],
             `Gene count` = length(genes), Genes = paste(genes, collapse = "; "))
}))
# 复合评分
cort <- unique(get("GENE_CORTICAL", envir = sig_env))
medu <- unique(get("GENE_MEDULLARY", envir = sig_env))
comp <- rbindlist(list(
  data.table(Signature = "CORTICAL (composite)", `Cell state / programme` = "Maturation arrest = ETP + DN + DP",
             `Gene count` = length(cort), Genes = paste(cort, collapse = "; ")),
  data.table(Signature = "MEDULLARY (composite)", `Cell state / programme` = "Negative-selection integrity = mTEC + TRA + Treg + SP-CD4 + SP-CD8",
             `Gene count` = length(medu), Genes = paste(medu, collapse = "; "))
))
s1 <- rbindlist(list(s1, comp))
fwrite(s1, file.path(OUT, "Table_S1_prelocked_gene_sets.csv"))
cat("[S1] rows:", nrow(s1), "| cortical:", length(cort), "medullary:", length(medu), "\n")

# ---- Table S2: 样本级数据 ------------------------------------------------------
tcga <- as.data.table(readRDS(file.path(DIR_PROC, "merged_defect_clinical.rds")))
gse  <- as.data.table(readRDS(file.path(DIR_PROC, "external_defect.rds")))
zs   <- readRDS(file.path(DIR_PROC, "scores_zscore.rds"))   # 120 x 11 (主方法)
gs   <- readRDS(file.path(DIR_PROC, "scores_gsva.rds"))     # 11 x 120 (交叉验证)

# ssGSEA 口径的 DefectScore（交叉验证）
ss <- data.table(
  sample          = colnames(gs),
  CORTICAL_ssGSEA = as.numeric(gs["CORTICAL", ]),
  MEDULLARY_ssGSEA= as.numeric(gs["MEDULLARY", ]),
  DefectScore_ssGSEA = as.numeric(gs["CORTICAL", ]) - as.numeric(gs["MEDULLARY", ])
)
tcga2 <- merge(tcga, ss, by = "sample")
stopifnot(nrow(tcga2) == 120)

tcga_out <- tcga2[, .(
  sample, cohort = "TCGA-THYM",
  WHO, WHO_worst, is_mixed, MG, masaoka_stage, Masaoka_grp,
  age_at_diagnosis, gender,
  CORTICAL_z = CORTICAL, MEDULLARY_z = MEDULLARY,
  DefectScore = DEFECTSCORE, ActIdx,
  CORTICAL_ssGSEA, MEDULLARY_ssGSEA, DefectScore_ssGSEA,
  rfs_time_days, rfs_event, vital_status, death_days_to, last_contact_days_to
)]
setnames(tcga_out, c("CORTICAL_z", "MEDULLARY_z"), c("CORTICAL (mean z)", "MEDULLARY (mean z)"))

gse_out <- gse[, .(
  sample, cohort = "GSE29695",
  WHO, WHO_worst, masaoka_stage = masaoka,
  DefectScore = DEFECT, ActIdx = ACTIDX, relapse
)]
fwrite(tcga_out, file.path(OUT, "Table_S2_TCGA-THYM_sample_scores.csv"))
fwrite(gse_out,  file.path(OUT, "Table_S2_GSE29695_sample_scores.csv"))
cat("[S2] TCGA rows:", nrow(tcga_out), "| GSE29695 rows:", nrow(gse_out), "\n")

# ---- Figure S1: 方法一致性散点 -------------------------------------------------
d <- data.table(
  sample = colnames(gs),
  Cortical_ssGSEA   = as.numeric(gs["CORTICAL", ]),
  Medullary_ssGSEA  = as.numeric(gs["MEDULLARY", ]),
  Cortical_z        = as.numeric(zs[, "CORTICAL"]),
  Medullary_z       = as.numeric(zs[, "MEDULLARY"])
)
r_c <- cor(d$Cortical_ssGSEA, d$Cortical_z, method = "pearson")
r_m <- cor(d$Medullary_ssGSEA, d$Medullary_z, method = "pearson")
cat(sprintf("[S1] concordance: Cortical r=%.2f, Medullary r=%.2f\n", r_c, r_m))

panels <- rbindlist(list(
  data.table(x = d$Cortical_ssGSEA,  y = d$Cortical_z,  panel = "Cortical composite"),
  data.table(x = d$Medullary_ssGSEA, y = d$Medullary_z, panel = "Medullary composite")
))
ann <- data.table(
  panel = c("Cortical composite", "Medullary composite"),
  rr    = c(r_c, r_m),
  xx    = c(min(panels$x) + 0.08 * diff(range(panels$x)), min(panels$x) + 0.08 * diff(range(panels$x))),
  yy    = c(max(panels$y) * 0.92, max(panels$y) * 0.92)
)
p <- ggplot(panels, aes(x = x, y = y)) +
  geom_point(color = "#2E86C1", alpha = 0.65, size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#B2561F", linewidth = 0.7, se = FALSE) +
  geom_text(data = ann, aes(x = xx, y = yy, label = sprintf("Pearson r = %.2f", rr)),
            hjust = 0, size = 4.2, fontface = "bold", color = "#1B4F72") +
  facet_wrap(~ panel, scales = "free") +
  labs(x = "ssGSEA enrichment score (GSVA)", y = "Gene-set mean z-score") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#EAF2F8"), strip.text = element_text(face = "bold"))
ggsave(file.path(OUT, "Figure_S1_method_concordance.png"), p,
       width = 9, height = 4.6, dpi = 300, bg = "white")
cat("[S1] figure saved\n")

# ---- 方法一致性摘要（写入补充说明） -------------------------------------------
info <- c(
  "Supplementary materials — generated from pipeline outputs (2026-10-23)",
  "Table S1: pre-locked thymopoiesis cell-state gene sets (nine signatures plus two composite scores).",
  "Table S2: sample-level DefectScore, ActIdx and clinical annotations for TCGA-THYM (n=120) and GSE29695 (n=37).",
  "  - Primary method: gene-set mean z-score on log2-CPM; DefectScore = z(CORTICAL) - z(MEDULLARY); ActIdx = (z(CORTICAL)+z(MEDULLARY))/2.",
  "  - Cross-check method: ssGSEA (GSVA, ssgseaParam, alpha=0.25); ssGSEA-based composite scores included for TCGA-THYM.",
  "Figure S1: method-concordance scatter (ssGSEA vs gene-set mean z-score) for the Cortical and Medullary composites (TCGA-THYM, n=120).",
  sprintf("Concordance: Cortical r=%.2f; Medullary r=%.2f.", r_c, r_m),
  "RFS columns: rfs_time_days (days to relapse or last contact), rfs_event (1=relapse).",
  "WHO_worst: worst-component assignment used in sensitivity analyses (mixed histologies).",
  "is_mixed: TRUE for composite WHO annotations excluded from primary analyses."
)
writeLines(info, file.path(OUT, "README_supplementary.txt"))
cat("[done] all supplementary materials written to", OUT, "\n")
