#!/usr/bin/env Rscript
# =============================================================================
# 04_survival.R —— 复发-free 生存（RFS）分析（唯一可行生存终点）
#   背景：TCGA-THYM 死亡事件仅 6-9 例 → OS 不可行，RFS（NTE 复发）唯一终点
#   克制原则：仅单变量 Cox + DefectScore 二分 KM；禁 ML 预后模型（事件数不足）
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
load_libs(c("survival", "survminer", "data.table", "ggplot2"))

m <- readRDS(file.path(DIR_PROC, "merged_defect_clinical.rds"))
cat(sprintf("[surv] 载入 %d 样本, 复发事件 %d\n", nrow(m), sum(m$rfs_event, na.rm=TRUE)))

# ---- 构造 RFS 时间（年）---------------------------------------------------
m$rfs_time_yrs <- suppressWarnings(as.numeric(m$rfs_time_days)) / 365.25
m$rfs_time_yrs[is.infinite(m$rfs_time_yrs)] <- NA
# 仅对非缺失值截断 >20 年的可疑随访（避免 pmin(NA,20,na.rm=TRUE) 把 NA 变成 20）
m$rfs_time_yrs[!is.na(m$rfs_time_yrs) & m$rfs_time_yrs > 20] <- 20

use <- m[!is.na(m$DEFECTSCORE) & !is.na(m$rfs_time_yrs) & m$rfs_time_yrs > 0, ]
cat(sprintf("[surv] 可分析 %d 样本 (事件 %d, 删失 %d)\n",
            nrow(use), sum(use$rfs_event), sum(use$rfs_time_yrs>0 & use$rfs_event==0)))

# ---- ① DefectScore 二分 KM -------------------------------------------------
use$defect_grp <- ifelse(use$DEFECTSCORE >= median(use$DEFECTSCORE, na.rm=TRUE),
                          "High defect", "Low defect")
fit <- survfit(Surv(rfs_time_yrs, rfs_event) ~ defect_grp, data = use)
logrank <- survdiff(Surv(rfs_time_yrs, rfs_event) ~ defect_grp, data = use)
pval <- 1 - pchisq(logrank$chisq, 1)

km_fig <- ggsurvplot(fit, data = use, pval = TRUE, risk.table = TRUE,
                     palette = c("#2166ac","#b2182b"),
                     xlab = "Time (years)", ylab = "Relapse-free survival",
                     title = "RFS by Thymopoiesis Defect Score")
ggsave(file.path(DIR_FIG, "fig8_rfs_km.png"), km_fig$plot,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI)

# ---- ② 单变量 Cox（DefectScore 主 + ActIdx 探索）-----------------------------
cox1 <- coxph(Surv(rfs_time_yrs, rfs_event) ~ DEFECTSCORE, data = use)
cox_sum <- summary(cox1)
res <- list(
  n = nrow(use), events = sum(use$rfs_event),
  km_logrank_p = pval,
  cox_hr = exp(cox_sum$coefficients[, "coef"]),
  cox_p = cox_sum$coefficients[, "Pr(>|z|)"],
  cox_ci = cox_sum$conf.int[, c("lower .95","upper .95")]  # conf.int 已是 HR 尺度，勿再 exp
)
cat(sprintf("[surv] KM log-rank p=%.3g | Cox HR=%.2f (95%%CI %.2f-%.2f) p=%.3g\n",
            pval, res$cox_hr, res$cox_ci[1], res$cox_ci[2], res$cox_p))

# 探索性: ActIdx（胸腺生成总活性）与 RFS
if ("ActIdx" %in% colnames(use)) {
  use$act_grp <- ifelse(use$ActIdx >= median(use$ActIdx, na.rm=TRUE),
                        "High activity", "Low activity")
  fit2 <- survfit(Surv(rfs_time_yrs, rfs_event) ~ act_grp, data = use)
  lr2 <- survdiff(Surv(rfs_time_yrs, rfs_event) ~ act_grp, data = use)
  p2 <- 1 - pchisq(lr2$chisq, 1)
  cox2 <- coxph(Surv(rfs_time_yrs, rfs_event) ~ ActIdx, data = use)
  s2 <- summary(cox2)
  km2 <- ggsurvplot(fit2, data = use, pval = TRUE, risk.table = TRUE,
                    palette = c("#b2182b","#2166ac"),
                    xlab = "Time (years)", ylab = "Relapse-free survival",
                    title = "RFS by thymopoietic activity index (exploratory)")
  ggsave(file.path(DIR_FIG, "fig8b_rfs_actidx.png"), km2$plot,
         width = FIG_W, height = FIG_H, dpi = FIG_DPI)
  res$actidx_km_p <- p2
  res$actidx_hr <- exp(s2$coefficients[, "coef"])
  res$actidx_p <- s2$coefficients[, "Pr(>|z|)"]
  res$actidx_ci <- s2$conf.int[, c("lower .95","upper .95")]  # 同上
  cat(sprintf("[surv] 探索 ActIdx: KM p=%.3g | HR=%.2f (95%%CI %.2f-%.2f) p=%.3g\n",
              p2, res$actidx_hr, res$actidx_ci[1], res$actidx_ci[2], res$actidx_p))
}
saveRDS(res, file.path(DIR_RES, "survival_rfs.rds"))
cat(sprintf("[surv] KM log-rank p=%.3g | Cox HR=%.2f (95%%CI %.2f-%.2f) p=%.3g\n",
            pval, res$cox_hr, res$cox_ci[1], res$cox_ci[2], res$cox_p))
cat("[surv] 写出 fig8_rfs_km.png + survival_rfs.rds\n")
