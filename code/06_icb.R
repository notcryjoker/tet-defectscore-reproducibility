#!/usr/bin/env Rscript
# =============================================================================
# 06_icb.R —— GSE181815 免疫治疗探索层（pembrolizumab, 胸腺癌 5R/5NR）
#   数据形态：GEO 可能仅 SRA 类型 → 先尝试取作者上传的表达矩阵+样本特征表
#   增强：从 characteristics 表识别应答标签；可行时输出 Fig 9；失败优雅降级
#   定位：探索性假设生成（非确证），讨论衔接 irAE/疗效
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
load_libs(c("GEOquery", "data.table", "edgeR", "org.Hs.eg.db", "ggplot2"))

gse <- "GSE181815"
res <- list(feasible = FALSE, note = "未获取到表达矩阵")
suppressWarnings({
  tabs <- tryCatch(getGSEDataTables(gse, destdir = DIR_RAW, do.log = FALSE),
                   error = function(e) NULL)
})

expr_tbl <- NULL
if (!is.null(tabs) && length(tabs) >= 1) {
  for (t in tabs) {
    if (ncol(t) > 15 && nrow(t) > 100) { expr_tbl <- t; break }
  }
}

if (!is.null(expr_tbl)) {
  idref <- as.character(expr_tbl[, 1])
  em <- as.matrix(expr_tbl[, -1, drop = FALSE]); mode(em) <- "numeric"
  rownames(em) <- idref
  # 平台注释桥接（ILMN 探针 → symbol）
  if (grepl("^ILMN|^A_|^[0-9]+$", idref[1])) {
    gpl <- tryCatch(getGPL("GPL17302", destdir = DIR_RAW, parseCharacteristics = FALSE), error = function(e) NULL)
    if (!is.null(gpl)) {
      ann <- as.data.frame(gpl@dataTable@table)
      sc <- grep("Symbol", colnames(ann), ignore.case = TRUE, value = TRUE)[1]
      idc <- grep("ID_REF", colnames(ann), ignore.case = TRUE, value = TRUE)[1]
      map_sym <- setNames(as.character(ann[[sc]]), as.character(ann[[idc]]))
      rownames(em) <- map_sym[rownames(em)]
    }
  }
  em <- em[!is.na(rownames(em)) & rownames(em) != "", ]
  em <- em[!duplicated(rownames(em)), ]
  lcpm <- cpm(em, log = TRUE, prior.count = 1)
  sig <- readRDS(file.path(DIR_CODE, "thymo_signatures.rds"))
  zs <- sapply(sig, function(g) {
    x <- intersect(g, rownames(lcpm)); if (length(x) < 3) return(NA)
    m <- colMeans(lcpm[x, , drop = FALSE], na.rm = TRUE)
    (m - mean(m, na.rm = TRUE)) / sd(m, na.rm = TRUE)
  })
  defect <- zs["CORTICAL", ] - zs["MEDULLARY", ]
  res$feasible <- TRUE
  res$defect <- defect
  res$n_samples <- length(defect)

  # ---- 应答标签识别：从样本名或特征表 ------------------------------------------
  labels <- rep(NA_character_, length(defect)); names(labels) <- names(defect)
  char_hit <- NULL
  if (!is.null(tabs)) {
    for (t in tabs) {
      cn <- colnames(t)
      if (any(grepl("response|outcome|characteristics", cn, ignore.case = TRUE)) &&
          nrow(t) <= 50 && nrow(t) >= 5) {
        # 找含应答信息的列
        rc <- grep("response|outcome|pfs|progression", cn, ignore.case = TRUE, value = TRUE)
        if (length(rc) >= 1) {
          key_col <- cn[1]
          vals <- as.character(t[[rc[1]]])
          if (any(grepl("responder|non-responder|respond|resist|CR|PR|SD|PD", vals, ignore.case = TRUE))) {
            char_hit <- data.frame(sample = as.character(t[[key_col]]), resp = vals,
                                   stringsAsFactors = FALSE)
            break
          }
        }
      }
    }
  }
  if (!is.null(char_hit)) {
    r <- ifelse(grepl("non-?respond|resist|^PD$|progressive", char_hit$resp, ignore.case = TRUE), "NR", "R")
    names(r) <- char_hit$sample
    m <- match(names(defect), names(r))
    labels <- r[m]
    res$labels <- labels
    res$n_labeled <- sum(!is.na(labels))
    cat(sprintf("[icb] 应答标签识别: R=%d NR=%d\n", sum(labels == "R", na.rm = TRUE), sum(labels == "NR", na.rm = TRUE)))
  } else {
    res$note <- paste(res$note, "；未识别到应答标签，仅报告 DefectScore 分布")
    cat("[icb] 未识别应答标签——降级为描述性\n")
  }

  # ---- Fig 9（仅在有标签时生成）-------------------------------------------------
  if (!is.null(char_hit) && sum(!is.na(labels)) >= 6) {
    df <- data.frame(
      sample = names(defect),
      DefectScore = as.numeric(defect),
      Response = factor(labels, levels = c("R", "NR"))
    )
    df <- df[!is.na(df$Response), ]
    p <- ggplot(df, aes(x = Response, y = DefectScore, fill = Response)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55) +
      geom_jitter(width = 0.12, size = 2.6, color = "grey20") +
      stat_summary(fun = mean, geom = "point", shape = 4, size = 3.5, color = "black") +
      scale_fill_manual(values = c("R" = "#1E8449", "NR" = "#943126")) +
      labs(x = "Pembrolizumab response", y = "DefectScore",
           title = "Exploratory: DefectScore and ICB response in thymic carcinoma (GSE181815)",
           subtitle = "Exploratory layer only — hypothesis-generating, not confirmatory") +
      theme_classic(base_size = 13) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 10, color = "grey40", face = "italic"))
    # Wilcoxon p 标注
    w <- tryCatch(wilcox.test(DefectScore ~ Response, data = df), error = function(e) NULL)
    if (!is.null(w)) p <- p + annotate("text", x = 1.5, y = max(df$DefectScore) * 1.02,
                                       label = sprintf("Wilcoxon p = %.3f", w$p.value), size = 4.5)
    ggsave(file.path(DIR_FIG, "fig9_icb_exploratory.png"), p,
           width = FIG_W * 0.85, height = FIG_H * 0.95, dpi = FIG_DPI)
    cat("[icb] Fig 9 已生成（fig9_icb_exploratory.png）\n")
  } else {
    cat("[icb] 标签不足，Fig 9 不生成——按计划降级为讨论衔接\n")
  }
  cat(sprintf("[icb] DefectScore 范围 %.2f ~ %.2f (n=%d)\n",
              min(defect, na.rm = TRUE), max(defect, na.rm = TRUE), length(defect)))
} else {
  res$note <- "GSEDataTables 无表达矩阵（SRA-only），需 SRA 重定量"
  cat("[icb] 仅 SRA 类型或无表——按预案降级为讨论衔接\n")
}

saveRDS(res, file.path(DIR_RES, "icb_gse181815.rds"))
cat("[icb] 结果已写 icb_gse181815.rds\n")
