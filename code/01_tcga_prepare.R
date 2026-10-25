#!/usr/bin/env Rscript
# =============================================================================
# 01_tcga_prepare.R —— TCGA-THYM 主队列制备
#   1) TCGAbiolinks 拉取 STAR counts（RNA-seq）
#   2) 合并 GDC 临床 biotab（MG / WHO / Masaoka）+ NTE（复发）表
#   3) 过滤非 TET 诊断，输出 expr_counts.rds + clinical_master.rds
# =============================================================================
source("/root/.codebuddy/artifact/2f8a7d30-9c21-4f0b-1d53-0f0283e72000/code/00_config.R")
load_libs(c("TCGAbiolinks", "data.table", "SummarizedExperiment"))

clin_dir <- file.path(DIR_RAW, "TCGA-THYM", "clinical")

# ---- 1. 表达矩阵 ------------------------------------------------------------
# 数据来源：GDC API 逐文件下载的 STAR-Counts TSV（见 code/dl_tcga_api.sh）
#   布局: data/raw/TCGA-THYM/gdc/<file_uuid>/<file_name>
#   策略: 直接读取 TSV 构建矩阵（绕开 TCGAbiolinks 分块下载的不稳定性）；
#         122 文件 / 120 病例 → 按病例聚合（同病例多文件取均值），保证患者级唯一
gdc_dir <- file.path(DIR_RAW, "TCGA-THYM", "gdc")
cache   <- file.path(DIR_PROC, "tcga_thym_se.rds")

if (file.exists(cache)) {
  cobj <- readRDS(cache)
  cat("[prepare] 读取缓存矩阵\n")
  if (is.list(cobj) && !is.null(cobj$counts)) {
    counts <- cobj$counts; map <- cobj$map
  } else {                                  # 兼容旧版 SE 缓存
    counts <- assay(cobj)
    rownames(counts) <- make.names(rownames(counts), unique = TRUE)
    map <- as.data.frame(colData(cobj))
  }
  rownames(counts) <- make.names(rownames(counts), unique = TRUE)
} else {
  file_list <- jsonlite::fromJSON("/tmp/thym_star_files.json", simplifyVector = FALSE)
  # 病例 -> 文件路径
  case_files <- split(file_list, factor(sapply(file_list, function(x) x$case)))
  tsvi <- gdc_file <- function(f) file.path(gdc_dir, f$file_id, f$file_name)
  # 读取单个 TSV：gene_id + unstranded counts + tpm（首行 # 注释需跳过；gene_id 剥版本号）
  read_one <- function(f) {
    d <- fread(tsvi(f), sep = "\t", header = TRUE, data.table = FALSE,
               comment.char = "#")
    d <- d[!is.na(d$unstranded), ]
    d$gene_id <- sub("\\.[0-9]+$", "", d$gene_id)   # ENSG... .15 -> ENSG...
    d <- d[grepl("^ENSG", d$gene_id), ]             # 去 N_unmapped 等汇总行
    d <- d[!duplicated(d$gene_id), ]
    data.frame(gene_id = d$gene_id,
               unstranded = suppressWarnings(as.numeric(d$unstranded)),
               tpm = suppressWarnings(as.numeric(d$tpm_unstranded)),
               stringsAsFactors = FALSE)
  }
  gene_order <- NULL
  per_case <- list()
  n_missing <- 0
  for (cs in names(case_files)) {
    fl <- case_files[[cs]]
    # 跳过未下载完成的文件（注意：fl 是 list，须逐元素判断，不能整体传入 file.path）
    keep <- vapply(fl, function(f)
      file.exists(file.path(gdc_dir, f$file_id, f$file_name)), logical(1))
    fl <- fl[keep]
    if (length(fl) == 0) { n_missing <- n_missing + 1; next }
    mats <- lapply(fl, read_one)
    if (is.null(gene_order)) gene_order <- mats[[1]]$gene_id
    # 同病例多文件取均值
    cnt <- do.call(cbind, lapply(mats, `[[`, "unstranded"))
    tpm <- do.call(cbind, lapply(mats, `[[`, "tpm"))
    per_case[[cs]] <- list(counts = rowMeans(cnt, na.rm = TRUE),
                           tpm = rowMeans(tpm, na.rm = TRUE))
  }
  if (n_missing > 0) cat(sprintf("[prepare] 警告: %d 个病例文件缺失被跳过（下载未完成）\n", n_missing))
  cat(sprintf("[prepare] 构建 %d 个病例的表达矩阵\n", length(per_case)))
  counts <- do.call(cbind, lapply(per_case, `[[`, "counts"))
  colnames(counts) <- names(per_case)
  rownames(counts) <- gene_order
  # 保存中间缓存（matrix + 病例注释），供后续直接复用
  map <- data.frame(submitter_id = names(per_case),
                    patient = substr(names(per_case), 1, 12),
                    stringsAsFactors = FALSE)
  saveRDS(list(counts = counts, map = map), cache)
  cat("[prepare] 本地构建完成并缓存（tcga_thym_se.rds）\n")
}
cat(sprintf("[prepare] 表达矩阵 %d 基因 x %d 样本\n", nrow(counts), ncol(counts)))

# ---- 2. 临床 biotab（MG / WHO / Masaoka） -----------------------------------
pat <- fread(file.path(clin_dir, "patient_clinical_clean.csv"), data.table = FALSE)
pat$patient <- substr(pat$bcr_patient_barcode, 1, 12)
# 统一 MG 字段
mg_raw <- tolower(as.character(pat$history_myasthenia_gravis))
pat$MG <- ifelse(grepl("yes", mg_raw), "MG+",
                 ifelse(grepl("no", mg_raw), "MG-", NA))

# NTE（新肿瘤事件 = 复发）表
nte_file <- file.path(clin_dir, "nationwidechildrens.org_clinical_nte_thym.txt")
if (file.exists(nte_file)) {
  # biotab 表头结构：第1行=真实列名(含 new_tumor_event_dx_days_to)，
  # 第2行=旧列名，第3行=CDE_ID 元数据；第2/3行作数据行读入后
  # barcode 为 "bcr_patient_barcode"/"CDE_ID:" 不会匹配真实病例，dx_days 为 NA，无害
  nte <- fread(nte_file, data.table = FALSE, fill = TRUE, skip = 0)
  nte$patient <- substr(nte$bcr_patient_barcode, 1, 12)
  # 复发事件时间（GDC biotab 实际列名：new_tumor_event_dx_days_to）
  dx_col <- grep("dx_days_to", colnames(nte), value = TRUE)[1]
  nte$dx_days <- suppressWarnings(as.numeric(nte[[dx_col]]))
  rfs <- aggregate(dx_days ~ patient, data = nte, FUN = function(x) min(x, na.rm = TRUE))
  rfs$has_recurrence <- 1
  colnames(rfs) <- c("patient", "rfs_time_days", "rfs_event")
} else {
  rfs <- data.frame(patient = character(), rfs_time_days = numeric(), rfs_event = numeric())
}

# ---- 3. 合并 -----------------------------------------------------------------
clin <- map
clin <- merge(clin, pat[, c("patient", "MG", "histologic_diagnosis",
                            "masaoka_stage", "vital_status",
                            "death_days_to", "last_contact_days_to",
                            "age_at_diagnosis", "gender")],
              by = "patient", all.x = TRUE)
clin <- merge(clin, rfs, by = "patient", all.x = TRUE)
clin$rfs_event[is.na(clin$rfs_event)] <- 0
clin$rfs_time_days[is.na(clin$rfs_time_days) & clin$vital_status == "Alive"] <-
  suppressWarnings(as.numeric(clin$last_contact_days_to))[is.na(clin$rfs_time_days) & clin$vital_status == "Alive"]

# ---- 4. 过滤非 TET 诊断 ------------------------------------------------------
# 保留 WHO A/AB/B1/B2/B3 与胸腺癌；排除基底样鳞癌/腺癌/囊腺癌等明确非 TET
keep_types <- c("type a", "type ab", "type b1", "type b2", "type b3",
                "thymic carcinoma", "carcinoma", "nos", "not reported")
clin$who_clean <- tolower(as.character(clin$histologic_diagnosis))
clin$is_tet <- grepl("a|ab|b1|b2|b3|thymic|carcinoma|nos", clin$who_clean) &
  !grepl("basaloid|adenocarcinoma|cystadenocarcinoma|squamous cell carcinoma", clin$who_clean)
cat(sprintf("[prepare] TET 保留 %d / %d 样本\n", sum(clin$is_tet), nrow(clin)))

# ---- 5. 输出 -----------------------------------------------------------------
expr_out <- counts[, as.character(clin$submitter_id[clin$is_tet])]
saveRDS(expr_out, file.path(DIR_PROC, "expr_counts.rds"))
saveRDS(clin[clin$is_tet, ], file.path(DIR_PROC, "clinical_master.rds"))
cat("[prepare] 已写出 expr_counts.rds + clinical_master.rds\n")
cat("MG 分布:", table(clin$MG[clin$is_tet]), "\n")
cat("WHO 分布:", table(clin$who_clean[clin$is_tet]), "\n")
