#!/usr/bin/env Rscript
# =============================================================================
# signature_definition.R
# 胸腺上皮肿瘤（TET）课题 —— 胸腺发育程序细胞状态基因集（预锁定常量）
#
# 重要：本文件在"接触任何 TCGA 关联检验之前"即固定，作为候选基因预锁定
#      （candidate pre-locking）的证据，杜绝数据窥探（呼应 ESCC 稿方法学资产）。
#      基因集来源：人胸腺发育生物学共识 + GSE284479 scRNA 解析（M3.2 校准，
#      但本常量本身不因 TCGA 关联结果调整）。
# =============================================================================

# ---- 0. 胸腺发育各细胞状态基因集（HUGO 符号，人） --------------------------
# 参考：Human thymopoiesis staging (Cortés & Rocha 2023; Park et al. 2020;
#       Basse et al. Cancer Immunol Res 2026; GSE284479 作者注释)

# 早期胸腺祖细胞 ETP (CD34+ KIT+ CD1a-)
GENE_ETP <- c("CD34","KIT","MME","PROM1","LEF1","TCF7","CCR9","DNTT",
              "FLT3","HOXA9","MEF2C","LYL1","LMO2","GATA3","CD45RA","IL7R")

# 双阴性 DN (CD4- CD8-; CD44/CD25 动态)
GENE_DN <- c("CD44","IL2RA","CD7","CD5","BCL11B","TCF7","LCK","ZAP70",
             "CD3D","CD3E","TRDC","NOTCH1","PTCRA","RAG1","RAG2","DNTT","CD1A")

# 双阳性 DP (CD4+ CD8+; 皮质主体)
GENE_DP <- c("CD4","CD8A","CD8B","CD1A","TRAC","TRBC1","TRBC2","CD3D","CD3E","CD3G",
             "CORO1A","LCK","ZAP70","PRF1","LAT")

# 单阳性 SP-CD4 (髓质输出)
GENE_SP_CD4 <- c("CD4","CD3D","CD3E","CD3G","CCR7","CD40LG","IL7R","LEF1","TCF7","CD28","SAT1")

# 单阳性 SP-CD8
GENE_SP_CD8 <- c("CD8A","CD8B","CD3D","CD3E","CD3G","CCR7","GZMK","GZMA","CTSW","CD7","LEF1","CCL5")

# 调节性 T 细胞 Treg (FoxP3+; 中枢耐受核心)
GENE_TREG <- c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","TNFRSF4","IL10","TGFB1",
               "RTKN2","BATF","TIGIT","CCR8","LAYN","SAT1","LAIR2")

# 髓质胸腺上皮 mTEC (AIRE+; 阴性选择机器)
GENE_MTEC <- c("AIRE","KRT14","KRT5","CLDN3","CLDN4","CCR7","PSMB8","PSMB9","TAP1","TAP2",
               "CD40","CCL19","CCL21","PRSS16","LOR","IVL","KRT10","CD80","CD86","LTB")

# AIRE 依赖的组织限制性抗原 (TRA) 面板 —— 阴性选择完整性代理
GENE_TRA <- c("TPO","CHRNA1","S100B","AQP4","MBP","INS","TYR","SAG","OMP","PRL","GH1",
              "GAD2","CACNA1A","GRIN1","CACNA1D","SLC17A6","SLC6A2","TSHR","GAD1","RYR1")

# 皮质胸腺上皮 cTEC (DLL4+; 阳性选择机器)
GENE_CTEC <- c("DLL4","CORIN","LOXL2","PSMB11","CTSL","CLPTM1L","LTB","DSG1","KRT8",
               "KRT18","CAV1","CD205","CDR1","TBX1","PRSS16")

# ---- 1. 复合评分构造（bulk 可计算） -----------------------------------------
# Medullary（高 = 阴性选择完整）= mTEC + Treg + SP-CD4 + SP-CD8
# Cortical （高 = 成熟阻滞）    = ETP + DN + DP
# DefectScore = z(Cortical) - z(Medullary)  （连续，越高越缺陷）

GENE_MEDULLARY <- unique(c(GENE_MTEC, GENE_TRA, GENE_TREG, GENE_SP_CD4, GENE_SP_CD8))
GENE_CORTICAL  <- unique(c(GENE_ETP, GENE_DN, GENE_DP))

# 顶层导出列表（供 02_scores.R 使用）
THYMO_SIGNATURES <- list(
  ETP       = GENE_ETP,
  DN        = GENE_DN,
  DP        = GENE_DP,
  SP_CD4    = GENE_SP_CD4,
  SP_CD8    = GENE_SP_CD8,
  TREG      = GENE_TREG,
  MTEC      = GENE_MTEC,
  TRA       = GENE_TRA,
  CTEC      = GENE_CTEC,
  MEDULLARY = GENE_MEDULLARY,   # 复合：阴性选择完整度
  CORTICAL  = GENE_CORTICAL     # 复合：成熟阻滞
)

# 写入 results，供下游脚本 source
if (!exists("SIG_OUT", mode="numeric")) SIG_OUT <- NULL
saveRDS(THYMO_SIGNATURES,
        file = "file.path(Sys.getenv("TET_ROOT", "."), "code", "thymo_signatures.rds")")

cat("预锁定基因集就绪:\n")
for (nm in names(THYMO_SIGNATURES))
  cat(sprintf("  %-10s n=%d\n", nm, length(THYMO_SIGNATURES[[nm]])))
cat("MEDULLARY 复合 =", length(GENE_MEDULLARY), "| CORTICAL 复合 =", length(GENE_CORTICAL), "\n")
