## Comprehensive install: CRAN (PPM binary) + Bioconductor
options(repos = c(PPM = "https://packagemanager.posit.co/cran/__linux__/jammy/latest",
                  CRAN = "https://cloud.r-project.org"))
options(Ncpus = 16, timeout = 1800)
suppressMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(version = "3.23", ask = FALSE, update = FALSE)
})
cran_pkgs <- c("data.table","pheatmap","ComplexHeatmap","ggplot2","ggpubr","patchwork","cowplot",
  "viridis","RColorBrewer","corrplot","factoextra","gridExtra","scales","ggraph","igraph",
  "UpSetR","ggalluvial","circlize","survival","survminer","rms","timeROC","glmnet","maxstat",
  "clusterProfiler","immunedeconv","presto","Seurat","harmony","monocle3","SeuratObject","uwot",
  "RcppAnnoy","RcppHNSW","irlba","tidyverse","Matrix","tibble","dplyr","stringr","readr","purrr",
  "forcats","doParallel","foreach","pbapply","XML","RCurl","httr","jsonlite","openxlsx","writexl",
  "readxl","fmsb","naniar","VIM","mice","pROC","caret","randomForest","xgboost","rpart","ranger",
  "glmnetUtils","Hmisc","rstatix","ggsci","wesanderson","ggrepel","ggforce","tidygraph","graphlayouts",
  "networkD3","visNetwork","threejs","plotly","htmlwidgets","reshape2","Rtsne","umap","devtools",
  "remotes","bio3d","genefilter","sva","limma")
bioc_pkgs <- c("limma","DESeq2","edgeR","GSVA","TCGAbiolinks","GEOquery","SummarizedExperiment",
  "S4Vectors","Biobase","GenomicRanges","SingleCellExperiment","scater","scran","SingleR","celldex",
  "org.Hs.eg.db","ReactomePA","enrichplot","AUCell","RcppML","batchelor","bluster","scDblFinder",
  "dittoSeq","tradeSeq","slingshot","Biostrings","GenomicFeatures","AnnotationDbi","multtest",
  "GSVA","GSEABase","fgsea","DOSE")
installed <- rownames(installed.packages())
need_cran <- setdiff(cran_pkgs, installed)
cat("CRAN to install:", length(need_cran), "\n")
if (length(need_cran)) install.packages(need_cran)
need_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
cat("Bioconductor to install:", length(need_bioc), "\n")
if (length(need_bioc)) BiocManager::install(need_bioc, ask = FALSE, update = FALSE)
cat("DONE_ALL_INSTALL\n")
