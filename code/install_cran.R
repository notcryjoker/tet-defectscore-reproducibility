options(repos = c(PPM = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
options(Ncpus = 16, timeout = 1200)
pkgs <- c("data.table","pheatmap","ComplexHeatmap","ggplot2","ggpubr","patchwork","cowplot",
"viridis","RColorBrewer","corrplot","factoextra","gridExtra","scales","ggraph","igraph",
"UpSetR","ggalluvial","circlize","survival","survminer","rms","timeROC","glmnet","maxstat",
"clusterProfiler","immunedeconv","presto","Seurat","harmony","monocle3","SeuratObject","uwot",
"RcppAnnoy","RcppHNSW","irlba","tidyverse","Matrix","tibble","dplyr","stringr","readr","purrr",
"forcats","doParallel","foreach","pbapply","XML","RCurl","httr","jsonlite","openxlsx","writexl",
"readxl","BiocManager","fmsb","naniar","VIM","mice","pROC","caret","randomForest","xgboost",
"rpart","ranger","glmnetUtils","Hmisc","rstatix","ggsci","wesanderson","ggrepel","ggforce",
"tidygraph","graphlayouts","networkD3","visNetwork","threejs","plotly","htmlwidgets")
need <- pkgs[!pkgs %in% rownames(installed.packages())]
cat("to install:", length(need), "\n")
if (length(need)) install.packages(need)
cat("DONE CRAN INSTALL\n")
