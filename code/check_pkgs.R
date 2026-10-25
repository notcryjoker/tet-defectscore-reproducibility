pkgs <- c("DESeq2","edgeR","GSVA","TCGAbiolinks","GEOquery","org.Hs.eg.db",
          "survival","survminer","clusterProfiler","ComplexHeatmap","limma")
ok <- all(vapply(pkgs, function(p) requireNamespace(p, quietly=TRUE), logical(1)))
if (!ok) {
  missing <- pkgs[!vapply(pkgs, function(p) requireNamespace(p, quietly=TRUE), logical(1))]
  cat("MISSING:", paste(missing, collapse=", "), "\n")
}
quit(status = ifelse(ok, 0, 1))
