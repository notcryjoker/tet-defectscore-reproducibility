#!/bin/bash
# 全流程驱动：R 包就绪后依次执行 01→08 + 03b(scRNA) + 09(Fig1)
# 任一步失败不中断，全部记录于 results/pipeline.log
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(dirname "$SCRIPT_DIR")"   # deposit 根目录（code/ 的父目录）
export TET_ROOT="${TET_ROOT:-$BASE}"
PL=$BASE/results/pipeline.log
echo "[driver] start $(date)" > $PL

pkgs_ok() { Rscript $BASE/code/check_pkgs.R >> $PL 2>&1; }

for i in $(seq 1 60); do
  if pkgs_ok; then echo "[driver] R READY at iter $i ($(date))" >> $PL; break; fi
  echo "[driver] iter $i waiting ($(date))" >> $PL
  sleep 60
done
if ! pkgs_ok; then echo "[driver] FATAL TIMEOUT ($(date))" >> $PL; exit 1; fi

cd $BASE
for step in 01_tcga_prepare 02_scores 03_associations 04_survival 05_external 06_icb 07_figures 09_fig1_design 08_results_fill; do
  echo "===== $step START $(date) =====" >> $PL
  if Rscript code/$step.R >> $PL 2>&1; then
    echo "[driver] $step OK $(date)"" >> $PL
  else
    rc=$?
    echo "[driver] $step FAILED rc=$rc ($(date)) -- continuing" >> $PL
  fi
done
echo "[driver] PIPELINE_DONE $(date)" >> $PL
