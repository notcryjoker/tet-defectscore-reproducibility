#!/bin/bash
# 并行下载 TCGA-THYM STAR-Counts（6 并发；每文件独立 wget -c 可断点续传）
OUT=${TET_ROOT:-.}/data/raw/TCGA-THYM/gdc
mkdir -p "$OUT"
python3 - << 'PY' > /tmp/dl_jobs.txt
import json
files = json.load(open("/tmp/thym_star_files.json"))
for f in files:
    print(f"{f['file_id']}\t{f['file_name']}")
PY
echo "jobs: $(wc -l < /tmp/dl_jobs.txt)"
cat /tmp/dl_jobs.txt | xargs -P 6 -I{} bash -c '
  IFS=$'"'"'\t'"'"' read -r fid fname <<< "{}"
  d="'"$OUT"'/$fid"; p="$d/$fname"
  mkdir -p "$d"
  if [ -s "$p" ]; then exit 0; fi
  wget -q -c --tries=8 --timeout=90 -O "$p" "https://api.gdc.cancer.gov/data/$fid" || { [ -s "$p" ] || rm -f "$p"; exit 1; }
'
echo "=== final: $(find '"$OUT"' -name "*.tsv" | wc -l)/122 files ==="
du -sh "$OUT"
