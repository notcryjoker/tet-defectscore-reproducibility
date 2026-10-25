#!/bin/bash
set -e
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284479/suppl"
OUT="/root/.codebuddy/artifact/2f8a7d30-9c21-4f0b-1d53-0f0283e72000/data/raw/GSE284479/scRNA"
mkdir -p "$OUT"
cd "$OUT"
# parse filelist: download MTX, barcodes, features, and TCR contig annotations
python3 - << 'PY'
import urllib.request, os
base="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284479/suppl"
out="/root/.codebuddy/artifact/2f8a7d30-9c21-4f0b-1d53-0f0283e72000/data/raw/GSE284479/scRNA"
os.makedirs(out, exist_ok=True)
for line in open("/tmp/gse284479_fl.txt"):
    if not line.startswith("File"): continue
    parts=line.rstrip("\n").split("\t")
    fn=parts[1]
    ext=fn.split(".")[-1]
    # keep MTX matrix, barcodes, features, and TCR contig csv
    if any(k in fn for k in ["matrix.mtx.gz","barcodes.tsv.gz","features.tsv.gz"]) or "contig_annotations" in fn:
        url=f"{base}/{fn}"
        if not os.path.exists(os.path.join(out,fn)):
            print("GET", fn)
PY
# now actually download (re-run python to fetch)
python3 - << 'PY'
import urllib.request, os, sys
base="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284479/suppl"
out="/root/.codebuddy/artifact/2f8a7d30-9c21-4f0b-1d53-0f0283e72000/data/raw/GSE284479/scRNA"
os.makedirs(out, exist_ok=True)
n=0
for line in open("/tmp/gse284479_fl.txt"):
    if not line.startswith("File"): continue
    parts=line.rstrip("\n").split("\t")
    fn=parts[1]
    if any(k in fn for k in ["matrix.mtx.gz","barcodes.tsv.gz","features.tsv.gz"]) or "contig_annotations" in fn:
        path=os.path.join(out,fn)
        if os.path.exists(path) and os.path.getsize(path)>0: 
            continue
        url=f"{base}/{fn}"
        try:
            urllib.request.urlretrieve(url, path)
            n+=1
        except Exception as e:
            print("FAIL", fn, e)
print("downloaded", n, "new files")
PY
echo "=== total size ==="; du -sh "$OUT"
