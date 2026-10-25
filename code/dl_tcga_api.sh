#!/bin/bash
set -u
OUT=/root/.codebuddy/artifact/2f8a7d30-9c21-4f0b-1d53-0f0283e72000/data/raw/TCGA-THYM/gdc
mkdir -p "$OUT"
python3 - << 'PY'
import json, subprocess, os
files = json.load(open("/tmp/thym_star_files.json"))
out = "/root/.codebuddy/artifact/2f8a7d30-9c21-4f0b-1d53-0f0283e72000/data/raw/TCGA-THYM/gdc"
ok = skip = fail = 0
for i, f in enumerate(files, 1):
    d = os.path.join(out, f["file_id"])
    path = os.path.join(d, f["file_name"])
    os.makedirs(d, exist_ok=True)
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        skip += 1; continue
    url = f"https://api.gdc.cancer.gov/data/{f['file_id']}"
    r = subprocess.run(["wget", "-q", "-c", "--tries=6", "--timeout=90", "-O", path, url])
    if r.returncode == 0 and os.path.exists(path) and os.path.getsize(path) > 1000:
        ok += 1
        if i % 10 == 0: print(f"[{i}/{len(files)}] ok={ok} skip={skip} fail={fail}", flush=True)
    else:
        fail += 1
        if os.path.exists(path) and os.path.getsize(path) <= 1000: os.remove(path)
        print(f"FAIL {f['file_name']} rc={r.returncode}", flush=True)
print(f"DONE ok={ok} skip={skip} fail={fail}", flush=True)
PY
echo "=== total size ==="; du -sh "$OUT"
