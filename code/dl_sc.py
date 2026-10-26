import urllib.request, os, time, sys
base="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284479/suppl"
out=os.path.join(os.environ.get("TET_ROOT", "."), "data/raw/GSE284479/scRNA")
os.makedirs(out, exist_ok=True)
keep=lambda fn: (any(k in fn for k in ["matrix.mtx.gz","barcodes.tsv.gz","features.tsv.gz"]) or "contig_annotations" in fn)
targets=[]
for line in open("/tmp/gse284479_fl.txt"):
    if not line.startswith("File"): continue
    fn=line.rstrip("\n").split("\t")[1]
    if keep(fn): targets.append(fn)
print(f"targets: {len(targets)}", flush=True)
ok=fail=skip=0
for i,fn in enumerate(targets,1):
    path=os.path.join(out,fn)
    if os.path.exists(path) and os.path.getsize(path)>0:
        skip+=1; continue
    url=f"{base}/{fn}"
    try:
        urllib.request.urlretrieve(url, path)
        ok+=1
        if i%5==0: print(f"[{i}/{len(targets)}] ok={ok} fail={fail} skip={skip} last={fn} size={os.path.getsize(path)}", flush=True)
    except Exception as e:
        fail+=1; print(f"FAIL {fn}: {e}", flush=True)
print(f"DONE ok={ok} fail={fail} skip={skip}", flush=True)
