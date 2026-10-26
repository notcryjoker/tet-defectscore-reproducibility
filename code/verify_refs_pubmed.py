#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
撤稿库普查 + PMID/DOI 补全核验（NCBI E-utilities）
对稿件文献表中缺 DOI/PMID 的条目逐条检索 PubMed：
  1) esearch 定位 PMID（标题词 + 第一作者 + 年份）
  2) esummary 取回 DOI、卷期页、出版类型
  3) 交叉验证：返回记录的标题须与稿件条目标题一致（防止张冠李戴）
  4) 检查 Retracted Publication / Expression of Concern 标记
输出: project/audit_report/ref_verification_pubmed.md
"""
import json
import re
import time
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
MANUSCRIPT = os.environ.get("MS_PATH", "manuscript_draft.md")
OUT = os.environ.get("REF_OUT", "project/audit_report/ref_verification_pubmed.md")

# 稿件条目 -> 检索词（标题词 + 第一作者姓 + 年份）；仅收录缺 DOI/PMID 的条目
TARGETS = [
    ("engels2010",       "Engels EA", "Epidemiology of thymoma and associated malignancies", 2010),
    ("diagnostics2023",  "", "Thymoma an overview Diagnostics", 2023),
    ("meriggioli2009",   "Meriggioli MN", "Autoimmune myasthenia gravis emerging clinical and biological heterogeneity", 2009),
    ("gilhus2016",       "Gilhus NE", "Myasthenia gravis autoantibody characteristics and their implications for therapy", 2016),
    ("berrihaknin2014",  "Berrih-Aknin S", "Diagnostic and clinical classification of autoimmune myasthenia gravis", 2014),
    ("suzuki2013",       "Suzuki S", "Overlooked non-motor symptoms in myasthenia gravis", 2013),
    ("takaba2015",       "Takaba H", "Fezf2 orchestrates a thymic program of self-antigen expression for immune tolerance", 2015),
    ("yasumizu2022",     "Yasumizu Y", "Single-cell atlas of thymic epithelial tumors and the thymic ecosystem", 2022),
    ("benitez2021",      "Benitez JC", "Cancer activation pathways of thymic epithelial tumors", 2021),
    ("bornstein2018",    "Bornstein C", "Single-cell mapping of the thymic stroma identifies IL-25-producing tuft epithelial cells", 2018),
    ("hao2021",          "Hao Y", "Integrated analysis of multimodal single-cell data", 2021),
    ("huber2015",        "Huber W", "Orchestrating high-throughput genomic analysis with Bioconductor", 2015),
    ("iorio2024",        "Iorio R", "Myasthenia gravis the changing treatment landscape in the era of molecular therapies", 2024),
    ("jclinpathol2023",  "", "Unique correlation between GTF2I mutation and spindle cell morphology in thymomas", 2023),
    ("kashima2021",      "Kashima J", "CD70 in thymic squamous cell carcinoma potential diagnostic markers", 2021),
    ("koga1994",         "Koga K", "A review of 79 thymomas modification of staging system", 1994),
    ("lungcancer2017",   "", "GTF2I mutation frequently occurs in more indolent thymic epithelial tumors", 2017),
    ("masaoka1981",      "Masaoka A", "Follow-up study of thymomas with special reference to their clinical stages", 1981),
    ("miller2018",       "Miller CN", "Thymic tuft cells promote an IL-4-enriched medulla and shape thymocyte development", 2018),
    ("robinson2010",     "Robinson MD", "edgeR a Bioconductor package for differential expression analysis", 2010),
    ("shelly2011",       "Shelly S", "Thymoma and autoimmunity", 2011),
    ("thongprayoon2013", "Thongprayoon C", "Associated immunological disorders and cellular immune dysfunction in thymoma", 2013),
    ("yamada2021",       "Yamada Y", "A tuft cell-like signature is highly prevalent in thymic squamous cell carcinoma", 2021),
]


def get(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception:
            time.sleep(2)
    return None


def esearch(term):
    q = urllib.parse.quote(term)
    d = get(BASE + f"esearch.fcgi?db=pubmed&term={q}&retmode=json&retmax=3")
    if not d:
        return []
    return d.get("esearchresult", {}).get("idlist", [])


def esummary(pmids):
    d = get(BASE + f"esummary.fcgi?db=pubmed&id={','.join(pmids)}&retmode=json")
    if not d:
        return {}
    return d.get("result", {})


def norm(s):
    return re.sub(r"[^a-z0-9 ]", " ", s.lower())


def title_match(record_title, expect):
    """返回记录的标题是否与期望标题一致（关键词覆盖率）"""
    rt = set(norm(record_title).split())
    et = [w for w in norm(expect).split() if len(w) > 3]
    if not et:
        return True
    hit = sum(1 for w in et if w in rt)
    return hit / len(et) >= 0.8


# 读取稿件条目文本用于交叉验证
text = open(MANUSCRIPT, encoding="utf-8").read()
ref_sec = text[text.find("## References"):]
entries = dict(re.findall(r"^(\d+)\.\s+(.*)$", ref_sec, re.M))
key_by_num = {}
for num, body in entries.items():
    m = re.match(r"\[@([A-Za-z0-9_\-]+)\]", body)
    if m:
        key_by_num[num] = m.group(1)
entry_text = {}
for num, body in entries.items():
    entry_text[key_by_num.get(num, num)] = body

rows = []
for key, author, title_words, year in TARGETS:
    time.sleep(0.4)
    term = title_words + (f" AND {author}[Author]" if author else "") + f" AND {year}[pdat]"
    ids = esearch(term)
    result = {"key": key, "pmid": None, "doi": None, "match": False, "flag": "", "found_title": "", "vol": "", "issue": "", "pages": ""}
    if ids:
        s = esummary(ids)
        for pid in ids:
            rec = s.get(pid, {})
            t = rec.get("title", "")
            if title_match(t, title_words):
                result["pmid"] = pid
                result["found_title"] = t
                result["match"] = True
                for aid in rec.get("articleids", []):
                    if aid.get("idtype") == "doi":
                        result["doi"] = aid.get("value")
                result["vol"] = rec.get("source", "")
                pts = rec.get("pubtype", [])
                if any("Retracted" in p for p in pts):
                    result["flag"] = "RETRACTED PUBLICATION!"
                elif any("Expression of Concern" in p for p in pts):
                    result["flag"] = "EXPRESSION OF CONCERN!"
                else:
                    result["flag"] = "clean"
                break
        else:
            # 取第一条做失配报告
            rec = s.get(ids[0], {})
            result["found_title"] = rec.get("title", "")
            result["flag"] = "NO MATCH (需人工核对)"
    else:
        result["flag"] = "NOT FOUND (需人工核对)"
    time.sleep(0.4)
    rows.append(result)

lines = [
    "# 文献核验记录 — PubMed E-utilities 撤稿库普查 + PMID/DOI 补全",
    "",
    "**日期**：2026-10-23  |  **工具**：NCBI E-utilities（esearch + esummary）",
    "**范围**：稿件文献表中缺 DOI/PMID 的 23 条条目（终审遗留事项：撤稿库普查收尾）",
    "**方法**：标题词+第一作者+年份检索 → PMID → DOI/出版类型回取 → 标题交叉验证（防张冠李戴）→ Retracted/EOI 标记检查",
    "",
    "| # | 引用键 | PMID | DOI | 标题交叉验证 | 撤稿/关注状态 |",
    "|---|---|---|---|---|---|",
]
n_ok = n_flag = 0
for i, r in enumerate(rows, 1):
    ok = r["match"] and r["flag"] == "clean"
    n_ok += ok
    if r["flag"] not in ("clean",):
        n_flag += 1
    lines.append("| %d | %s | %s | %s | %s | %s |" % (
        i, r["key"],
        r["pmid"] or "—",
        r["doi"] or "—",
        "一致" if r["match"] else "**失配**",
        r["flag"],
    ))
lines += [
    "",
    "**汇总**：%d/%d 条核验通过（PMID 定位 + 标题一致 + 无撤稿/关注标记）。" % (n_ok, len(rows)),
    "",
    "说明：",
    "- 标题交叉验证阈值：期望标题关键词（长度>3）在返回记录标题中的覆盖率 ≥80% 方判一致。",
    "- 'clean' = PubMed 出版类型中无 Retracted Publication / Expression of Concern 标记。",
    "- 失配或 NOT FOUND 条目需人工在 PubMed 官网复核，不得凭记忆填写。",
    "- 已核验 DOI/PMID 的既有条目（radovich2018、malchow2016、fenioux2023 等）不在本次范围；其更正记录见重审报告 N2（均为良性更正，非撤稿）。",
]
open(OUT, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n".join(lines[:12]))
print("...")
print("汇总: %d/%d 通过, %d 条需关注" % (n_ok, len(rows), n_flag))
for r in rows:
    if r["flag"] != "clean" or not r["match"]:
        print("  [%s] pmid=%s flag=%s found=%s" % (r["key"], r["pmid"], r["flag"], r["found_title"][:80]))
