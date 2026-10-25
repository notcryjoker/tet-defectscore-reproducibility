#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 Markdown 稿件中的 [@key] 引用转换为 MDPI 数字顺序编码 [n]，
并按首次引用顺序重排参考文献表。

用法: python3 convert_citations.py <input.md> <output.md>

规则:
- 引用组 [@a; @b] -> [1, 2]（按组内出现顺序编号）
- 单条引用 [@a]   -> [1]
- 文献表条目 [@key] text -> n. text（按编号升序）
- 删除"按引用键字母序排列"的说明性引用块
"""
import re
import sys

REF_HEAD = re.compile(r'^## (References|参考文献)\s*$', re.M)
ENTRY = re.compile(r'^\[@([A-Za-z0-9_\-]+)\]\s+(.*)$', re.M)
CITE_GROUP = re.compile(r'\[([^\[\]]*?@[A-Za-z0-9_\-][^\[\]]*?)\]')


def split_sections(text):
    m = REF_HEAD.search(text)
    if not m:
        raise SystemExit('未找到参考文献节标题')
    head = text[:m.start()]
    rest = text[m.end():]
    # 参考文献节从标题后到下一个 --- 分隔线或文末
    m2 = re.search(r'^---\s*$', rest, re.M)
    if m2:
        ref_body = rest[:m2.start()]
        tail = rest[m2.start():]
    else:
        ref_body, tail = rest, ''
    return head, m.group(1), ref_body, tail


def parse_entries(ref_body):
    entries = {}
    order = []
    for m in ENTRY.finditer(ref_body):
        key, body = m.group(1), m.group(2).strip()
        if key in entries:
            raise SystemExit('文献表重复定义: @%s' % key)
        entries[key] = body
        order.append(key)
    if not entries:
        raise SystemExit('文献表未解析到任何条目')
    return entries, order


def convert(path_in, path_out):
    text = open(path_in, encoding='utf-8').read()
    head, ref_title, ref_body, tail = split_sections(text)
    entries, alpha_order = parse_entries(ref_body)

    # 1) 按正文首次出现顺序编号
    numbering = {}
    next_n = 1
    warnings = []

    def repl(m):
        nonlocal next_n
        inner = m.group(1)
        keys = re.findall(r'@([A-Za-z0-9_\-]+)', inner)
        nums = []
        for k in keys:
            if k not in entries:
                warnings.append('正文引用 @%s 在文献表中不存在' % k)
                continue
            if k not in numbering:
                numbering[k] = next_n
                next_n += 1
            nums.append(numbering[k])
        if not nums:
            return m.group(0)
        return '[' + ', '.join(str(n) for n in nums) + ']'

    new_head = CITE_GROUP.sub(repl, head)

    # 2) 重建文献表（按编号升序；未被引用的条目排在末尾并告警）
    unreferenced = [k for k in alpha_order if k not in numbering]
    for k in unreferenced:
        warnings.append('文献 @%s 未被正文引用（保留在表末）' % k)
    cited_sorted = sorted(numbering, key=lambda k: numbering[k])
    lines = ['## ' + ref_title, '']
    for k in cited_sorted:
        lines.append('%d. %s' % (numbering[k], entries[k]))
        lines.append('')
    for k in unreferenced:
        lines.append('%d. %s' % (next_n, entries[k]))
        lines.append('')
        next_n += 1
    new_ref = '\n'.join(lines).rstrip() + '\n'

    out = new_head.rstrip() + '\n\n---\n\n' + new_ref + tail
    open(path_out, 'w', encoding='utf-8').write(out)

    # 3) 校验
    leftover = re.findall(r'\[@[^\]]*\]', out)
    if leftover:
        raise SystemExit('转换后仍残留 [@key] 引用: %s' % leftover[:5])
    if warnings:
        print('警告:', file=sys.stderr)
        for w in warnings:
            print(' -', w, file=sys.stderr)
    print('完成: %s' % path_out)
    print('  引用文献 %d 条（编号 1-%d）' % (len(cited_sorted), len(cited_sorted)))
    if unreferenced:
        print('  未引用保留 %d 条' % len(unreferenced))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    convert(sys.argv[1], sys.argv[2])
