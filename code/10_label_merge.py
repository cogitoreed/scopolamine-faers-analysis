#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10_label_merge.py —— 把 v2 主分析（22b 核心 PT 总体表）的 a/ROR 并入 label 编码草表
输出：材料包_EN/tab_label.csv（PT｜FAERS ROR(95%CI)｜a｜Label mention｜Label source｜snippet）
"""
import csv

draft = list(csv.DictReader(open("材料包_EN/tab_label_draft.csv", encoding="utf-8-sig")))
core = {}
for row in csv.DictReader(open("output_v2/22b_核心PT_总体.csv", encoding="utf-8")):
    core[row["pt"].strip().upper()] = row

with open("材料包_EN/tab_label.csv", "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.writer(fh)
    w.writerow(["PT", "FAERS a", "FAERS ROR (95% CI)", "Label mention (Y/N/U)", "Label source", "Label text snippet", "note"])
    merged, label_only = 0, []
    for r in draft:
        pt = r["PT"].strip().upper()
        c = core.get(pt)
        if c:
            a = c["a"]
            ror = f"{float(c['ROR']):.2f} ({float(c['lcl']):.2f}-{float(c['ucl']):.2f})"
            merged += 1
        else:
            a, ror = "", ""
            label_only.append(pt)
        w.writerow([pt, a, ror, r["Label mention (Y/N/U)"], "DailyMed TRANSDERM SCŌP SPL b877a694-a1d0-4280-937a-a06820b12a88", r["Label text snippet"], r["note"]])
print(f"[完成] 材料包_EN/tab_label.csv（并入 {merged} 条 PT 的 a/ROR）")
if label_only:
    print("仅 label 表有（未进核心 PT 集）:", ", ".join(label_only))
