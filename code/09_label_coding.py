#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09_label_coding.py —— S0.8C Label comparison 三档编码（v3.1 D 组）
对照 Transderm Scōp FDA 说明书原文（DailyMed SPL，已抓取存档），
对核心 PT 逐条编码 Yes/No/Unclear + 原文摘录。
【人核验点】输出为 tab_label_draft.csv，每条 snippet 须人工对照原文确认后才可入稿。
Yes     = 说明书明确列出该 PT 或其直接等价表述
Unclear = 仅提及上游症状/机制性表述（如"胃肠动力下降"），未明确该 PT
No      = 说明书未检索到相关表述
"""
import re, csv

SRC = "材料包_EN/label/transderm_scop_plain.txt"
OUT = "材料包_EN/tab_label_draft.csv"

text = open(SRC, encoding="utf-8").read()
# 规范空白，保留句读
text_norm = re.sub(r"\s+", " ", text)

# PT → 匹配词（正则，忽略大小写）；unclear_pat = 机制性/上游表述（命中则 Unclear）
SPEC = [
    # (PT, yes_pat, unclear_pat)
    ("DRY MOUTH",              r"dry mouth|xerostomia", None),
    ("SOMNOLENCE",             r"drowsiness|somnolence|sleepiness", None),
    ("DISORIENTATION",         r"disorientation|disoriented", None),
    ("CONFUSIONAL STATE",      r"confusion(?! of)|confusional|confused", None),
    ("HALLUCINATION",          r"hallucinat", None),
    ("DELIRIUM",               r"delirium|acute toxic psychosis|toxic psychosis", r"psychosis"),
    ("AGITATION",              r"agitation", None),
    ("AMNESIA",                r"amnesia|memory (loss|impairment)", None),
    ("DIZZINESS",              r"dizziness|dizzy|lightheaded", None),
    ("VISION BLURRED",         r"blurred vision|blurry vision", None),
    ("MYDRIASIS",              r"mydriasis|dilat(ion|ed) (of )?the pupil|pupillary dilation", None),
    ("URINARY RETENTION",      r"urinary retention|difficulty in urination|urine flow", r"urinary"),
    ("TACHYCARDIA",            r"tachycardia|palpitation", r"heart rate"),
    ("DRY EYE",                r"dry eye", None),
    ("ACCOMMODATION DISORDER", r"accommodation", None),
    ("ANGLE CLOSURE GLAUCOMA", r"angle closure glaucoma", r"glaucoma"),
    ("SEIZURE",                r"seizure", None),
    ("HYPERTHERMIA",           r"hyperthermia", r"heat"),
    ("HEAT STROKE",            r"heat stroke", r"heat"),
    ("HEAT EXHAUSTION",        r"heat exhaustion", r"heat"),
    ("HEAT CRAMPS",            r"heat cramps", r"heat"),
    ("PYREXIA",                r"pyrexia|fever", None),
    ("ANGIOEDEMA",             r"angioedema", r"swell"),
    ("HYPERSENSITIVITY",       r"hypersensitivity|allergic reaction|anaphyla", None),
    ("ANAPHYLACTIC REACTION",  r"anaphyla", None),
    ("URTICARIA",              r"urticaria|hives", r"rash"),
    ("RASH",                   r"rash(?!-)", None),
    ("PRURITUS",               r"pruritus|itch", None),
    ("ERYTHEMA",               r"erythema|redness", None),
    ("DYSPHONIA",              r"dysphonia|hoarseness|voice disorder", None),
    ("SWOLLEN TONGUE",         r"swollen tongue|tongue swelling", r"swell"),
    ("SWELLING FACE",          r"swelling of the face|facial swelling|swollen face", r"swell"),
    ("NAUSEA",                 r"nausea", None),
    ("VOMITING",               r"vomit", None),
    ("HEADACHE",               r"headache", None),
    ("DIARRHOEA",              r"diarrhea|diarrhoea", None),
    ("ABDOMINAL PAIN",         r"abdominal pain|stomach pain", r"gastrointestinal motility|decrease.{0,30}motility"),
    ("CONSTIPATION",           r"constipation", r"gastrointestinal motility|decrease.{0,30}motility"),
    ("DYSuria",                r"dysuria", None),
    ("APPLICATION SITE REACTION", r"application site|(skin )?site reaction", None),
    ("PRODUCT ADHESION ISSUE", r"adhesion|does not stick|falls off", None),
    ("DERMATITIS CONTACT",     r"contact dermatitis|skin (irritation|reaction)", r"skin"),
    ("MALAISE",                r"malaise|fatigue|weakness", None),
    ("DYSPNOEA",               r"dyspnea|dyspnoea|shortness of breath|breathing difficult", None),
    ("DEATH",                  r"death|fatal", None),
    ("PNEUMONIA",              r"pneumonia", None),
    ("FALL",                   r"\bfalls?\b", None),
    ("ANXIETY",                r"anxiety", None),
    ("DRUG INTOLERANCE",       r"drug intolerance|intolerant", r"withdrawal|post-removal"),
    ("DRUG WITHDRAWAL SYMPTOMS", r"withdrawal|post-removal", None),
    ("WEIGHT DECREASED",       r"weight (loss|decrease)", None),
]

def snippet(pat, width=170):
    m = re.search(pat, text_norm, re.I)
    if not m: return None, None
    s = max(0, m.start() - width//2); e = min(len(text_norm), m.end() + width//2)
    return text_norm[s:e].strip(), m.start()

rows = []
for pt, yes_pat, unc_pat in SPEC:
    sn, pos = snippet(yes_pat)
    if sn:
        rows.append([pt.upper(), "Yes", yes_pat, sn, "明确提及，待人核对原文"])
        continue
    if unc_pat:
        sn2, _ = snippet(unc_pat)
        if sn2:
            rows.append([pt.upper(), "Unclear", unc_pat, sn2, "仅机制性/上游表述，待人裁定 Yes 或 No"])
            continue
    rows.append([pt.upper(), "No", "", "", "说明书未检索到相关表述（待人复核，可改用其他同义词再检索）"])

with open(OUT, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.writer(fh)
    w.writerow(["PT", "Label mention (Y/N/U)", "matched_term", "Label text snippet", "note"])
    w.writerows(rows)

n = {v: sum(1 for r in rows if r[1] == v) for v in ("Yes", "Unclear", "No")}
print(f"[完成] {OUT}  Yes={n['Yes']}  Unclear={n['Unclear']}  No={n['No']}")
for r in rows:
    if r[1] != "Yes":
        print(f"  {r[1]:8s} {r[0]}")
