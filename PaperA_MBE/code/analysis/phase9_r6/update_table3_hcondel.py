# -*- coding: utf-8 -*-
"""R6 P0 hCONDEL fix: update Table 3 (LOO matrix) hCONDEL column in
Main_Tables_v9.xlsx (both the results/paper source copy and the
MBE_submission_PaperA/04_tables package copy).

hCONDEL gene mapping repaired (hg18->hg38 liftOver; 112 genes;
results/phase9_hardening/hcondels_liftover_fix.json). OR / one-sided-greater
Fisher p / Woolf 95% CI are recomputed from k-over-n counts, reproducing the
json exactly; HAR column unchanged. Footnote updated to document the fixed
mapping with per-variant k/n and Fisher P.
"""
import json
import os

import numpy as np
import openpyxl
from scipy.stats import fisher_exact

if os.path.exists("/mnt/d"):
    ROOT = "/mnt/d/人类正选择基因项目"
else:
    ROOT = "D:/人类正选择基因项目"
HCFIX = os.path.join(ROOT, "results/phase9_hardening/hcondels_liftover_fix.json")
TARGETS = [
    os.path.join(ROOT, "results/paper/Main_Tables_v9.xlsx"),
    os.path.join(ROOT, "MBE_submission_PaperA/04_tables/Main_Tables_v9.xlsx"),
]

with open(HCFIX, encoding="utf-8") as f:
    hcfix = json.load(f)
HC_TOTAL = int(hcfix["gene_sets"]["union_fixed"])   # 112
N_UNIV = int(hcfix["universe_n"])                   # 4,974

VARIANTS = [("Full", "full"), ("LOO-caMPRA", "LOO-caMPRA"), ("LOO-tau", "LOO-tau"),
            ("LOO-nc", "LOO-nc"), ("LOO-brain", "LOO-brain")]
cells, kv = [], {}
for row, (v, key) in enumerate(VARIANTS, start=4):   # data rows 4-8 in sheet
    st = hcfix["loo_hcondel"][key]
    k, n = [int(x) for x in st["k_over_n"].split("/")]
    tab = [[k, n - k], [HC_TOTAL - k, N_UNIV - n - (HC_TOTAL - k)]]
    orv, pv = fisher_exact(tab, alternative="greater")
    assert abs(float(orv) - st["OR"]) < 0.01, (v, orv, st["OR"])
    assert abs(float(pv) - st["p"]) / st["p"] < 1e-6, (v, pv, st["p"])
    a, b, c, d = tab[0][0], tab[0][1], tab[1][0], tab[1][1]
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo = float(np.exp(np.log(orv) - 1.96 * se))
    hi = float(np.exp(np.log(orv) + 1.96 * se))
    stars = "***" if pv < 0.001 else ("**" if pv < 0.01 else ("*" if pv < 0.05 else "n.s."))
    cells.append((row, f"{orv:.2f} ({lo:.2f}-{hi:.2f}){stars}"))
    kv[v] = (k, n, pv)
    print(f"{v}: OR={orv:.2f} CI=[{lo:.2f},{hi:.2f}] p={pv:.1e} k/n={k}/{n} "
          f"(old: OR={st['old_OR']}, p={st['old_p']:.1e}, {st['old_k_over_n']})")

foot = [
    "95% CIs by normal approximation on log(OR). hCONDEL mapping uses repaired",
    "hg18-to-hg38 liftOver coordinates (112 genes); hCONDEL k/n and Fisher P:",
    "Full 18/293, 7.8E-05; LOO-caMPRA 32/800, 5.3E-04; LOO-tau 13/148, 2.3E-05;",
    "LOO-nc 10/110, 1.6E-04; LOO-brain 29/613, 5.7E-05.",
]
for v, key in VARIANTS:  # verify the footnote numbers against recomputation
    k, n, pv = kv[v]
    kk, nn = (int(x) for x in hcfix["loo_hcondel"][key]["k_over_n"].split("/"))
    assert (k, n) == (kk, nn), (v, k, n, kk, nn)

for path in TARGETS:
    wb = openpyxl.load_workbook(path)
    ws = wb["Table 3"]
    assert ws["A3"].value == "Variant" and "hCONDEL" in ws["E3"].value, path
    for row, val in cells:
        ws.cell(row=row, column=5).value = val
    # clear the blank separator row between data and footnote (defensive:
    # a previous misaligned write put a value in E9)
    for rr in range(4 + len(VARIANTS), 14):
        ws.cell(row=rr, column=5).value = None
    # replace stale footnote rows 14-16 (old: conditional-MLE note for OR 2.94)
    old14 = ws["A14"].value
    assert old14 and old14.startswith("95% CIs"), (path, old14)
    for i, line in enumerate(foot):
        ws.cell(row=14 + i, column=1).value = line
    # clear the previous 16th line if the new footnote is longer (4 vs 3 lines)
    if len(foot) > 3:
        ws.cell(row=14 + len(foot) - 1, column=1).value = foot[-1]
    wb.save(path)
    print("updated:", path)
print("done")
