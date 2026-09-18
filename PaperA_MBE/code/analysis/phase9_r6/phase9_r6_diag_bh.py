# -*- coding: utf-8 -*-
import json
from pathlib import Path
import pandas as pd

BASE = Path(r"d:\人类正选择基因项目")
PH9 = BASE / "results/phase9_hardening"
master = pd.read_csv(BASE / "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv").set_index("gene_id")
print("bh_fdr_sig dtype:", master.bh_fdr_sig.dtype, "values:", master.bh_fdr_sig.unique()[:5])

rows = []
for fp in (PH9 / "busted_rerun_baseline").glob("*.json"):
    try:
        d = json.load(open(fp, encoding="utf-8"))
        rows.append({"gene_id": fp.stem, "rerun_p": float(d["test results"]["p-value"])})
    except Exception:
        pass
rr = pd.DataFrame(rows).set_index("gene_id")
m = master[["busted_p", "bh_fdr_sig"]].join(rr, how="inner")
print("n:", len(m))
print("prod sig:", m.bh_fdr_sig.astype(bool).sum(), " rerun p<0.05:", (m.rerun_p < 0.05).sum())

import numpy as np
def bh(pvals, alpha=0.05):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals, kind="mergesort")
    ranked = pvals[order]
    q = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.clip(q, 0, 1)
    return q

q400 = bh(m.rerun_p.values)
m["q400"] = q400
m["sig400"] = q400 < 0.05
print("rerun BH400 sig:", m.sig400.sum())
print("prod sig:", m.bh_fdr_sig.astype(bool).sum())
print("agreement:", (m.sig400 == m.bh_fdr_sig.astype(bool)).mean())
# cross-tab
print(pd.crosstab(m.bh_fdr_sig.astype(bool), m.sig400))
# distribution of rerun p for prod-sig genes
ps = m[m.bh_fdr_sig.astype(bool)]
print("prod-sig genes rerun p quantiles:", np.percentile(ps.rerun_p, [50, 75, 90, 95, 99]).round(5))
print("prod-sig genes rerun p>0.025:", (ps.rerun_p > 0.025).sum())
