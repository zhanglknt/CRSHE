# -*- coding: utf-8 -*-
"""Diagnostic: per-gene p agreement between production and baseline rerun."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(r"d:\人类正选择基因项目")
PH9 = BASE / "results/phase9_hardening"
master = pd.read_csv(BASE / "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv").set_index("gene_id")
subset = pd.read_csv(PH9 / "stratified_400_genes.csv")

rows = []
for fp in (PH9 / "busted_rerun_baseline").glob("*.json"):
    try:
        d = json.load(open(fp, encoding="utf-8"))
        p = d["test results"]["p-value"]
        rows.append({"gene_id": fp.stem, "rerun_p": float(p)})
    except Exception:
        pass
rr = pd.DataFrame(rows).set_index("gene_id")
sub = subset.drop(columns=[c for c in ["busted_p", "bh_fdr_sig", "busted_lrt", "bh_fdr"] if c in subset.columns]).set_index("gene_id")
m = sub.join(rr, how="inner").join(master[["busted_p", "bh_fdr_sig"]], how="inner")
print(f"n={len(m)}")
sig = m[m.bh_fdr_sig.astype(bool)]
nonsig = m[~m.bh_fdr_sig.astype(bool)]
print(f"sig stratum n={len(sig)}: rerun p<0.05: {(sig.rerun_p<0.05).sum()}, rerun p<0.01: {(sig.rerun_p<0.01).sum()}")
print(f"  prod p quantiles: {np.percentile(sig.busted_p,[10,50,90]).round(5)}")
print(f"  rerun p quantiles: {np.percentile(sig.rerun_p,[10,50,90]).round(5)}")
print(f"  spearman within sig stratum: {stats.spearmanr(sig.busted_p, sig.rerun_p)[0]:.3f}")
print(f"  pearson log10 within sig stratum: {stats.pearsonr(np.log10(sig.busted_p+1e-12), np.log10(sig.rerun_p+1e-12))[0]:.3f}")
print(f"  median abs log10 shift: {np.median(np.abs(np.log10(sig.rerun_p+1e-12)-np.log10(sig.busted_p+1e-12))):.3f}")
print(f"nonsig stratum n={len(nonsig)}: rerun p<0.05: {(nonsig.rerun_p<0.05).sum()}, p<0.01: {(nonsig.rerun_p<0.01).sum()}, p==0.5: {(nonsig.rerun_p==0.5).sum()}")
# where do flips happen
m["prod_sig"] = m.bh_fdr_sig.astype(bool)
both = m.dropna()
lo = both[both.busted_p < 0.5]
print(f"\nall p<0.05 genes (prod): {(both.busted_p<0.05).sum()}, of which rerun p<0.05: {((both.busted_p<0.05)&(both.rerun_p<0.05)).sum()}")
print(f"prod p in [0.05,0.5): {((both.busted_p>=0.05)&(both.busted_p<0.5)).sum()}, of which rerun p<0.05: {((both.busted_p>=0.05)&(both.busted_p<0.5)&(both.rerun_p<0.05)).sum()}")
