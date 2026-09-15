#!/usr/bin/env python3
"""
Paper A (MBE) - Figure 2 (BUSTED p-value distribution) as PDF for submission.
Mirrors the Fig 2 block in scripts/paperA_figures.py (lines 100-127).
Output: results/paper/figures_v9/Figure2_BUSTED_pvalue_distribution.{png,pdf}
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if os.path.exists("/mnt/d"):
    BASE = "/mnt/d/人类正选择基因项目"
else:
    BASE = "D:/人类正选择基因项目"

OUT = f"{BASE}/results/paper/figures_v9"
os.makedirs(OUT, exist_ok=True)

v7 = pd.read_csv(
    f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv",
    low_memory=False)

p = pd.to_numeric(v7["busted_p"], errors="coerce").dropna().values
q = pd.to_numeric(v7["bh_fdr"], errors="coerce").dropna().values
p_max_sig = p[q < 0.05].max() if (q < 0.05).any() else np.nan

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
ax = axes[0]
ax.hist(p, bins=50, color="#2166ac", edgecolor="white", linewidth=0.3)
ax.axvline(p_max_sig, color="#d6604d", ls="--", lw=1.3,
           label=f"max p with FDR<0.05 ({p_max_sig:.2e})")
ax.set_xlabel("BUSTED p-value (bounded at 0.5)")
ax.set_ylabel("Genes")
ax.set_title(f"(a) BUSTED p-value distribution (n={len(p):,})")
ax.legend(fontsize=8)

ax = axes[1]
qc = np.sort(q)
ax.plot(np.arange(1, len(qc) + 1), qc, color="#2166ac", lw=1.2, label="sorted q-values")
ax.axhline(0.05, color="#d6604d", ls="--", lw=1.3, label="FDR = 0.05")
ax.set_xlabel("Gene rank (sorted by q)")
ax.set_ylabel("BH FDR q-value")
ax.set_title(f"(b) FDR q-values: {(q < 0.05).sum():,} / {len(q):,} significant")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/Figure2_BUSTED_pvalue_distribution.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Figure2_BUSTED_pvalue_distribution.pdf", bbox_inches="tight")
plt.close(fig)
print("saved Figure2_BUSTED_pvalue_distribution.png/.pdf")
