#!/usr/bin/env python3
"""
Paper A (MBE) — P0: Four supplementary figures
===============================================
Fig 1: Analysis pipeline flowchart (Phase 1-8)
Fig 2: BUSTED p-value distribution (with FDR threshold)
Fig 3: GDS vs RDS scatter (classification + threshold bands)
Fig 4: v5 vs v7 classification comparison (confusion + counts)

Output: results/paper/supplementary/figures/

Run (Windows): python.exe scripts/paperA_figures.py
"""
import os
import numpy as np
import pandas as pd
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = os.environ.get("HSD_BASE") or (
    "/mnt/d/人类正选择基因项目" if os.path.exists("/mnt/d") else "D:/人类正选择基因项目")

OUT_DIR = f"{BASE}/results/paper/supplementary/figures"
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8})

v7 = pd.read_csv(f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv",
                 low_memory=False)
v5 = pd.read_csv(f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v5/gene_classification_v5.csv",
                 low_memory=False)

print("=" * 70)
print("Paper A — Four supplementary figures")
print("=" * 70)

# ==================================================================
# Fig 1: Pipeline flowchart
# ==================================================================
print("\n--- Fig 1: pipeline flowchart ---")
fig, ax = plt.subplots(figsize=(10, 7.2))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, text, fc="#dbe9f6", ec="#2166ac", fs=8.5, lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                                fc=fc, ec=ec, lw=lw))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, wrap=True)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=12, color="#444444", lw=1.0))

# Row 1: input
box(0.3, 8.7, 2.9, 0.95, "10 primate genomes\n13,509–23,485 CDS/species\n(Ensembl)", fc="#e8f4e8", ec="#33a02c")
box(3.6, 8.7, 2.9, 0.95, "Strict 1:1 orthologs\n6,796 genes\n(Ensembl Compara)", fc="#e8f4e8", ec="#33a02c")
box(6.9, 8.7, 2.9, 0.95, "MAFFT MSA + PHYLIP\n6,871 alignments\n(quality-filtered)", fc="#e8f4e8", ec="#33a02c")
arrow(3.2, 9.18, 3.6, 9.18); arrow(6.5, 9.18, 6.9, 9.18)

# Row 2: selection scans
box(0.3, 7.0, 4.4, 1.05, "Coding selection — HyPhy BUSTED v3\n4,974 genes tested; 1,690 BH FDR<0.05\n(PAML branch-site cross-validation)", fc="#dbe9f6", ec="#2166ac")
box(5.2, 7.0, 4.4, 1.05, "Selection intensity — HyPhy RELAX\n1,660 informative genes\n604 intensified (K>1) / 117 relaxed (K<1)", fc="#dbe9f6", ec="#2166ac")
arrow(3.0, 8.7, 2.5, 8.05); arrow(8.0, 8.7, 7.4, 8.05)

# Row 3: regulatory evidence
box(0.3, 5.3, 4.4, 1.05, "Regulatory evidence\nHARs (3,171) / caMPRA active HARs (508)\nhCONDELs · phyloP/phastCons", fc="#fde3d8", ec="#d6604d")
box(5.2, 5.3, 4.4, 1.05, "Expression evidence\nGTEx v11 tau (30 tissues)\nBrainSpan tau (26 regions)", fc="#fde3d8", ec="#d6604d")
arrow(2.5, 7.0, 2.5, 6.35); arrow(7.4, 7.0, 7.4, 6.35)

# Row 4: scoring + classification
box(1.2, 3.6, 7.6, 1.05, "Composite scoring & classification (Phase 7G v7)\nGDS = 0.35·p + 0.30·LRT + 0.20·RELAX + 0.15·Selectome (CDS-residualized)\n"
    "RDS = 0.30·caMPRA + 0.25·tau + 0.25·nc-accel + 0.20·brain-tau   (threshold = 0.15)", fc="#f5e6f5", ec="#8e44ad")
arrow(2.5, 5.3, 3.5, 4.65); arrow(7.4, 5.3, 6.5, 4.65)

# Row 5: classification results
box(0.3, 1.9, 4.4, 1.05, "Classification (4,974-gene universe)\nGD 1,214 + GD-relaxed 52 | RD 293\nDual 11 | Neutral 3,404", fc="#f5e6f5", ec="#8e44ad")
box(5.2, 1.9, 4.4, 1.05, "Robustness & non-circular validation\n1000× weight perturbation (94.6% agreement)\nLOO-caMPRA / LOO-tau / LOO-brain", fc="#f5e6f5", ec="#8e44ad")
arrow(4.0, 3.6, 2.5, 2.95); arrow(6.0, 3.6, 7.4, 2.95)

# Row 6: downstream
box(1.2, 0.3, 7.6, 1.0, "Downstream analyses: GO/KEGG/Reactome/SynGO enrichment · conservation (phyloP/phastCons) ·\n"
    "GWAS enrichment · tissue & cell-type resolution (GTEx/HPA/BrainSpan) · developmental trajectories", fc="#fff7d6", ec="#e6a817")
arrow(2.5, 1.9, 4.0, 1.3); arrow(7.4, 1.9, 6.0, 1.3)

ax.set_title("Analysis pipeline: gene-driven vs regulation-driven positive selection in the human lineage",
             fontsize=11, fontweight="bold", pad=10)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/Fig_pipeline_flowchart.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  Saved Fig_pipeline_flowchart.png")

# ==================================================================
# Fig 2: BUSTED p-value distribution
# ==================================================================
print("\n--- Fig 2: BUSTED p-value distribution ---")
p = pd.to_numeric(v7["busted_p"], errors="coerce").dropna().values
q = pd.to_numeric(v7["bh_fdr"], errors="coerce").dropna().values
q05 = np.sort(q[q < 0.05])
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
ax.plot(np.arange(1, len(qc)+1), qc, color="#2166ac", lw=1.2, label="sorted q-values")
ax.axhline(0.05, color="#d6604d", ls="--", lw=1.3, label="FDR = 0.05")
ax.set_xlabel("Gene rank (sorted by q)")
ax.set_ylabel("BH FDR q-value")
ax.set_title(f"(b) FDR q-values: {(q<0.05).sum():,} / {len(q):,} significant")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/Fig_pvalue_distribution.png", dpi=300)
plt.close(fig)
print("  Saved Fig_pvalue_distribution.png")

# ==================================================================
# Fig 3: GDS vs RDS scatter
# ==================================================================
print("\n--- Fig 3: GDS vs RDS scatter ---")
gds = pd.to_numeric(v7["gds_v7"], errors="coerce")
rds = pd.to_numeric(v7["rds_v7"], errors="coerce")
cls = v7["classification_v7"].astype(str)

COLORS = {"gene-driven": "#2166ac", "gene-driven (relaxed)": "#92c5de",
          "regulation-driven": "#d6604d", "dual-driven": "#e08214",
          "neutral": "#c8c8c8"}
THR = 0.15

fig, ax = plt.subplots(figsize=(6.4, 5.6))
lims = [-0.02, 1.02]
xx = np.linspace(lims[0], lims[1], 100)
ax.fill_between(xx, xx - THR, xx + THR, color="#f0f0f0", alpha=0.8, zorder=0,
                label="neutral band (|GDS−RDS| ≤ 0.15)")
ax.plot(xx, xx + THR, color="#999999", lw=0.9, ls="--")
ax.plot(xx, xx - THR, color="#999999", lw=0.9, ls="--")
ax.plot(xx, xx, color="#666666", lw=0.7, ls=":")

for c in ["neutral", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "gene-driven"]:
    m = (cls == c) & gds.notna() & rds.notna()
    ax.scatter(gds[m], rds[m], s=9 if c != "neutral" else 5, alpha=0.65,
               c=COLORS[c], label=f"{c} (n={m.sum():,})", edgecolors="none", zorder=2)

ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel("GDS (gene-driven score, v7)")
ax.set_ylabel("RDS (regulation-driven score, v7)")
ax.set_title("GDS vs RDS classification landscape (4,974 genes)")
ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/Fig_GDS_RDS_scatter.png", dpi=300)
plt.close(fig)
print("  Saved Fig_GDS_RDS_scatter.png")

# ==================================================================
# Fig 4: v5 vs v7 classification comparison
# ==================================================================
print("\n--- Fig 4: v5 vs v7 comparison ---")
m = v7[["gene_id", "classification_v7"]].merge(
    v5[["gene_id", "classification_v5"]], on="gene_id", how="inner")
print(f"  merged genes: {len(m)}")

CLASSES = ["gene-driven", "gene-driven (relaxed)", "gene-driven (dual)", "dual-driven",
           "regulation-driven", "neutral"]
lab_v5 = ["gene-driven", "gene-driven (relaxed)", "gene-driven (dual)", "dual-driven",
          "regulation-driven", "neutral"]
lab_v7 = ["gene-driven", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "neutral"]

M = np.zeros((len(lab_v5), len(lab_v7)), dtype=int)
for i, a in enumerate(lab_v5):
    for j, b in enumerate(lab_v7):
        M[i, j] = ((m["classification_v5"] == a) & (m["classification_v7"] == b)).sum()

fig = plt.figure(figsize=(11, 4.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1])

# (a) confusion matrix
ax = fig.add_subplot(gs[0, 0])
Mnorm = M / M.sum(axis=1, keepdims=True) * 100
im = ax.imshow(Mnorm, cmap="Blues", vmin=0, vmax=100)
ax.set_xticks(range(len(lab_v7))); ax.set_xticklabels(lab_v7, rotation=35, ha="right", fontsize=7.5)
ax.set_yticks(range(len(lab_v5))); ax.set_yticklabels(lab_v5, fontsize=7.5)
for i in range(len(lab_v5)):
    for j in range(len(lab_v7)):
        if M[i, j] > 0:
            ax.text(j, i, f"{M[i, j]:,}", ha="center", va="center", fontsize=6.8,
                    color="white" if Mnorm[i, j] > 50 else "black")
ax.set_xlabel("Classification v7 (current)")
ax.set_ylabel("Classification v5 (pre-revision)")
ax.set_title("(a) Classification transition: v5 → v7")
fig.colorbar(im, ax=ax, shrink=0.75, label="% of v5 class")

# (b) class counts
ax = fig.add_subplot(gs[0, 1])
cats = ["gene-driven", "gene-driven\n(relaxed)", "dual-driven", "regulation-driven", "neutral"]
v5_counts = [1620, 0, 19 + 57, 125, 3153]  # v5 dual split: dual 19 + gene-driven(dual) 57
v7_counts = [1214, 52, 11, 293, 3404]
xpos = np.arange(len(cats)); w = 0.38
ax.bar(xpos - w/2, v5_counts, w, label="v5 (pre-revision)", color="#92c5de", edgecolor="white")
ax.bar(xpos + w/2, v7_counts, w, label="v7 (current)", color="#2166ac", edgecolor="white")
for x, (a, b) in enumerate(zip(v5_counts, v7_counts)):
    ax.text(x - w/2, a + 30, f"{a:,}", ha="center", fontsize=7)
    ax.text(x + w/2, b + 30, f"{b:,}", ha="center", fontsize=7)
ax.set_xticks(xpos); ax.set_xticklabels(cats, fontsize=8)
ax.set_ylabel("Number of genes")
ax.set_title("(b) Class composition: v5 vs v7")
ax.legend(fontsize=8)
ax.set_ylim(0, 3900)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/Fig_v5_v7_comparison.png", dpi=300)
plt.close(fig)
print("  Saved Fig_v5_v7_comparison.png")

print("\n" + "=" * 70)
print("All 4 figures complete")
print("=" * 70)
