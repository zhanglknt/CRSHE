# -*- coding: utf-8 -*-
"""Round-4 figure fixes (three-reviewer merged list):
Fig1: full redraw per R1's 7-point spec (reciprocal BLAST dual-track, 6,908->6,871->6,840->4,974,
      no PAML phantom, unclassified, 4-component LOO, no unreported downstream analyses)
Fig3: legend neutral -> unclassified; axis labels GDS/RDS (no internal version)
Fig4: (b) retitled 'RD enrichment under LOO variants' + 95% CI error bars
FigS1: de-versioned labels (previous/current framework, unclassified)
FigS4: add (b) k-means 100-seed recovery distribution panel
Outputs -> results/paper/figures_v9/
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

BASE = os.environ.get("HSD_BASE") or "D:/人类正选择基因项目"
OUT = f"{BASE}/results/paper/figures_v9"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8})

v7 = pd.read_csv(f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv", low_memory=False)
v5 = pd.read_csv(f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v5/gene_classification_v5.csv", low_memory=False)
loo = pd.read_csv(f"{BASE}/results/paper/revision_v7/loo_full_matrix.csv")

# ==================================================================
# Figure 1: pipeline flowchart (REDRAWN)
# ==================================================================
print("--- Figure 1: pipeline flowchart (redrawn) ---")
fig, ax = plt.subplots(figsize=(10, 7.2))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, text, fc="#dbe9f6", ec="#2166ac", fs=8.5, lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec=ec, lw=lw))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, wrap=True)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12, color="#444444", lw=1.0))

# Row 1: input & orthology (reciprocal BLAST, dual track)
box(0.3, 8.7, 2.7, 0.95, "10 primate genomes\n13,509–23,485 CDS/species\n(Ensembl)", fc="#e8f4e8", ec="#33a02c")
box(3.35, 8.7, 3.1, 0.95, "Reciprocal-BLAST orthology\nstrict 1:1 set (6,796)\npermissive shared set (6,908)", fc="#e8f4e8", ec="#33a02c")
box(6.8, 8.7, 2.9, 0.95, "MAFFT MSA (6,871)\nPAML-format QC (6,840)\n= 5,717 strict + 1,123 permissive", fc="#e8f4e8", ec="#33a02c")
arrow(3.0, 9.18, 3.35, 9.18); arrow(6.45, 9.18, 6.8, 9.18)

# Row 2: selection scans
box(0.3, 7.0, 4.4, 1.05, "Coding selection — HyPhy BUSTED\n4,974 genes tested; 1,690 BH FDR < 0.05\n(human branch as foreground)", fc="#dbe9f6", ec="#2166ac")
box(5.2, 7.0, 4.4, 1.05, "Selection intensity — HyPhy RELAX\n1,701 candidates → 1,660 valid\n604 intensified (K>1) / 117 relaxed (0<K<1)", fc="#dbe9f6", ec="#2166ac")
arrow(3.0, 8.7, 2.5, 8.05); arrow(8.0, 8.7, 7.4, 8.05)

# Row 3: regulatory & expression evidence
box(0.3, 5.3, 4.4, 1.05, "Regulatory evidence\nHARs (3,171) / caMPRA-active HARs (508)\nhCONDELs (583) · phyloP/phastCons", fc="#fde3d8", ec="#d6604d")
box(5.2, 5.3, 4.4, 1.05, "Expression evidence\nGTEx v11 tau (30 tissues)\nBrainSpan tau (26 brain regions)", fc="#fde3d8", ec="#d6604d")
arrow(2.5, 7.0, 2.5, 6.35); arrow(7.4, 7.0, 7.4, 6.35)

# Row 4: scoring & classification
box(1.2, 3.6, 7.6, 1.05, "Comparability-corrected scoring and classification\nGDS = 0.35·p + 0.30·LRT + 0.20·RELAX + 0.15·Selectome (CDS-length-residualized)\n"
    "RDS = 0.30·caMPRA + 0.25·tau + 0.25·nc-divergence + 0.20·brain-tau   (margin = 0.15)", fc="#f5e6f5", ec="#8e44ad")
arrow(2.5, 5.3, 3.5, 4.65); arrow(7.4, 5.3, 6.5, 4.65)

# Row 5: classification & validation
box(0.3, 1.9, 4.4, 1.05, "Classification (4,974-gene universe)\nGD 1,214 + GD (relaxed) 52 | RD 293\nDual 11 | Unclassified 3,404", fc="#f5e6f5", ec="#8e44ad")
box(5.2, 1.9, 4.4, 1.05, "Leave-one-out validation (4 RDS components)\nLOO-caMPRA / LOO-tau / LOO-nc / LOO-brain\nHAR & hCONDEL enrichment (Table 3)", fc="#f5e6f5", ec="#8e44ad")
arrow(4.0, 3.6, 2.5, 2.95); arrow(6.0, 3.6, 7.4, 2.95)

# Row 6: sensitivity & cross-reference (only analyses reported in the paper)
box(1.2, 0.3, 7.6, 1.0, "Sensitivity and cross-reference: threshold/weight/covariate robustness · k-means unsupervised control ·\n"
    "Selectome and Shao et al. cross-reference · MEME site-level evidence (GARD non-convergent)", fc="#fff7d6", ec="#e6a817")
arrow(2.5, 1.9, 4.0, 1.3); arrow(7.4, 1.9, 6.0, 1.3)

ax.set_title("Analysis pipeline: gene-driven versus regulation-driven selection in the human lineage",
             fontsize=11, fontweight="bold", pad=10)
fig.tight_layout()
fig.savefig(f"{OUT}/Figure1_analysis_pipeline.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Figure1_analysis_pipeline.pdf", bbox_inches="tight")
plt.close(fig)
print("  saved Figure1_analysis_pipeline.png")

# ==================================================================
# Figure 3: GDS vs RDS scatter (unclassified label, de-versioned axes)
# ==================================================================
print("--- Figure 3: GDS-RDS scatter ---")
gds = pd.to_numeric(v7["gds_v7"], errors="coerce")
rds = pd.to_numeric(v7["rds_v7"], errors="coerce")
cls = v7["classification_v7"].astype(str)
COLORS = {"gene-driven": "#2166ac", "gene-driven (relaxed)": "#92c5de",
          "regulation-driven": "#d6604d", "dual-driven": "#e08214", "neutral": "#c8c8c8"}
DISP = {"neutral": "unclassified"}
THR = 0.15
fig, ax = plt.subplots(figsize=(6.4, 5.6))
lims = [-0.02, 1.02]
xx = np.linspace(lims[0], lims[1], 100)
ax.fill_between(xx, xx - THR, xx + THR, color="#f0f0f0", alpha=0.8, zorder=0,
                label="unclassified band (|GDS−RDS| ≤ 0.15)")
ax.plot(xx, xx + THR, color="#999999", lw=0.9, ls="--")
ax.plot(xx, xx - THR, color="#999999", lw=0.9, ls="--")
ax.plot(xx, xx, color="#666666", lw=0.7, ls=":")
for c in ["neutral", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "gene-driven"]:
    m = (cls == c) & gds.notna() & rds.notna()
    ax.scatter(gds[m], rds[m], s=9 if c != "neutral" else 5, alpha=0.65,
               c=COLORS[c], label=f"{DISP.get(c, c)} (n={m.sum():,})", edgecolors="none", zorder=2)
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel("GDS (gene-driven score)")
ax.set_ylabel("RDS (regulation-driven score)")
ax.set_title("GDS versus RDS classification landscape (4,974 genes)")
ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
fig.tight_layout()
fig.savefig(f"{OUT}/Figure3_GDS_RDS_scatter.png", dpi=300)
fig.savefig(f"{OUT}/Figure3_GDS_RDS_scatter.pdf", bbox_inches="tight")
plt.close(fig)
print("  saved Figure3_GDS_RDS_scatter.png")

# ==================================================================
# Figure 4: LOO enrichment matrix (CI error bars, retitled)
# ==================================================================
print("--- Figure 4: LOO enrichment matrix ---")
variants = loo["variant"].tolist()
gd_tot = (loo["gene-driven"] + loo["gene-driven (relaxed)"]).values
rd = loo["regulation-driven"].values

def parse_ci(s):
    return [float(x) for x in s.strip("[]").split(",")]

har_or = loo["HAR_OR"].values; har_ci = np.array([parse_ci(s) for s in loo["HAR_CI"]])
hc_or = loo["hCONDEL_OR"].values; hc_ci = np.array([parse_ci(s) for s in loo["hCONDEL_CI"]])

fig = plt.figure(figsize=(11, 4.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15])
ax = fig.add_subplot(gs[0, 0])
xpos = np.arange(len(variants)); w = 0.38
ax.bar(xpos - w/2, gd_tot, w, label="GD (+relaxed)", color="#2166ac", edgecolor="white")
ax.bar(xpos + w/2, rd, w, label="RD", color="#d6604d", edgecolor="white")
for x, (a, b) in enumerate(zip(gd_tot, rd)):
    ax.text(x - w/2, a + 18, f"{a:,}", ha="center", fontsize=7)
    ax.text(x + w/2, b + 18, f"{b:,}", ha="center", fontsize=7)
ax.set_xticks(xpos); ax.set_xticklabels(variants, fontsize=8, rotation=15)
ax.set_ylabel("Genes"); ax.set_title("(a) Class sizes by LOO variant"); ax.legend(fontsize=8)

ax = fig.add_subplot(gs[0, 1])
xp = np.arange(len(variants))
ax.errorbar(xp - 0.12, har_or, yerr=[har_or - har_ci[:, 0], har_ci[:, 1] - har_or],
            fmt="o", color="#2166ac", ms=6, capsize=4, lw=1.2, label="HAR OR (95% CI)")
ax.errorbar(xp + 0.12, hc_or, yerr=[hc_or - hc_ci[:, 0], hc_ci[:, 1] - hc_or],
            fmt="s", color="#d6604d", ms=6, capsize=4, lw=1.2, label="hCONDEL OR (95% CI)")
for x, (o1, o2) in enumerate(zip(har_or, hc_or)):
    ax.text(x - 0.12, o1 + 0.45, f"{o1:.2f}", ha="center", fontsize=7, color="#2166ac")
    ax.text(x + 0.12, o2 - 0.75, f"{o2:.2f}", ha="center", fontsize=7, color="#d6604d")
ax.axhline(1.0, color="#888888", lw=0.8, ls="--")
ax.set_xticks(xp); ax.set_xticklabels(variants, fontsize=8, rotation=15)
ax.set_ylabel("Odds ratio (RD vs rest, Fisher exact)")
ax.set_title("(b) RD enrichment under LOO variants (95% CI)")
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig(f"{OUT}/Figure4_LOO_enrichment_matrix.png", dpi=300)
fig.savefig(f"{OUT}/Figure4_LOO_enrichment_matrix.pdf", bbox_inches="tight")
plt.close(fig)
print("  saved Figure4_LOO_enrichment_matrix.png")

# ==================================================================
# Figure S1: previous vs current framework comparison (de-versioned)
# ==================================================================
print("--- Figure S1: framework comparison ---")
m = v7[["gene_id", "classification_v7"]].merge(v5[["gene_id", "classification_v5"]], on="gene_id", how="inner")
lab_v5 = ["gene-driven", "gene-driven (relaxed)", "gene-driven (dual)", "dual-driven", "regulation-driven", "neutral"]
lab_v7 = ["gene-driven", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "neutral"]
disp_v5 = ["gene-driven", "gene-driven (relaxed)", "gene-driven (dual)", "dual-driven", "regulation-driven", "unclassified"]
disp_v7 = ["gene-driven", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "unclassified"]
M = np.zeros((len(lab_v5), len(lab_v7)), dtype=int)
for i, a in enumerate(lab_v5):
    for j, b in enumerate(lab_v7):
        M[i, j] = ((m["classification_v5"] == a) & (m["classification_v7"] == b)).sum()
fig = plt.figure(figsize=(11, 4.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1])
ax = fig.add_subplot(gs[0, 0])
Mnorm = M / M.sum(axis=1, keepdims=True) * 100
im = ax.imshow(Mnorm, cmap="Blues", vmin=0, vmax=100)
ax.set_xticks(range(len(lab_v7))); ax.set_xticklabels(disp_v7, rotation=35, ha="right", fontsize=7.5)
ax.set_yticks(range(len(lab_v5))); ax.set_yticklabels(disp_v5, fontsize=7.5)
for i in range(len(lab_v5)):
    for j in range(len(lab_v7)):
        if M[i, j] > 0:
            ax.text(j, i, f"{M[i, j]:,}", ha="center", va="center", fontsize=6.8,
                    color="white" if Mnorm[i, j] > 50 else "black")
ax.set_xlabel("Current framework"); ax.set_ylabel("Previous framework")
ax.set_title("(a) Classification transition: previous → current")
fig.colorbar(im, ax=ax, shrink=0.75, label="% of previous-framework class")
ax = fig.add_subplot(gs[0, 1])
cats = ["gene-driven", "gene-driven\n(relaxed)", "dual-driven", "regulation-driven", "unclassified"]
v5c = [1620, 0, 19 + 57, 125, 3153]
v7c = [1214, 52, 11, 293, 3404]
xpos = np.arange(len(cats)); w = 0.38
ax.bar(xpos - w/2, v5c, w, label="previous framework", color="#92c5de", edgecolor="white")
ax.bar(xpos + w/2, v7c, w, label="current framework", color="#2166ac", edgecolor="white")
for x, (a, b) in enumerate(zip(v5c, v7c)):
    ax.text(x - w/2, a + 30, f"{a:,}", ha="center", fontsize=7)
    ax.text(x + w/2, b + 30, f"{b:,}", ha="center", fontsize=7)
ax.set_xticks(xpos); ax.set_xticklabels(cats, fontsize=8)
ax.set_ylabel("Number of genes")
ax.set_title("(b) Class composition (previous dual pools 19 + 57)")
ax.legend(fontsize=8); ax.set_ylim(0, 3900)
fig.tight_layout()
fig.savefig(f"{OUT}/FigureS1_framework_comparison.png", dpi=300)
fig.savefig(f"{OUT}/FigureS1_framework_comparison.pdf", bbox_inches="tight")
plt.close(fig)
print("  saved FigureS1_framework_comparison.png")

# ==================================================================
# Figure S4: PCA + k-means recovery panel
# ==================================================================
print("--- Figure S4: PCA + k-means ---")
feats = ["gds_p_pct_resid", "gds_lrt_pct_resid", "gds_relax_pct", "gds_selectome",
         "rds_doan", "rds_tau", "rds_nc", "rds_brain"]
X = v7[feats].astype(float).values
Xz = (X - X.mean(0)) / X.std(0)
rd_mask = (v7["classification_v7"] == "regulation-driven").values
n_rd = int(rd_mask.sum())

pca = PCA(n_components=2, random_state=42)
PC = pca.fit_transform(Xz)

rec = []
for seed in range(100):
    km = KMeans(n_clusters=4, n_init=10, random_state=seed).fit(Xz)
    lab = km.labels_
    counts = np.array([(lab[rd_mask] == c).sum() for c in range(4)])
    rec.append(int(counts.max()) / n_rd)
rec = np.array(rec) * 100
print(f"  k-means recovery: mean {rec.mean():.1f}% sd {rec.std():.1f}% median {np.median(rec):.1f}% ge80 {(rec>=80).mean():.2f}")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axes[0]
for c in ["neutral", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "gene-driven"]:
    mm = (cls == c)
    ax.scatter(PC[mm, 0], PC[mm, 1], s=8 if c != "neutral" else 4, alpha=0.6,
               c=COLORS[c], label=DISP.get(c, c), edgecolors="none")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("(a) PCA of the eight component scores")
ax.legend(fontsize=7.5, markerscale=1.6)
ax = axes[1]
ax.hist(rec, bins=20, color="#2166ac", edgecolor="white", linewidth=0.4)
ax.axvline(rec.mean(), color="#d6604d", ls="--", lw=1.2,
           label=f"mean {rec.mean():.1f}% ± {rec.std():.1f}%")
ax.axvline(80, color="#888888", ls=":", lw=1.0, label="80% recovery")
ax.set_xlabel("Best-cluster recovery of RD genes (%)")
ax.set_ylabel("Seeds (of 100 × 10 restarts)")
ax.set_title("(b) k-means (k=4) recovery of the RD class")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/FigureS4_unsupervised_PCA_kmeans.png", dpi=300)
fig.savefig(f"{OUT}/FigureS4_unsupervised_PCA_kmeans.pdf", bbox_inches="tight")
plt.close(fig)
print("  saved FigureS4_unsupervised_PCA_kmeans.png")
print("ALL ROUND-4 FIGURES DONE")
