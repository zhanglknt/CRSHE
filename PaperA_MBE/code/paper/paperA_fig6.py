# -*- coding: utf-8 -*-
"""
Paper A (MBE) — Figure 6: Functional content of the classes
(6 panels, 2x3). Output: results/paper/figures_v2/Figure6_class_functional_content.

Data source: results/phase6_functional_annotation/{class}_{lib}.csv
(Enrichr reports; columns: Gene_set, Term, Overlap, P-value,
Adjusted P-value, Odds Ratio, ...). All statistics read dynamically.

NOTE on significance scale: Enrichr BH-adjusted P-values are highly
conservative for these gene sets (most libraries have 0 terms at FDR<0.05;
signal concentrates in the nominal layer). Dot plots therefore use
x = -log10(nominal P) with filled markers for FDR<0.05 and open markers
for nominal-only terms; panel (f) shows FDR / nominal counts per cell.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

if os.path.exists("/mnt/d"):
    ROOT = "/mnt/d/人类正选择基因项目"
else:
    ROOT = "D:/人类正选择基因项目"
D = os.path.join(ROOT, "results", "phase6_functional_annotation")
OUTDIR = os.path.join(ROOT, "results", "paper", "figures_v2")
os.makedirs(OUTDIR, exist_ok=True)

GD = "#2166ac"; GDR = "#92c5de"; RD = "#d6604d"; NHI = "#b2182b"
DUAL = "#7b3294"; ALL = "#969696"; NEUT = "#bdbdbd"; GREY = "#737373"

plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
})

CLS = ["gene_driven", "regulation_driven", "dual_driven", "all_candidates"]
LIBS = ["GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "SynGO"]
CLS_COLOR = {"gene_driven": GD, "regulation_driven": RD,
             "dual_driven": DUAL, "all_candidates": ALL}
CLS_SHORT = {"gene_driven": "GD", "regulation_driven": "RD",
             "dual_driven": "dual", "all_candidates": "all cand."}

# ------------------------------------------------------------- load all ----
data = {}
for c in CLS:
    for l in LIBS:
        f = os.path.join(D, f"{c}_{l}.csv")
        data[(c, l)] = pd.read_csv(f).sort_values("P-value")

sig_counts = {(c, l): (int((data[(c, l)]["P-value"] < 0.05).sum()),
                       int((data[(c, l)]["Adjusted P-value"] < 0.05).sum()))
              for c in CLS for l in LIBS}
summary = []


def panel_label(ax, letter):
    ax.set_title(f"({letter})", loc="left", fontweight="bold", fontsize=11,
                 pad=6)


def trunc(t, n=34):
    return t if len(t) <= n else t[:n - 1] + "…"


# ---------------------------------------------------------------- figure ---
fig = plt.figure(figsize=(11, 7.6))
gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1.15, 1.0],
                       hspace=0.52, wspace=1.3,
                       left=0.185, right=0.975, top=0.94, bottom=0.08)


def dot_panel(ax, lib, topn, letter, title):
    """GD vs RD top-N dot plot, x = -log10(nominal P)."""
    panel_label(ax, letter)
    ax.set_title(title, loc="center", pad=6, fontsize=9)
    rows = []  # (y, term, x, color, fdr_sig, cls, overlap)
    y = 2 * topn
    for c in ["gene_driven", "regulation_driven"]:
        df = data[(c, lib)].head(topn)
        for _, r in df.iterrows():
            x = -np.log10(r["P-value"])
            fdr = r["Adjusted P-value"] < 0.05
            rows.append((y, r["Term"], x, CLS_COLOR[c], fdr, c,
                         r["Overlap"]))
            y -= 1
        y -= 1  # gap between groups
    for yy, term, x, color, fdr, c, ov in rows:
        ax.scatter(x, yy, s=30, facecolor=color if fdr else "white",
                   edgecolor=color, linewidth=1.1, zorder=3)
    ax.set_yticks([r[0] for r in rows])
    ax.set_yticklabels([trunc(r[1]) for r in rows], fontsize=6.5)
    # colour the y-tick labels by class instead of in-axes group tags
    for tick, r in zip(ax.get_yticklabels(), rows):
        tick.set_color(r[3])
    ax.axvline(-np.log10(0.05), color=GREY, ls="--", lw=0.9, zorder=1)
    r_max = max(r[2] for r in rows)
    xmax = max(r_max * 1.15, 4.05)  # room for the legend note on the right
    ax.set_xlim(0, xmax)
    ax.set_ylim(-2.4, 2 * topn + 0.6)
    ax.set_xlabel("-log$_{10}$(nominal P)")
    ax.grid(axis="x", color="#e9e9e9", lw=0.6)
    ax.set_axisbelow(True)
    ax.text(-np.log10(0.05) + 0.15, -1.75,
            "filled: FDR<0.05 · open: nominal",
            fontsize=6.2, color=GREY, style="italic", va="center")
    ngd, fgd = sig_counts[("gene_driven", lib)]
    nrd, frd = sig_counts[("regulation_driven", lib)]
    return (f"{lib}: GD {ngd} nominal/{fgd} FDR; RD {nrd} nominal/{frd} FDR")


# (a) GO-BP top-8
axa = fig.add_subplot(gs[0, 0])
summary.append("a | " + dot_panel(axa, "GO_BP", 8, "a", "GO-BP (top 8)"))

# (b) KEGG top-6
axb = fig.add_subplot(gs[0, 1])
summary.append("b | " + dot_panel(axb, "KEGG", 6, "b", "KEGG (top 6)"))

# (c) Reactome top-6
axc = fig.add_subplot(gs[0, 2])
summary.append("c | " + dot_panel(axc, "Reactome", 6, "c", "Reactome (top 6)"))

# =============================================== panel d: SynGO bars
axd = fig.add_subplot(gs[1, 0])
panel_label(axd, "d")
axd.set_title("SynGO terms", loc="center", pad=6, fontsize=9)
nom = [sig_counts[(c, "SynGO")][0] for c in CLS]
fdr = [sig_counts[(c, "SynGO")][1] for c in CLS]
xd = np.arange(len(CLS))
cols = [CLS_COLOR[c] for c in CLS]
axd.bar(xd, nom, width=0.62, color=cols, alpha=0.38, edgecolor=cols,
        linewidth=0.8, label="nominal P<0.05")
axd.bar(xd, fdr, width=0.62, color=cols, label="FDR<0.05")
for xi, (n, f) in enumerate(zip(nom, fdr)):
    axd.text(xi, n + 0.12, f"{n} ({f})", ha="center", fontsize=7.6,
             color="#333333")
axd.set_xticks(xd)
axd.set_xticklabels([CLS_SHORT[c] for c in CLS], fontsize=7.8)
axd.set_ylabel("Significant SynGO terms")
axd.set_ylim(0, max(nom) * 1.3 + 1)
axd.yaxis.set_major_locator(plt.MaxNLocator(5, integer=True))
axd.grid(axis="y", color="#e9e9e9", lw=0.6)
axd.set_axisbelow(True)
axd.legend(frameon=False, loc="upper right", fontsize=6.8, handlelength=1.1)
summary.append(f"d | SynGO sig terms (nominal/FDR): GD {nom[0]}/{fdr[0]}, "
               f"RD {nom[1]}/{fdr[1]}, dual {nom[2]}/{fdr[2]}, "
               f"all {nom[3]}/{fdr[3]}")

# =============================================== panel e: GO-MF top-6
axe = fig.add_subplot(gs[1, 1])
summary.append("e | " + dot_panel(axe, "GO_MF", 6, "e", "GO-MF (top 6)"))

# =============================================== panel f: overview matrix
axf = fig.add_subplot(gs[1, 2])
panel_label(axf, "f")
axf.set_title("Libraries x classes", loc="center", pad=6, fontsize=9)
mat_fdr = np.array([[sig_counts[(c, l)][1] for c in CLS] for l in LIBS],
                   dtype=float)
mat_nom = np.array([[sig_counts[(c, l)][0] for c in CLS] for l in LIBS],
                   dtype=float)
im = axf.imshow(mat_fdr, cmap="Blues", aspect="auto", vmin=0,
                vmax=max(mat_fdr.max(), 1))
for i in range(len(LIBS)):
    for j in range(len(CLS)):
        fv, nv = int(mat_fdr[i, j]), int(mat_nom[i, j])
        txt = f"{fv} / {nv}"
        color = "white" if fv > mat_fdr.max() * 0.55 else "#333333"
        axf.text(j, i, txt, ha="center", va="center", fontsize=7.2,
                 color=color,
                 fontweight="bold" if fv > 0 else "normal")
axf.set_xticks(range(len(CLS)))
axf.set_xticklabels([CLS_SHORT[c] for c in CLS], fontsize=7.6)
axf.set_yticks(range(len(LIBS)))
axf.set_yticklabels(LIBS, fontsize=7.6)
axf.set_xlabel("cell: FDR<0.05 / nominal P<0.05 terms")
for spine in axf.spines.values():
    spine.set_visible(False)
cb = fig.colorbar(im, ax=axf, fraction=0.046, pad=0.03)
cb.set_label("FDR<0.05 terms", fontsize=7)
cb.ax.tick_params(labelsize=6.5)
summary.append("f | matrix FDR/nominal per library x class (see stdout)")

# ----------------------------------------------------------------- save ----
fig.savefig(os.path.join(OUTDIR, "Figure6_class_functional_content.png"),
            dpi=300)
fig.savefig(os.path.join(OUTDIR, "Figure6_class_functional_content.pdf"),
            bbox_inches="tight")
plt.close(fig)

print("PaperA Figure 6 written to", OUTDIR)
for s in summary:
    print("  " + s)
print("  full matrix (rows=libs, cols=GD/RD/dual/all), nominal/FDR:")
for l in LIBS:
    print("   ", l, {CLS_SHORT[c]: sig_counts[(c, l)] for c in CLS})
