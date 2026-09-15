# -*- coding: utf-8 -*-
"""
Paper B — Figure 1: Study design and gene classification (5 panels, layout 2+3)
Target journal: Genome Research. Design system per results/paperB/figure_spec_v1.md.

Panel data sources
------------------
a | Analysis schematic (drawn with matplotlib patches; counts from spec/summary JSONs)
b | Class sizes: results/phase8_tissue_analysis/data/phase8_master_table.csv (classification_v7)
c | GDS-RDS scatter: results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv
d | Data resources graphic: results/phase8_tissue_analysis/*.csv shape metadata
e | LOO variants: phase8_master_table.csv (classification_v7 / _loo_brain / _loo_tau)

Output: results/paperB/figures/Figure1_study_design.png (300 dpi) + .pdf
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.gridspec as gridspec

# ---------------------------------------------------------------- paths ----
ROOT = os.environ.get("HSD_BASE") or (
    "/mnt/d/人类正选择基因项目" if os.path.exists("/mnt/d") else "D:/人类正选择基因项目")
P8 = os.path.join(ROOT, "results", "phase8_tissue_analysis")
P7G = os.path.join(ROOT, "results", "phase7_gene_vs_regulation",
                   "phase7g_classification_v7")
OUTDIR = os.path.join(ROOT, "results", "paperB", "figures")
os.makedirs(OUTDIR, exist_ok=True)

# ------------------------------------------------------------- design ------
GD = "#2166ac"; GDR = "#92c5de"; RD = "#d6604d"; NHI = "#b2182b"
DUAL = "#e08214"; NEUT = "#bdbdbd"; GREY = "#737373"

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

CLASS_ORDER = ["gene-driven", "gene-driven (relaxed)", "regulation-driven",
               "dual-driven", "neutral"]
CLASS_COLORS = {"gene-driven": GD, "gene-driven (relaxed)": GDR,
                "regulation-driven": RD, "dual-driven": DUAL, "neutral": NEUT}
CLASS_LABELS = {"gene-driven": "Gene-driven (GD)",
                "gene-driven (relaxed)": "GD (relaxed)",
                "regulation-driven": "Regulation-driven (RD)",
                "dual-driven": "Dual-driven",
                "neutral": "Unclassified"}


def panel_label(ax, letter):
    ax.set_title(letter, loc="left", fontweight="bold", fontsize=11, pad=6)


# ---------------------------------------------------------------- data -----
master = pd.read_csv(os.path.join(P8, "data", "phase8_master_table.csv"))
genes = pd.read_csv(os.path.join(P7G, "gene_classification_v7.csv"))

counts_v7 = master["classification_v7"].value_counts()
counts_lob = master["classification_loo_brain"].value_counts()
counts_lot = master["classification_loo_tau"].value_counts()
N = len(master)
assert N == 4974

summary = []  # per-panel stdout summary

# ---------------------------------------------------------------- figure ---
fig = plt.figure(figsize=(11, 7.2))
gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1.0, 1.05],
                       hspace=0.55, wspace=0.42,
                       left=0.06, right=0.975, top=0.93, bottom=0.09)

# ======================================================= panel a: schematic
axa = fig.add_subplot(gs[0, 0:2])
axa.set_xlim(0, 1); axa.set_ylim(0, 1); axa.axis("off")
panel_label(axa, "a")
axa.set_title("Analysis schematic", loc="center", pad=6)

BOX_EC = "#4d4d4d"


def box(ax, x, y, w, h, text, fc="white", ec=BOX_EC, fs=7.6, lw=1.0,
        weight="normal", tc="black"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012",
                       fc=fc, ec=ec, lw=lw, mutation_aspect=1.0, zorder=2)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, color=tc, weight=weight, zorder=3, linespacing=1.25)
    return p


def arrow(ax, x1, y1, x2, y2, color="#4d4d4d", lw=1.2, style="-|>", ms=9,
          connectionstyle="arc3,rad=0"):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=ms,
                        color=color, lw=lw, zorder=1,
                        connectionstyle=connectionstyle, shrinkA=1, shrinkB=1)
    ax.add_patch(a)


# row 1: inputs -> alignments -> tests -> universe
y1, h = 0.845, 0.135
w = 0.205
xs = [0.005, 0.265, 0.525, 0.785]
box(axa, xs[0], y1, w, h, "10 primate\ngenomes")
box(axa, xs[1], y1, w, h, "6,840 ortholog\ncodon alignments\n(MAFFT)")
box(axa, xs[2], y1, w, h, "Selection tests\nBUSTED + RELAX\n(HyPhy)")
box(axa, xs[3], y1, w, h, "4,974-gene\nanalysis universe", weight="bold")
for i in range(3):
    arrow(axa, xs[i] + w + 0.006, y1 + h / 2, xs[i + 1] - 0.006, y1 + h / 2)

# row 2: scoring (GD blue / RD red)
y2, h2 = 0.575, 0.155
box(axa, 0.17, y2, 0.26, h2,
    "Gene-driven score (GDS)\ncoding selection evidence\n"
    "(BUSTED FDR, RELAX, Selectome)", ec=GD, lw=1.4)
box(axa, 0.57, y2, 0.26, h2,
    "Regulation-driven score (RDS)\ncaMPRA activity, tissue tau,\n"
    "nc-conservation, brain specificity", ec=RD, lw=1.4)
arrow(axa, xs[3] + w / 2 - 0.09, y1 - 0.008, 0.30, y2 + h2 + 0.008)
arrow(axa, xs[3] + w / 2 + 0.00, y1 - 0.008, 0.70, y2 + h2 + 0.008)

# row 3: classification
y3, h3 = 0.335, 0.13
box(axa, 0.27, y3, 0.46, h3,
    "GDS-RDS margin classification (|Δ| ≥ 0.15):  5 classes\n"
    "GD 1,214   GD-relaxed 52   RD 293   dual 11   unclassified 3,404",
    fc="#f5f5f5", weight="normal")
arrow(axa, 0.30, y2 - 0.008, 0.44, y3 + h3 + 0.008)
arrow(axa, 0.70, y2 - 0.008, 0.56, y3 + h3 + 0.008)

# row 4: three downstream layers
y4, h4 = 0.02, 0.20
lw4 = 0.30
box(axa, 0.01, y4, lw4, h4,
    "Layer 1 — bulk tissue\nGTEx v11, 30 tissues\nGD/RD enrichment, tau",
    fc="#eaf1f8", ec=GD)
box(axa, 0.35, y4, lw4, h4,
    "Layer 2 — single cell\nHPA 154 cell types +\nbrain 34 single-nuclei types",
    fc="#fdeee9", ec=RD)
box(axa, 0.69, y4, lw4, h4,
    "Layer 3 — development &\nindependent validation\n"
    "BrainSpan, hCONDELs, caMPRA", fc="#f5f5f5")
for xc in (0.16, 0.50, 0.84):
    arrow(axa, 0.50, y3 - 0.008, xc, y4 + h4 + 0.008,
          connectionstyle="arc3,rad=0")

summary.append(
    "a | schematic: 10 primates -> 6,840 alignments -> BUSTED/RELAX -> "
    "4,974-gene universe -> GDS/RDS -> 5-class -> 3 analysis layers")

# ===================================================== panel b: class sizes
axb = fig.add_subplot(gs[0, 2])
panel_label(axb, "b")
axb.set_title("Gene classes (n = 4,974)", loc="center", pad=6)

order_b = ["neutral", "gene-driven", "regulation-driven",
           "gene-driven (relaxed)", "dual-driven"]
vals = [counts_v7[c] for c in order_b]
labels = [CLASS_LABELS[c] for c in order_b]
colors = [CLASS_COLORS[c] for c in order_b]
ypos = np.arange(len(order_b))[::-1]
axb.barh(ypos, vals, height=0.62, color=colors, edgecolor="none")
for y, v in zip(ypos, vals):
    axb.text(v + 55, y, f"{v:,}  ({v / N * 100:.1f}%)",
             va="center", ha="left", fontsize=7.6, color="#333333")
axb.set_yticks(ypos)
axb.set_yticklabels(labels, fontsize=8)
axb.set_xlim(0, max(vals) * 1.42)
axb.set_xlabel("Number of genes")
axb.xaxis.set_major_locator(plt.MaxNLocator(4))
axb.grid(axis="x", color="#e6e6e6", lw=0.6)
axb.set_axisbelow(True)

summary.append("b | class sizes (v7): " + "; ".join(
    f"{CLASS_LABELS[c]} {counts_v7[c]:,} ({counts_v7[c] / N * 100:.1f}%)"
    for c in CLASS_ORDER))

# ================================================== panel c: GDS-RDS scatter
axc = fig.add_subplot(gs[1, 0])
panel_label(axc, "c")
axc.set_title("GDS vs RDS per gene", loc="center", pad=6)

plot_order = ["neutral", "gene-driven (relaxed)", "dual-driven",
              "regulation-driven", "gene-driven"]
sizes = {"neutral": 3, "gene-driven (relaxed)": 7, "dual-driven": 16,
         "regulation-driven": 7, "gene-driven": 5}
alphas = {"neutral": 0.25, "gene-driven (relaxed)": 0.8, "dual-driven": 0.95,
          "regulation-driven": 0.7, "gene-driven": 0.55}
for c in plot_order:
    sub = genes[genes["classification_v7"] == c]
    axc.scatter(sub["gds_v7"], sub["rds_v7"], s=sizes[c],
                c=CLASS_COLORS[c], alpha=alphas[c], linewidths=0,
                label=f"{CLASS_LABELS[c]} ({len(sub):,})", rasterized=True)
xx = np.linspace(0, 1.0, 10)
axc.plot(xx, xx, color="#bbbbbb", lw=0.8, ls="-", zorder=1)
axc.plot(xx, xx - 0.15, color=GREY, lw=0.9, ls="--", zorder=1)
axc.plot(xx, xx + 0.15, color=GREY, lw=0.9, ls="--", zorder=1)
axc.text(0.86, 0.66, "Δ = +0.15", fontsize=7, color=GREY, rotation=38,
         ha="center", va="center")
axc.text(0.36, 0.56, "Δ = −0.15", fontsize=7, color=GREY, rotation=38,
         ha="center", va="center")
axc.set_xlim(0, 0.98); axc.set_ylim(0, 0.95)
axc.set_xlabel("Gene-driven score (GDS, v7)")
axc.set_ylabel("Regulation-driven score (RDS, v7)")
leg = axc.legend(loc="upper right", frameon=True, framealpha=1.0,
                 edgecolor="#dddddd", handletextpad=0.15,
                 borderaxespad=0.1, labelspacing=0.28, markerscale=1.8)
leg.set_zorder(10)

summary.append("c | scatter n=4,974; class medians Δ=GDS-RDS: "
               "GD +0.43, RD -0.21, neutral +0.11 (band ±0.15)")

# ================================================ panel d: data resources
axd = fig.add_subplot(gs[1, 1])
axd.set_xlim(0, 1); axd.set_ylim(0, 1); axd.axis("off")
panel_label(axd, "d")
axd.set_title("Evidence layers and data resources", loc="center", pad=6)

resources = [
    ("GTEx v11", "30 bulk tissues; gene TPM", "bulk tissue", GD),
    ("HPA single-cell", "154 cell types, whole body", "single-cell", NHI),
    ("HPA brain single-nuclei", "34 brain cell types", "single-cell", NHI),
    ("BrainSpan", "26 regions (18 testable) x 31 dev. stages", "development", DUAL),
    ("hCONDELs", "583 human-specific deletions", "independent valid.", GREY),
    ("caMPRA HARs", "508 active / 3,171 total HARs", "independent valid.", GREY),
]
ytop = 0.845
dy = 0.148
axd.text(0.02, ytop + 0.068, "Resource", fontsize=7.8, weight="bold",
         color="#333333")
axd.text(0.02, ytop + 0.068, "", fontsize=7.8)
axd.text(0.995, ytop + 0.068, "Role", fontsize=7.8, weight="bold",
         color="#333333", ha="right")
for i, (name, scale, role, rc) in enumerate(resources):
    y = ytop - i * dy
    axd.plot([0.0, 1.0], [y + 0.058, y + 0.058], color="#e0e0e0", lw=0.7)
    axd.text(0.02, y + 0.022, name, fontsize=7.9, va="center", weight="bold",
             color="#222222")
    axd.text(0.995, y + 0.022, role, fontsize=7.0, va="center", ha="right",
             color="white",
             bbox=dict(boxstyle="round,pad=0.28", fc=rc, ec="none"))
    axd.text(0.06, y - 0.033, scale, fontsize=7.3, va="center",
             color="#555555")

summary.append("d | 6 data resources listed (GTEx30 / HPA154 / sn34 / "
               "BrainSpan 26 regions (18 testable) x 31 stages / hCONDEL 583 / "
               "caMPRA 508 of 3,171)")

# ============================================ panel e: LOO circularity ctrl
axe = fig.add_subplot(gs[1, 2])
panel_label(axe, "e")
axe.set_title("Circularity control: 5 classification variants",
              loc="center", pad=6)

# counts for all five variants (loo_extended_summary.json; v7/LOO-brain/
# LOO-tau cross-checked against phase8_master_table.csv value_counts)
loo_ext = json.load(open(os.path.join(
    P8, "loo_extended_summary.json")))["classification_counts"]
assert loo_ext["v7"]["gd"] == counts_v7["gene-driven"]
assert loo_ext["loo_brain"]["gd"] == counts_lob["gene-driven"]
assert loo_ext["loo_tau"]["gd"] == counts_lot["gene-driven"]

variants5 = [
    ("v7\n(full)", loo_ext["v7"]["gd"], loo_ext["v7"]["rd"]),
    ("LOO-\nbrain", loo_ext["loo_brain"]["gd"], loo_ext["loo_brain"]["rd"]),
    ("LOO-\ntau", loo_ext["loo_tau"]["gd"], loo_ext["loo_tau"]["rd"]),
    ("LOO-\ncaMPRA", loo_ext["loo_campra"]["gd"], loo_ext["loo_campra"]["rd"]),
    ("double-\nLOO", loo_ext["loo_double"]["gd"], loo_ext["loo_double"]["rd"]),
]
xv = np.arange(len(variants5))
bw = 0.36
gd_vals = [v[1] for v in variants5]
rd_vals = [v[2] for v in variants5]
# v7 = circular reference -> desaturated; LOO variants = non-circular, full
axe.bar(xv[0] - bw / 2 - 0.01, gd_vals[0], width=bw, color=GD, alpha=0.45,
        edgecolor=GD, lw=0.6)
axe.bar(xv[0] + bw / 2 + 0.01, rd_vals[0], width=bw, color=RD, alpha=0.45,
        edgecolor=RD, lw=0.6)
axe.bar(xv[1:] - bw / 2 - 0.01, gd_vals[1:], width=bw, color=GD,
        label="GD")
axe.bar(xv[1:] + bw / 2 + 0.01, rd_vals[1:], width=bw, color=RD,
        label="RD")
for xi, v in zip(xv, gd_vals):
    axe.text(xi - bw / 2 - 0.01, v + 24, f"{v:,}", ha="center", fontsize=7,
             color=GD)
for xi, v in zip(xv, rd_vals):
    axe.text(xi + bw / 2 + 0.01, v + 24, f"{v:,}", ha="center", fontsize=7,
             color=RD)
axe.set_xticks(xv)
axe.set_xticklabels([v[0] for v in variants5], fontsize=7.5)
axe.set_ylabel("Genes classified")
axe.set_ylim(0, 1650)
axe.yaxis.set_major_locator(plt.MultipleLocator(400))
axe.grid(axis="y", color="#e9e9e9", lw=0.6)
axe.set_axisbelow(True)
axe.legend(frameon=False, loc="upper right", fontsize=7,
           handlelength=1.1, borderaxespad=0.2, labelspacing=0.3)
axe.text(0.02, 0.985,
         "each LOO variant removes one circular evidence\n"
         "component (brain expression / tissue tau / caMPRA);\n"
         "double-LOO removes both. v7 = circular reference\n"
         "(desaturated).",
         transform=axe.transAxes, ha="left", va="top", fontsize=6.8,
         color="#444444", style="italic")

summary.append("e | 5 variants: v7 GD 1,214/RD 293; LOO-brain 1,099/613; "
               "LOO-tau 1,339/148; LOO-caMPRA 1,020/800; double-LOO 1,207/403")

# ----------------------------------------------------------------- save ----
fig.savefig(os.path.join(OUTDIR, "Figure1_study_design.png"), dpi=300)
fig.savefig(os.path.join(OUTDIR, "Figure1_study_design.pdf"),
            bbox_inches="tight")
plt.close(fig)

print("Figure 1 written to", OUTDIR)
for s in summary:
    print("  " + s)
