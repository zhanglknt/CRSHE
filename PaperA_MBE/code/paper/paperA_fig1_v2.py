# -*- coding: utf-8 -*-
"""
Paper A (MBE) — Figure 1 v2: Study design and classification framework
(5 panels, layout 2+3). Output: results/paper/figures_v2/ (new directory).

Panel data sources
------------------
a | compressed pipeline schematic (matplotlib patches; counts as in manuscript)
b | gene attrition funnel: 204,855 CDS -> 6,908/6,796 1:1 -> 6,871 MSA
    -> 6,840 QC -> 4,974 BUSTED universe (manuscript numbers)
c | class sizes: results/phase7_gene_vs_regulation/phase7g_classification_v7/
    gene_classification_v7.csv (classification_v7)
d | evidence availability: same CSV (busted_p / has_doan_campra / n_hars /
    tau / nc_conservation / has_selectome)
e | five classification variants:
    results/paper/revision_v7/loo_full_matrix.csv (full/LOO-caMPRA/LOO-tau/
    LOO-nc/LOO-brain)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

if os.path.exists("/mnt/d"):
    ROOT = "/mnt/d/人类正选择基因项目"
else:
    ROOT = "D:/人类正选择基因项目"
V7 = os.path.join(ROOT, "results", "phase7_gene_vs_regulation",
                  "phase7g_classification_v7")
OUTDIR = os.path.join(ROOT, "results", "paper", "figures_v2")
os.makedirs(OUTDIR, exist_ok=True)

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


def panel_label(ax, letter):
    ax.set_title(f"({letter})", loc="left", fontweight="bold", fontsize=11,
                 pad=6)


# ---------------------------------------------------------------- data -----
g = pd.read_csv(os.path.join(V7, "gene_classification_v7.csv"))
counts = g["classification_v7"].value_counts()
N = len(g)
lm = pd.read_csv(os.path.join(ROOT, "results", "paper", "revision_v7",
                              "loo_full_matrix.csv"))
summary = []

n_campra = int(g["has_doan_campra"].sum())
n_har = int((g["n_hars"] > 0).sum())
n_sel = int(g["has_selectome"].sum())

# ---------------------------------------------------------------- figure ---
fig = plt.figure(figsize=(11, 7.2))
gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1.05, 1.0],
                       hspace=0.55, wspace=0.52,
                       left=0.108, right=0.975, top=0.93, bottom=0.09)

# ======================================================== panel a: pipeline
axa = fig.add_subplot(gs[0, 0:2])
axa.set_xlim(0, 1); axa.set_ylim(0, 1); axa.axis("off")
panel_label(axa, "a")
axa.set_title("Analysis pipeline", loc="center", pad=6)


def box(ax, x, y, w, h, text, fc="white", ec="#4d4d4d", fs=7.4, lw=1.0,
        weight="normal"):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.006,rounding_size=0.012",
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color="#222222", weight=weight, zorder=3, linespacing=1.25)


def arrow(ax, x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=9, color="#4d4d4d", lw=1.1,
                                 zorder=1, shrinkA=1, shrinkB=1))


# row 1: input -> orthologs -> MSA -> QC
y1, h1 = 0.845, 0.135
w = 0.205
xs = [0.005, 0.265, 0.525, 0.785]
box(axa, xs[0], y1, w, h1, "10 primate\ngenomes (Ensembl)")
box(axa, xs[1], y1, w, h1, "1:1 orthologs\n6,908 permissive /\n6,796 strict")
box(axa, xs[2], y1, w, h1, "MAFFT MSA\n6,871 alignments")
box(axa, xs[3], y1, w, h1, "QC filtering\n6,840 alignments")
for i in range(3):
    arrow(axa, xs[i] + w + 0.006, y1 + h1 / 2, xs[i + 1] - 0.006, y1 + h1 / 2)

# row 2: coding scan (blue) | regulatory evidence (red)
y2, h2 = 0.565, 0.16
box(axa, 0.005, y2, 0.46, h2,
    "Coding selection scan — HyPhy BUSTED + RELAX\n"
    "4,974 tested · 1,690 BH FDR<0.05 ·\n"
    "1,660 RELAX-informative genes", ec=GD, lw=1.4)
box(axa, 0.53, y2, 0.46, h2,
    "Regulatory & expression evidence\n"
    "HARs · caMPRA · hCONDELs · phyloP/phastCons\n"
    "GTEx v11 tau · BrainSpan", ec=RD, lw=1.4)
arrow(axa, 0.25, y1 - 0.008, 0.24, y2 + h2 + 0.008)
arrow(axa, 0.75, y1 - 0.008, 0.75, y2 + h2 + 0.008)

# row 3: scoring -> classification -> validation
y3, h3 = 0.21, 0.19
box(axa, 0.005, y3, 0.30, h3,
    "GDS / RDS composite scores\n(CDS-length residualized,\n"
    "1000x weight perturbation)", fc="#eef4fa", ec=GD)
box(axa, 0.35, y3, 0.30, h3,
    "5-class classification\nGD 1,214 · GD-relaxed 52\n"
    "RD 293 · dual 11 · neutral 3,404", fc="#f5f5f5")
box(axa, 0.695, y3, 0.30, h3,
    "Non-circular validation\nLOO variants · hCONDELs\n"
    "tissue / cell resolution", fc="#fdeee9", ec=RD)
arrow(axa, 0.24, y2 - 0.008, 0.15, y3 + h3 + 0.008)
arrow(axa, 0.75, y2 - 0.008, 0.85, y3 + h3 + 0.008)
arrow(axa, 0.31, y3 + h3 / 2, 0.35 - 0.005, y3 + h3 / 2)
arrow(axa, 0.65, y3 + h3 / 2, 0.695 - 0.005, y3 + h3 / 2)

summary.append("a | pipeline: 10 primates -> 1:1 orthologs (6,908/6,796) -> "
               "MAFFT 6,871 -> QC 6,840 -> BUSTED 4,974 (1,690 FDR sig) + "
               "RELAX 1,660 -> GDS/RDS -> 5-class -> validation")

# ================================================== panel b: attrition funnel
axb = fig.add_subplot(gs[0, 2])
panel_label(axb, "b")
axb.set_title("Gene attrition", loc="center", pad=6)

stages = [
    ("Primate CDS (all)", 204855),
    ("1:1 orthologs (permissive)", 6908),
    ("1:1 orthologs (strict)", 6796),
    ("MAFFT alignments", 6871),
    ("QC-passed alignments", 6840),
    ("BUSTED universe", 4974),
]
yb = np.arange(len(stages))[::-1]
logw = [np.log10(v) for _, v in stages]
axb.barh(yb, logw, height=0.62, color="#92c5de", edgecolor="none")
axb.barh(yb[-1], logw[-1], height=0.62, color=GD, edgecolor="none")
for y, (name, v) in zip(yb, stages):
    axb.text(np.log10(v) + 0.04, y, f"{v:,}",
             va="center", ha="left", fontsize=7.4, color="#333333")
axb.set_yticks(yb)
axb.set_yticklabels([s for s, _ in stages], fontsize=7.3)
axb.set_xlim(3.4, 5.75)
axb.set_xticks([4, 5])
axb.set_xticklabels(["10$^4$", "10$^5$"])
axb.set_xlabel("Genes (log$_{10}$ scale)")
axb.grid(axis="x", color="#e9e9e9", lw=0.6)
axb.set_axisbelow(True)
axb.text(0.97, 0.04, "2.4% of initial CDS\nenter the universe",
         transform=axb.transAxes, fontsize=7, color=GD, style="italic",
         ha="right", va="bottom")

summary.append("b | attrition: 204,855 -> 6,908 -> 6,796 -> 6,871 -> 6,840 "
               "-> 4,974 (2.4%)")

# ==================================================== panel c: class sizes
axc = fig.add_subplot(gs[1, 0])
panel_label(axc, "c")
axc.set_title("Classification (v7, n = 4,974)", loc="center", pad=6)

order_c = ["neutral", "gene-driven", "regulation-driven",
           "gene-driven (relaxed)", "dual-driven"]
labels_c = {"gene-driven": "Gene-driven (GD)",
            "gene-driven (relaxed)": "GD (relaxed)",
            "regulation-driven": "Regulation-driven (RD)",
            "dual-driven": "Dual-driven", "neutral": "Unclassified"}
colors_c = {"gene-driven": GD, "gene-driven (relaxed)": GDR,
            "regulation-driven": RD, "dual-driven": DUAL, "neutral": NEUT}
vals = [counts[c] for c in order_c]
yc = np.arange(len(order_c))[::-1]
axc.barh(yc, vals, height=0.62, color=[colors_c[c] for c in order_c],
         edgecolor="none")
for y, v in zip(yc, vals):
    txt = f"{v:,} ({v / N * 100:.1f}%)"
    if v > max(vals) * 0.5:
        # long bar: label inside the bar, right-aligned
        axc.text(v - 45, y, txt, va="center", ha="right", fontsize=7.4,
                 color="white", fontweight="bold")
    else:
        axc.text(v + 55, y, txt, va="center", ha="left", fontsize=7.4,
                 color="#333333")
axc.set_yticks(yc)
axc.set_yticklabels([labels_c[c] for c in order_c], fontsize=7.8)
axc.set_xlim(0, max(vals) * 1.42)
axc.set_xlabel("Number of genes")
axc.xaxis.set_major_locator(plt.MaxNLocator(4))
axc.grid(axis="x", color="#e9e9e9", lw=0.6)
axc.set_axisbelow(True)

summary.append("c | classes: GD 1,214 (24.4%) / GD-relaxed 52 (1.0%) / "
               "RD 293 (5.9%) / dual 11 (0.2%) / neutral 3,404 (68.4%)")

# =========================================== panel d: evidence availability
axd = fig.add_subplot(gs[1, 1])
panel_label(axd, "d")
axd.set_title("Evidence availability", loc="center", pad=6)

ev = [
    ("BUSTED p (coding selection)", N, GD),
    ("RELAX K (intensity)", 1660, GD),
    ("Tissue tau (GTEx)", int(g["tau"].notna().sum()), RD),
    ("Non-coding conservation", int(g["nc_conservation"].notna().sum()), RD),
    ("HAR proximity", n_har, RD),
    ("caMPRA activity", n_campra, RD),
    ("Selectome", n_sel, GD),
]
pct = [v / N * 100 for _, v, _ in ev]
yd = np.arange(len(ev))[::-1]
axd.barh(yd, pct, height=0.62, color=[c for _, _, c in ev], edgecolor="none")
for y, v, p in zip(yd, [v for _, v, _ in ev], pct):
    txt = f"{v:,} ({p:.1f}%)" if p >= 0.5 else f"{v} ({p:.1f}%)"
    axd.text(min(p + 2, 88), y, txt, va="center", ha="left", fontsize=7.3,
             color="#333333")
axd.set_yticks(yd)
axd.set_yticklabels([n for n, _, _ in ev], fontsize=7.5)
axd.set_xlim(0, 132)
axd.set_xticks([0, 25, 50, 75, 100])
axd.set_xlabel("Genes with computable evidence (%)")
axd.grid(axis="x", color="#e9e9e9", lw=0.6)
axd.set_axisbelow(True)
axd.text(0.98, 0.03, "evidence-availability\nasymmetry",
         transform=axd.transAxes, ha="right", va="bottom", fontsize=7.4,
         color="#555555", style="italic")

summary.append(f"d | evidence availability: BUSTED 100% ({N}) / RELAX 33.4% "
               f"(1,660) / tau 100% / nc-cons 100% / HAR {n_har} "
               f"({n_har / N * 100:.1f}%) / caMPRA {n_campra} "
               f"({n_campra / N * 100:.1f}%) / Selectome {n_sel} "
               f"({n_sel / N * 100:.1f}%)")

# ============================================== panel e: 5 LOO variants
axe = fig.add_subplot(gs[1, 2])
panel_label(axe, "e")
axe.set_title("Classification variants", loc="center", pad=6)

var_rows = lm.set_index("variant")
order_e = ["full", "LOO-caMPRA", "LOO-tau", "LOO-nc", "LOO-brain"]
gd_e = [var_rows.loc[v, "gene-driven"] for v in order_e]
rd_e = [var_rows.loc[v, "regulation-driven"] for v in order_e]
xe = np.arange(len(order_e))
bw = 0.36
axe.bar(xe - bw / 2 - 0.01, gd_e, width=bw, color=GD, label="GD")
axe.bar(xe + bw / 2 + 0.01, rd_e, width=bw, color=RD, label="RD")
for xi, v in zip(xe, gd_e):
    axe.text(xi - bw / 2 - 0.01, v + 22, f"{v:,}", ha="center", fontsize=7,
             color=GD)
for xi, v in zip(xe, rd_e):
    axe.text(xi + bw / 2 + 0.01, v + 22, f"{v:,}", ha="center", fontsize=7,
             color=RD)
axe.set_xticks(xe)
axe.set_xticklabels(["full\n(v7)", "LOO-\ncaMPRA", "LOO-\ntau", "LOO-\nnc",
                    "LOO-\nbrain"], fontsize=7.5)
axe.set_ylabel("Genes classified")
axe.set_ylim(0, 1560)
axe.yaxis.set_major_locator(plt.MultipleLocator(400))
axe.grid(axis="y", color="#e9e9e9", lw=0.6)
axe.set_axisbelow(True)
axe.legend(frameon=False, loc="upper right", fontsize=7,
           handlelength=1.1, borderaxespad=0.2, labelspacing=0.3)

summary.append("e | 5 variants (loo_full_matrix.csv): full 1,214/293; "
               "LOO-caMPRA 1,020/800; LOO-tau 1,339/148; LOO-nc 1,369/110; "
               "LOO-brain 1,099/613")

# ----------------------------------------------------------------- save ----
fig.savefig(os.path.join(OUTDIR, "Figure1_design_classification.png"), dpi=300)
fig.savefig(os.path.join(OUTDIR, "Figure1_design_classification.pdf"),
            bbox_inches="tight")
plt.close(fig)

print("PaperA Figure 1 v2 written to", OUTDIR)
for s in summary:
    print("  " + s)
