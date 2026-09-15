# -*- coding: utf-8 -*-
"""
Paper A (MBE) — Figure 5: Independent evidence and cross-study context (6 panels, 2x3)
=======================================================================================
Design system per results/paperB/figure_spec_v1.md, identical to paperA_fig{1..4}_v2.py
(Arial; GD #2166ac, GD-relax #92c5de, RD #d6604d, dual #e08214, neutral #bdbdbd;
bold 11pt panel labels; PNG 300dpi + vector PDF).
NOTE: figure_spec_fig56.md mentions "dual=purple"; all existing Paper A scripts
(fig1-4 v2, paperA_figures.py) use dual #e08214 — kept for cross-figure consistency.

Panel data sources (all statistics computed on the fly with scipy, nothing hardcoded):
  a | gene_classification_v7.csv — phyloP CDS (x) vs phyloP gene-body (y):
      neutral as grey hexbin density, GD/RD foreground scatter, median crosses
  b | same — phyloP_cds by class (5 classes), two-sided MWU p (GD vs RD)
  c | same — nc_conservation by class (non-coding component mechanistic basis),
      two-sided MWU p (GD vs RD)
  d | same + data/selectome/selectome_primate_positive_selection.tsv —
      has_selectome proportion by class; Fisher exact GD vs non-GD (greater)
  e | data/downloads/shao_2023/table_s17_psgs.csv — Shao 2023 PSG (82), ENSG
      version-stripped match to the v7 universe: nested bars 82 -> in-universe ->
      BUSTED-FDR-sig; Fisher depletion OR/p conditioned on analyzable genes
  f | same as (a) — n_hars distribution among genes with >=1 HAR (GD/RD/neutral
      ECDF) + Fisher OR for any-HAR enrichment (RD vs non-RD); n_hars is NOT an
      RDS component, so this test is non-circular

Output: results/paper/figures_v2/Figure5_independent_evidence.png (300 dpi) + .pdf
Stdout prints per-panel key statistics (p / OR / n) for manuscript citation.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

if os.path.exists("/mnt/d"):
    ROOT = "/mnt/d/人类正选择基因项目"
else:
    ROOT = "D:/人类正选择基因项目"
V7 = os.path.join(ROOT, "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv")
SHAO = os.path.join(ROOT, "data/downloads/shao_2023/shao2023_psg_82_genes.txt")
SEL = os.path.join(ROOT, "data/selectome/selectome_primate_positive_selection.tsv")
OUTDIR = os.path.join(ROOT, "results/paper/figures_v2")
os.makedirs(OUTDIR, exist_ok=True)

C_GD, C_GD_RELAX, C_RD, C_DUAL, C_NEUT = "#2166ac", "#92c5de", "#d6604d", "#e08214", "#bdbdbd"
CLS_ORDER = ["gene-driven", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "neutral"]
CLS_COLOR = {"gene-driven": C_GD, "gene-driven (relaxed)": C_GD_RELAX, "dual-driven": C_DUAL,
             "regulation-driven": C_RD, "neutral": C_NEUT}
CLS_SHORT = {"gene-driven": "GD", "gene-driven (relaxed)": "GD-rel.", "dual-driven": "dual",
             "regulation-driven": "RD", "neutral": "neutral"}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.7, "axes.edgecolor": "#333333",
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.labelsize": 8.5,
    "axes.titlesize": 9.5, "legend.fontsize": 7.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
FS_LABEL, FS_TITLE = 11, 9.5
REPORT = []


def style_ax(ax, ygrid=False, xgrid=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if ygrid:
        ax.grid(axis="y", color="#e5e5e5", linewidth=0.6, zorder=0)
    if xgrid:
        ax.grid(axis="x", color="#e5e5e5", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def panel_label(ax, letter, dx=-0.16, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=FS_LABEL,
            fontweight="bold", va="top", ha="left")


def fmt_p(p):
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"{m}×10$^{{{int(e)}}}$"
    return f"{p:.3f}"


# ---------------------------------------------------------------- load data
gc = pd.read_csv(V7)
cls = gc["classification_v7"].to_numpy()
gd_mask = gc["classification_v7"] == "gene-driven"
rd_mask = gc["classification_v7"] == "regulation-driven"

fig = plt.figure(figsize=(11, 8))
gs = fig.add_gridspec(2, 3, left=0.06, right=0.99, top=0.94, bottom=0.075,
                      wspace=0.34, hspace=0.52, width_ratios=[1.12, 1.06, 1.0])

# ================================================================ panel a
# phyloP CDS (x) vs phyloP gene-body (y): neutral density + GD/RD scatter
axa = fig.add_subplot(gs[0, 0])
x_all, y_all = gc["phyloP_cds"].to_numpy(), gc["phyloP_gene"].to_numpy()
neu = gc["classification_v7"] == "neutral"
axa.hexbin(x_all[neu], y_all[neu], gridsize=30, cmap="Greys", mincnt=1, alpha=0.45,
           linewidths=0.1, edgecolors="#dddddd", zorder=1)
rng = np.random.default_rng(7)
for m, col, lab in [(gd_mask, C_GD, "GD"), (rd_mask, C_RD, "RD")]:
    idx = np.where(m)[0]
    axa.scatter(x_all[idx], y_all[idx], s=7, alpha=0.55, c=col, edgecolors="white",
                linewidths=0.12, zorder=3, label=f"{lab} (n = {len(idx):,})")
# median crosses
for m, col in [(gd_mask, C_GD), (rd_mask, C_RD), (neu.to_numpy(), "#7f7f7f")]:
    mx, my = np.median(x_all[m]), np.median(y_all[m])
    axa.scatter([mx], [my], marker="+", s=160, linewidths=1.8, color=col, zorder=4)
axa.set_xlabel("phyloP, CDS (coding constraint)")
axa.set_ylabel("phyloP, gene body")
axa.set_title("Coding vs gene-body conservation", loc="left", fontsize=FS_TITLE)
style_ax(axa)
panel_label(axa, "a", dx=-0.075)
handles_a = [Line2D([0], [0], marker="s", color="none", markerfacecolor="#cccccc",
                    markeredgecolor="#aaaaaa", ms=7, label=f"neutral density (n = {neu.sum():,})")]
handles_a += axa.get_legend_handles_labels()[0]
lega = axa.legend(handles=handles_a, loc="upper left", frameon=True, fontsize=6.8,
                  handletextpad=0.4, labelspacing=0.3, borderaxespad=0.1)
lega.get_frame().set_facecolor("white")
lega.get_frame().set_alpha(0.85)
lega.get_frame().set_edgecolor("none")
lega.set_zorder(5)
med = {k: (np.median(x_all[m]), np.median(y_all[m])) for k, m in
       [("GD", gd_mask), ("RD", rd_mask), ("neutral", neu.to_numpy())]}
REPORT.append("a | phyloP medians (CDS, gene-body): GD=(%.2f, %.2f), RD=(%.2f, %.2f), "
              "neutral=(%.2f, %.2f)" % (*med["GD"], *med["RD"], *med["neutral"]))

# ================================================================ panel b
# phyloP_cds by class boxplot + MWU GD vs RD
axb = fig.add_subplot(gs[0, 1])
data_b = [gc.loc[gc["classification_v7"] == c, "phyloP_cds"].to_numpy() for c in CLS_ORDER]
bp = axb.boxplot(data_b, positions=range(5), widths=0.55, patch_artist=True,
                 showfliers=False, medianprops=dict(color="black", linewidth=0.9),
                 whiskerprops=dict(linewidth=0.7), capprops=dict(linewidth=0.7),
                 boxprops=dict(linewidth=0.7))
for patch, c in zip(bp["boxes"], CLS_ORDER):
    patch.set_facecolor(CLS_COLOR[c])
    patch.set_alpha(0.9)
for i, d in enumerate(data_b):
    vs = rng.choice(d, size=min(len(d), 100), replace=False)
    axb.scatter(np.full(len(vs), i) + rng.uniform(-0.12, 0.12, len(vs)), vs,
                s=2.0, color=CLS_COLOR[CLS_ORDER[i]], alpha=0.35, edgecolors="none", zorder=3)
u_stat, p_gd_rd = stats.mannwhitneyu(data_b[0], data_b[3], alternative="two-sided")
axb.set_xticks(range(5))
axb.set_xticklabels([CLS_SHORT[c] for c in CLS_ORDER], fontsize=7.4)
axb.set_ylabel("phyloP, CDS")
axb.set_title("CDS constraint by class", loc="left", fontsize=FS_TITLE)
axb.set_ylim(-1.5, 8.6)
style_ax(axb, ygrid=True)
panel_label(axb, "b", dx=-0.135)
axb.text(0.03, 0.975, f"RD vs GD: two-sided MWU\np = {fmt_p(p_gd_rd)} (RD more constrained)",
         transform=axb.transAxes, ha="left", va="top", fontsize=6.8, color="#555555")
REPORT.append("b | phyloP_cds MWU GD vs RD: U=%.0f, two-sided p=%.3e; medians GD=%.2f, RD=%.2f"
              % (u_stat, p_gd_rd, np.median(data_b[0]), np.median(data_b[3])))

# ================================================================ panel c
# non-coding divergence (-nc_conservation) by class + MWU GD vs RD
# (rds_nc correlates -1.0 with raw nc_conservation: lower raw = more diverged;
#  plot -nc_conservation so the RD class is highest, matching rds_nc direction)
axc = fig.add_subplot(gs[0, 2])
data_c = [(-gc.loc[gc["classification_v7"] == c, "nc_conservation"]).to_numpy()
          for c in CLS_ORDER]
vp = axc.violinplot(data_c, positions=range(5), widths=0.7, showmedians=True,
                    showextrema=False)
for body, c in zip(vp["bodies"], CLS_ORDER):
    body.set_facecolor(CLS_COLOR[c])
    body.set_alpha(0.85)
    body.set_edgecolor("white")
    body.set_linewidth(0.4)
vp["cmedians"].set_color("black")
vp["cmedians"].set_linewidth(0.9)
_, p_nc = stats.mannwhitneyu(data_c[0], data_c[3], alternative="two-sided")
axc.set_xticks(range(5))
axc.set_xticklabels([CLS_SHORT[c] for c in CLS_ORDER], fontsize=7.4)
axc.set_ylabel("Non-coding divergence ($-$nc_conservation)")
axc.set_title("Non-coding component (RDS basis)", loc="left", fontsize=FS_TITLE)
axc.set_ylim(-1.0, 7.8)
style_ax(axc, ygrid=True)
panel_label(axc, "c", dx=-0.135)
axc.text(0.03, 0.975, f"RD vs GD: two-sided MWU p = {fmt_p(p_nc)}\n(RDS component — mechanistic basis)",
         transform=axc.transAxes, ha="left", va="top", fontsize=6.8, color="#555555")
REPORT.append("c | -nc_conservation MWU GD vs RD: two-sided p=%.3e; medians GD=%.3f, RD=%.3f"
              % (p_nc, np.median(data_c[0]), np.median(data_c[3])))

# ================================================================ panel d
# Selectome cross-reference: has_selectome proportion by class + Fisher GD vs rest
axd = fig.add_subplot(gs[1, 0])
sel = pd.read_csv(SEL, sep="\t")
n_sel_file = len(sel)
cls_d = CLS_ORDER  # all five classes
props, counts_d, ns_d = [], [], []
for c in cls_d:
    m = gc["classification_v7"] == c
    k = int(gc.loc[m, "has_selectome"].sum())
    counts_d.append(k)
    ns_d.append(int(m.sum()))
    props.append(k / m.sum() * 100)
cols_d = [CLS_COLOR[c] for c in cls_d]
bars = axd.bar(range(5), props, width=0.62, color=cols_d, edgecolor="white",
               linewidth=0.5, zorder=3)
for i, (k, n, p) in enumerate(zip(counts_d, ns_d, props)):
    axd.text(i, p + 0.015, f"{k}/{n}", ha="center", fontsize=6.6, color="#444444")
# Fisher: GD (incl. relaxed) vs rest — primary; strict GD vs non-GD — secondary
gdi_mask = gc["classification_v7"].isin(["gene-driven", "gene-driven (relaxed)"])
k_gdi, n_gdi = int(gc.loc[gdi_mask, "has_selectome"].sum()), int(gdi_mask.sum())
k_rest = int(gc.loc[~gdi_mask, "has_selectome"].sum())
n_rest = int((~gdi_mask).sum())
or_d, p_d = stats.fisher_exact([[k_gdi, n_gdi - k_gdi], [k_rest, n_rest - k_rest]],
                               alternative="greater")
k_gd_s, n_gd_s = counts_d[0], ns_d[0]
or_d2, p_d2 = stats.fisher_exact(
    [[k_gd_s, n_gd_s - k_gd_s],
     [int(gc.loc[~gd_mask, "has_selectome"].sum()), int((~gd_mask).sum()) -
      int(gc.loc[~gd_mask, "has_selectome"].sum())]], alternative="greater")
axd.set_xticks(range(5))
axd.set_xticklabels([CLS_SHORT[c] for c in cls_d], fontsize=7.2)
axd.set_ylabel("Selectome-positive (%)")
axd.set_title("Selectome cross-reference", loc="left", fontsize=FS_TITLE)
axd.set_ylim(0, max(props) * 1.32)
style_ax(axd, ygrid=True)
panel_label(axd, "d", dx=-0.075)
n_sel_pos = int(gc["has_selectome"].sum())
axd.text(0.97, 0.975, f"{k_gdi} of {n_sel_pos} Selectome-positive genes are GD (incl. relaxed)\n"
         f"GD vs rest: Fisher OR = {or_d:.2f}, p = {fmt_p(p_d)}\n"
         f"(strict GD only: OR = {or_d2:.2f}, p = {fmt_p(p_d2)}; Selectome set: {n_sel_file})",
         transform=axd.transAxes, ha="right", va="top", fontsize=6.6, color="#555555")
REPORT.append("d | has_selectome: GD %d/%d, GD-rel %d/%d, dual %d/%d, RD %d/%d, neutral %d/%d; "
              "Fisher GD(incl.rel) vs rest OR=%.2f p=%.3e; strict GD OR=%.2f p=%.3e (Selectome file n=%d)"
              % (counts_d[0], ns_d[0], counts_d[1], ns_d[1], counts_d[2], ns_d[2],
                 counts_d[3], ns_d[3], counts_d[4], ns_d[4], or_d, p_d, or_d2, p_d2, n_sel_file))

# ================================================================ panel e
# Shao 2023 PSG overlap — ENSG-based matching (table_s17_psgs.csv Ensembl_ID,
# version-stripped) against the 4,974-gene v7 universe.
# NOTE: the often-quoted "OR = 0.07" used all 82 Shao genes as the denominator,
# counting 67 non-analyzable genes as BUSTED-non-significant. Conditioned on the
# 18 genes actually present in the universe, the depletion is OR ~ 0.55 (n.s.).
axe = fig.add_subplot(gs[1, 1])
shao_tab = pd.read_csv(os.path.join(ROOT, "data/downloads/shao_2023/table_s17_psgs.csv"))
strip_ensg = lambda e: str(e).split(".")[0]
shao_ensg = set(strip_ensg(e) for e in shao_tab["Ensembl_ID"])
gc["ensg_stripped"] = gc["gene_id"].map(strip_ensg)
shao_univ = gc["ensg_stripped"].isin(shao_ensg).to_numpy()
n_univ = len(gc)
sig_mask = gc["bh_fdr_sig"].to_numpy()
n_sig = int(sig_mask.sum())
n_shao_univ = int(shao_univ.sum())
k_sig_shao = int((sig_mask & shao_univ).sum())
tab_e = [[k_sig_shao, n_shao_univ - k_sig_shao],
         [n_sig - k_sig_shao, n_univ - n_shao_univ - (n_sig - k_sig_shao)]]
or_e, p_e_less = stats.fisher_exact(tab_e, alternative="less")
_, p_e_two = stats.fisher_exact(tab_e, alternative="two-sided")
ov_rows = gc.loc[sig_mask & shao_univ]
ov_symbols = sorted(ov_rows["gene_symbol"].dropna().tolist())
exp_overlap = n_shao_univ * n_sig / n_univ
# nested horizontal bars: 82 Shao PSG -> in universe -> BUSTED-sig
levels = [("Shao 2023 PSGs", len(shao_ensg), "#bdbdbd"),
          ("in 4,974-gene universe", n_shao_univ, "#8c8c8c"),
          ("BUSTED-FDR significant", k_sig_shao, C_RD)]
for i, (lab, val, col) in enumerate(levels):
    axe.barh(2 - i, val, height=0.62, color=col, edgecolor="white", linewidth=0.5,
             zorder=3)
    axe.text(val + 1.2, 2 - i, f"{val}", va="center", fontsize=7.6,
             fontweight="bold" if i == 2 else "normal",
             color=C_RD if i == 2 else "#555555")
axe.set_yticks([2, 1, 0])
axe.set_yticklabels([l[0] for l in levels], fontsize=7.4)
axe.set_xlabel("Genes")
axe.set_title("Shao 2023 PSG cross-study context", loc="left", fontsize=FS_TITLE)
axe.set_xlim(0, 96)
style_ax(axe, xgrid=True)
panel_label(axe, "e", dx=-0.30)
axe.text(0.97, 0.40, f"conditioned on {n_shao_univ} analyzable genes:\n"
         f"overlap {k_sig_shao} vs {exp_overlap:.1f} expected; Fisher depletion\n"
         f"OR = {or_e:.2f}, p = {fmt_p(p_e_less)} (n.s.); "
         f"{len(shao_ensg) - n_shao_univ}/{len(shao_ensg)} Shao PSGs lack a\n"
         f"1:1 alignment in our conservative universe",
         transform=axe.transAxes, ha="right", va="top", fontsize=6.4, color="#555555")
REPORT.append("e | Shao2023 (ENSG match): %d/%d in universe; overlap BUSTED-sig=%d "
              "(expected %.1f); Fisher less OR=%.3f p=%.3e (two-sided %.3e); overlap genes: %s "
              "(+1 without symbol); naive 82-denominator OR=0.072 NOT reproducible as valid test "
              "(treats 67 non-analyzable genes as non-sig)"
              % (n_shao_univ, len(shao_ensg), k_sig_shao, exp_overlap, or_e, p_e_less,
                 p_e_two, ", ".join(ov_symbols)))

# ================================================================ panel f
# n_hars ECDF (genes with >=1 HAR) + Fisher any-HAR OR RD vs non-RD
axf = fig.add_subplot(gs[1, 2])
max_h = int(gc["n_hars"].max())
for c, col, lw in [("neutral", C_NEUT, 1.2), ("gene-driven", C_GD, 1.5),
                   ("regulation-driven", C_RD, 1.5)]:
    v = np.sort(gc.loc[(gc["classification_v7"] == c) & (gc["n_hars"] > 0), "n_hars"].to_numpy())
    y = np.arange(1, len(v) + 1) / len(v)
    axf.step(np.concatenate([[v[0]], v]), np.concatenate([[0], y]), color=col, lw=lw,
             label=f"{CLS_SHORT[c]} (n = {len(v)})", zorder=2, where="post")
k_rd_h = int((rd_mask & (gc["n_hars"] > 0)).sum())
n_rd = int(rd_mask.sum())
k_non_h = int((~rd_mask & (gc["n_hars"] > 0)).sum())
n_non = int((~rd_mask).sum())
or_f, p_f = stats.fisher_exact([[k_rd_h, n_rd - k_rd_h], [k_non_h, n_non - k_non_h]],
                               alternative="greater")
axf.set_xlabel("HARs per gene (among genes with ≥1 HAR)")
axf.set_ylabel("ECDF")
axf.set_title("HAR load per gene (non-circular)", loc="left", fontsize=FS_TITLE)
axf.set_xticks(range(1, max_h + 1))
style_ax(axf, ygrid=True)
panel_label(axf, "f", dx=-0.135)
axf.legend(loc="lower right", frameon=False, fontsize=6.8, handletextpad=0.4,
           labelspacing=0.3, borderaxespad=0.1)
axf.text(0.97, 0.30, f"any-HAR enrichment, RD vs non-RD:\nFisher OR = {or_f:.2f}, "
         f"p = {fmt_p(p_f)}\n({k_rd_h}/{n_rd} RD vs {k_non_h}/{n_non} non-RD)",
         transform=axf.transAxes, ha="right", va="top", fontsize=6.6, color="#555555")
REPORT.append("f | n_hars>0: RD %d/%d, non-RD %d/%d; Fisher greater OR=%.2f p=%.3e; "
              "median n_hars (HAR+ genes): GD=%.0f, RD=%.0f, neutral=%.0f"
              % (k_rd_h, n_rd, k_non_h, n_non, or_f, p_f,
                 np.median(gc.loc[gd_mask & (gc.n_hars > 0), "n_hars"]),
                 np.median(gc.loc[rd_mask & (gc.n_hars > 0), "n_hars"]),
                 np.median(gc.loc[neu & (gc.n_hars > 0), "n_hars"])))

out_png = os.path.join(OUTDIR, "Figure5_independent_evidence.png")
out_pdf = os.path.join(OUTDIR, "Figure5_independent_evidence.pdf")
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
print("saved:", out_png)
print("saved:", out_pdf)
print("\n=== per-panel key statistics ===")
for line in REPORT:
    print(line)
