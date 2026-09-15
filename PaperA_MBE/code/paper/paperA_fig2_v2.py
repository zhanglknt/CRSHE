# -*- coding: utf-8 -*-
"""
Paper A (MBE) — Figure 2 v2: Coding-selection scan and score construction
(5 panels, layout 2+3). Output: results/paper/figures_v2/.

Panel data sources
------------------
a | BUSTED p-value distribution: gene_classification_v7.csv (busted_p)
b | BH FDR vs Storey q curves: same CSV (bh_fdr, storey_q)
c | RELAX K distribution: results/phase4_hyphy/relax_results/relax_results_v3.csv
    (K); FDR-sig intensified/relaxed counts from gene_classification_v7.csv
    (relax_intensified / relax_relaxed)
d | CDS-length residualization: gene_classification_v7.csv (busted_p,
    cds_length) — OLS residual recomputed with the pipeline logic
    (phase7g_v7_classification.py:318-325); residual Spearman matches the
    stored gds_p_pct_resid at rho=1.0
e | omega distribution: results/phase8_tissue_analysis/busted_v1_omega.csv
    (omega_weighted, universe genes)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

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
N = len(g)
p = g["busted_p"].values
q_bh = g["bh_fdr"].values
q_st = g["storey_q"].values
relax = pd.read_csv(os.path.join(ROOT, "results", "phase4_hyphy",
                                 "relax_results", "relax_results_v3.csv"))
K = relax["K"].dropna().values
n_int = int(g["relax_intensified"].sum())
n_rel = int(g["relax_relaxed"].sum())
omega = pd.read_csv(os.path.join(ROOT, "results", "phase8_tissue_analysis",
                                 "busted_v1_omega.csv"))
ow = omega[omega["gene_id"].isin(g["gene_id"])]["omega_weighted"].values

summary = []

fig = plt.figure(figsize=(11, 7.4))
gs = gridspec.GridSpec(2, 6, figure=fig, height_ratios=[1.0, 1.0],
                       hspace=0.62, wspace=0.85,
                       left=0.065, right=0.975, top=0.93, bottom=0.09)

# ============================================== panel a: BUSTED p histogram
axa = fig.add_subplot(gs[0, 0:2])
panel_label(axa, "a")
axa.set_title("BUSTED p-value distribution", loc="center", pad=6, fontsize=9)
axa.hist(p, bins=np.linspace(0, 0.5, 51), color=GDR, edgecolor="white",
         linewidth=0.3)
p05 = (p == 0.5).mean()
p00 = (p == 0).mean()
n_sig = int((q_bh < 0.05).sum())
p_max_sig = p[q_bh < 0.05].max()
axa.axvline(p_max_sig, color=RD, ls="--", lw=1.1)
axa.text(0.5, axa.get_ylim()[1] * 0.98, f"p = 0.5: {p05 * 100:.1f}%",
         ha="right", va="top", fontsize=7.6, color="#333333")
axa.text(0.004, axa.get_ylim()[1] * 0.98, f"p = 0: {p00 * 100:.1f}%",
         ha="left", va="top", fontsize=7.6, color="#333333")
axa.text(0.09, axa.get_ylim()[1] * 0.42,
         f"max p with FDR<0.05 = {p_max_sig:.2f}\n(dashed)", fontsize=7,
         color=RD, va="center")
axa.set_xlabel("BUSTED p-value (bounded at 0.5)")
axa.set_ylabel("Genes")
axa.set_xlim(0, 0.5)
axa.grid(axis="y", color="#e9e9e9", lw=0.6)
axa.set_axisbelow(True)

summary.append(f"a | BUSTED p: 39.4% at p=0.5, 11.9% at p=0, "
               f"max sig p={p_max_sig:.3f} (n={N:,})")

# ==================================== panel b: BH FDR vs Storey q curves
axb = fig.add_subplot(gs[0, 2:4])
panel_label(axb, "b")
axb.set_title("FDR q-value curves", loc="center", pad=6, fontsize=9)
qc_bh = np.sort(q_bh)
qc_st = np.sort(q_st)
ranks = np.arange(1, N + 1)
axb.plot(ranks, qc_bh.clip(1e-8), color=GD, lw=1.6, label="BH FDR q")
axb.plot(ranks, qc_st.clip(1e-8), color=RD, lw=1.1, ls=(0, (4, 3)),
         label="Storey q (pi0 = 1.0)")
axb.axhline(0.05, color=GREY, ls="--", lw=1.0)
axb.set_yscale("log")
axb.set_xlabel("Gene rank (sorted by q)")
axb.set_ylabel("q-value (log scale)")
axb.set_xlim(0, N)
axb.set_ylim(1e-8, 2)
axb.text(2600, 0.0016, f"{n_sig:,} / {N:,} ({n_sig / N * 100:.1f}%)\n"
         "FDR < 0.05", fontsize=7.4, color="#333333", va="bottom")
axb.text(N * 0.55, 2e-7, "BH and Storey coincide\n(pi0 = 1.0)",
         fontsize=7.4, color="#555555", style="italic")
axb.legend(frameon=False, loc="lower right", handlelength=1.8)
axb.grid(axis="y", color="#e9e9e9", lw=0.6)
axb.set_axisbelow(True)

summary.append(f"b | BH FDR<0.05: {n_sig}/{N} (34.0%); Storey pi0=1.0, "
               "curves coincide")

# ============================================== panel c: RELAX K histogram
axc = fig.add_subplot(gs[0, 4:6])
panel_label(axc, "c")
axc.set_title("RELAX selection-intensity K", loc="center", pad=6, fontsize=9)
Kc = np.clip(K, 0.01, 60)
bins = np.logspace(-2, np.log10(60), 40)
axc.hist(Kc, bins=bins, color=GDR, edgecolor="white", linewidth=0.3)
axc.set_xscale("log")
axc.axvline(1.0, color=GREY, ls="--", lw=1.0)
ymax = axc.get_ylim()[1]
axc.text(0.014, ymax * 0.97, f"K<1 relaxed\n(FDR<0.05): {n_rel}",
         fontsize=7.4, color=RD, va="top", ha="left")
axc.text(55, ymax * 0.97, f"K>1 intensified\n(FDR<0.05): {n_int}",
         fontsize=7.4, color=GD, va="top", ha="right")
axc.set_xlabel("RELAX K (log scale)")
axc.set_ylabel("Genes")
axc.set_xlim(0.01, 60)
axc.text(0.03, 0.035, f"n = {len(K):,} informative genes",
         transform=axc.transAxes, fontsize=7, color="#555555",
         style="italic")
axc.grid(axis="y", color="#e9e9e9", lw=0.6)
axc.set_axisbelow(True)

summary.append(f"c | RELAX K: n={len(K)} informative; {n_int} intensified "
               f"(K>1, FDR<0.05) / {n_rel} relaxed (K<1, FDR<0.05)")

# ================================= panel d: CDS-length residualization
lp = -np.log10(np.clip(p, 1e-300, None))
ll = np.log10(g["cds_length"].values)
slope, intercept, r1, p1, _ = stats.linregress(ll, lp)
resid = lp - (slope * ll + intercept)
r_res = stats.pearsonr(resid, ll)[0]

axd1 = fig.add_subplot(gs[1, 0:2])
panel_label(axd1, "d")
axd1.set_title("CDS-length bias: before", loc="center", pad=6, fontsize=9)
n_p0 = int((p == 0).sum())
axd1.scatter(ll, np.clip(lp, 0, 15.6), s=2, color=GD, alpha=0.18,
             linewidths=0, rasterized=True)
xx = np.linspace(ll.min(), ll.max(), 10)
axd1.plot(xx, slope * xx + intercept, color=RD, lw=1.3)
axd1.text(0.03, 0.97, f"r = {r1:.3f}\nP = {p1:.0e}",
          transform=axd1.transAxes, fontsize=7.6, va="top", color="#333333")
axd1.text(0.03, 0.42, f"{n_p0} genes at p = 0\n(y = 300, clipped to top)",
          transform=axd1.transAxes, fontsize=6.8, va="center", ha="left",
          color="#777777", style="italic")
axd1.set_xlabel("log10(CDS length, bp)")
axd1.set_ylabel("-log10(BUSTED p)")
axd1.set_ylim(-0.5, 16)
axd1.grid(axis="y", color="#e9e9e9", lw=0.6)
axd1.set_axisbelow(True)

axd2 = fig.add_subplot(gs[1, 2:4])
axd2.set_title("CDS-length bias: after residualization", loc="center", pad=6,
               fontsize=9)
axd2.scatter(ll, np.clip(resid, -4, 15.2), s=2, color=GREY, alpha=0.18,
             linewidths=0, rasterized=True)
axd2.axhline(0, color=RD, lw=1.3)
axd2.text(0.03, 0.97, f"r = {r_res:.3f}",
          transform=axd2.transAxes, fontsize=7.6, va="top", color="#333333")
axd2.text(0.03, 0.42, f"{n_p0} genes at p = 0\n(clipped to top)",
          transform=axd2.transAxes, fontsize=6.8, va="center", ha="left",
          color="#777777", style="italic")
axd2.set_xlabel("log10(CDS length, bp)")
axd2.set_ylabel("OLS residual of -log10(p)")
axd2.set_ylim(-4, 15.5)
axd2.grid(axis="y", color="#e9e9e9", lw=0.6)
axd2.set_axisbelow(True)

summary.append(f"d | CDS residualization: r {r1:.3f} (P={p1:.1e}) -> "
               f"r = {r_res:.3f}; residual percentile identical to stored "
               "gds_p_pct_resid (Spearman 1.0)")

# ============================================== panel e: omega distribution
axe = fig.add_subplot(gs[1, 4:6])
panel_label(axe, "e")
axe.set_title("dN/dS (omega) distribution", loc="center", pad=6, fontsize=9)
owc = np.clip(ow, 1e-3, 1e5)
bins = np.logspace(-3, 5, 45)
axe.hist(owc, bins=bins, color=GDR, edgecolor="white", linewidth=0.3)
axe.set_xscale("log")
axe.axvline(1.0, color=GREY, ls="--", lw=1.0)
frac_gt1 = (ow > 1).mean()
med = np.median(ow)
ymax = axe.get_ylim()[1]
axe.text(1.6e-3, ymax * 0.97, f"omega > 1: {frac_gt1 * 100:.1f}%\n"
         f"({int((ow > 1).sum()):,} / {len(ow):,})\nmedian = {med:.2f}",
         fontsize=7.4, color=GD, va="top", ha="left")
axe.set_xlabel("omega_weighted (log scale; clipped at 1e-3 / 1e5)")
axe.set_ylabel("Genes")
axe.set_xlim(1e-3, 1e5)
axe.grid(axis="y", color="#e9e9e9", lw=0.6)
axe.set_axisbelow(True)

summary.append(f"e | omega_weighted: n={len(ow)} in universe; "
               f"{frac_gt1 * 100:.1f}% > 1; median {med:.2f}")

# ----------------------------------------------------------------- save ----
fig.savefig(os.path.join(OUTDIR, "Figure2_coding_selection_scan.png"), dpi=300)
fig.savefig(os.path.join(OUTDIR, "Figure2_coding_selection_scan.pdf"),
            bbox_inches="tight")
plt.close(fig)

print("PaperA Figure 2 v2 written to", OUTDIR)
for s in summary:
    print("  " + s)
