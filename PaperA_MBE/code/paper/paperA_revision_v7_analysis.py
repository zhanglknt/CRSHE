#!/usr/bin/env python3
"""
Paper A — Revision v7: Address simulated-review P0/P1 concerns
==============================================================
P0-D  : Complete leave-one-out matrix (caMPRA, tau, nc, brain) + HAR/hCONDEL
        enrichment re-estimated under each LOO variant
P1-D  : HAR-proximity vs caMPRA component correlation (nesting check)
P1-E  : Unsupervised control: PCA + k-means on 8 components; bootstrap class
        posteriors (from 1000-perturbation stability)
P0-B  : Evidence-availability quantification (coding vs regulatory coverage)
P0-A  : RELAX selection-bias quantification (RD/GD coverage numbers)
R3-M7 : 95% CIs for all key enrichment ORs

Output: results/paper/revision_v7/  (JSON + CSV + 2 figures)

Run: python.exe scripts/paperA_revision_v7_analysis.py
"""
import os
import json
import numpy as np
import pandas as pd
from collections import Counter
from scipy import stats
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.environ.get("HSD_BASE") or (
    "/mnt/d/人类正选择基因项目" if os.path.exists("/mnt/d") else "D:/人类正选择基因项目")
OUT = f"{BASE}/results/paper/revision_v7"
os.makedirs(OUT, exist_ok=True)

print("=" * 70)
print("Paper A — Revision v7 analyses (review P0/P1 fixes)")
print("=" * 70)

v7 = pd.read_csv(f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv",
                 low_memory=False)
stab = pd.read_csv(f"{BASE}/results/paper/sensitivity/per_gene_stability.csv")
n = len(v7)
print(f"Loaded {n} genes")

# Components
X_gds = v7[["gds_p_pct_resid", "gds_lrt_pct_resid", "gds_relax_pct", "gds_selectome"]].apply(
    pd.to_numeric, errors="coerce").fillna(0).values
rds_doan = pd.to_numeric(v7["rds_doan"], errors="coerce").fillna(0).values
rds_tau = pd.to_numeric(v7["rds_tau"], errors="coerce").fillna(0).values
rds_nc = pd.to_numeric(v7["rds_nc"], errors="coerce").fillna(0).values
rds_brain = pd.to_numeric(v7["rds_brain"], errors="coerce").fillna(0).values
sig = v7["bh_fdr_sig"].astype(str).str.strip().str.lower().eq("true").values
relaxed = v7["relax_relaxed"].astype(str).str.strip().str.lower().eq("true").values
har_arr = (pd.to_numeric(v7["n_hars"], errors="coerce").fillna(0) > 0).values
n_hars = pd.to_numeric(v7["n_hars"], errors="coerce").fillna(0).values
campra_arr = rds_doan.astype(bool)
baseline_cls = v7["classification_v7"].astype(str).values

GDS_W = np.array([0.35, 0.30, 0.20, 0.15])
gds = X_gds @ GDS_W

# hCONDEL gene set (from Phase 8 mapping — computed independently of RDS)
hc = pd.read_csv(f"{BASE}/results/phase8_tissue_analysis/hcondel_gene_mapping.csv")
hc_genes = set(hc["gene_id"]) if "gene_id" in hc.columns else set(hc.iloc[:, 0])
hcondel_arr = v7["gene_id"].isin(hc_genes).values

def classify(gds_v, rds_v, threshold=0.15):
    out = np.empty(len(gds_v), dtype=object)
    for i in range(len(gds_v)):
        g, r, s = gds_v[i], rds_v[i], sig[i]
        if g > r + threshold and s:
            c = "gene-driven"
        elif r > g + threshold and not s:
            c = "regulation-driven"
        elif g > r + threshold and r > 0.4 and s:
            c = "dual-driven"
        elif r > g + threshold and s and g > 0.3:
            c = "dual-driven"
        else:
            c = "neutral"
        if relaxed[i] and c == "gene-driven":
            c = "gene-driven (relaxed)"
        out[i] = c
    return out

def counts(c):
    cc = Counter(c)
    return {k: int(cc.get(k, 0)) for k in
            ["gene-driven", "gene-driven (relaxed)", "regulation-driven", "dual-driven", "neutral"]}

def fisher_ci(idx, feat):
    """Fisher exact one-sided + 95% CI for OR."""
    a = int(feat[idx].sum()); b = int(len(idx) - a)
    bg_a = int(feat.sum()) - a; bg_b = n - len(idx) - bg_a
    table = [[a, b], [bg_a, bg_b]]
    orv, p = stats.fisher_exact(table, alternative="greater")
    if a + b > 0 and bg_a + bg_b > 0 and min(a, b, bg_a, bg_b) >= 0 and a * b * bg_a * bg_b > 0:
        try:
            ci = scipy_odds_ratio(np.array(table), kind="sample").confidence_interval(0.95)
            lo, hi = float(np.nan_to_num(ci.low, nan=0.0)), float(np.nan_to_num(ci.high, nan=np.inf))
        except Exception:
            lo, hi = np.nan, np.nan
    else:
        lo, hi = np.nan, np.nan
    return float(orv), float(p), lo, hi, a, len(idx)

# Baseline
rds_full = 0.30 * rds_doan + 0.25 * rds_tau + 0.25 * rds_nc + 0.20 * rds_brain
cls_base = classify(gds, rds_full)
assert np.mean(cls_base == baseline_cls) == 1.0, "baseline reproduction failed"

# ------------------------------------------------------------------
# P0-D: complete LOO matrix
# ------------------------------------------------------------------
print("\n--- P0-D: complete leave-one-out matrix ---")
LOO = {
    "full":          dict(w=None, cls=cls_base),
    "LOO-caMPRA":    dict(w=dict(doan=None, tau=0.35, nc=0.35, brain=0.30)),
    "LOO-tau":       dict(w=dict(doan=0.40, tau=None, nc=0.33, brain=0.27)),
    "LOO-nc":        dict(w=dict(doan=0.40, tau=0.33, nc=None, brain=0.27)),
    "LOO-brain":     dict(w=dict(doan=0.375, tau=0.3125, nc=0.3125, brain=None)),
}
loo_rows = []
for name, spec in LOO.items():
    if spec["w"] is None:
        rds_v = rds_full
        cls = cls_base
    else:
        w = spec["w"]
        rds_v = (0.0 * rds_doan if w["doan"] is None else w["doan"] * rds_doan) + \
                (0.0 * rds_tau if w["tau"] is None else w["tau"] * rds_tau) + \
                (0.0 * rds_nc if w["nc"] is None else w["nc"] * rds_nc) + \
                (0.0 * rds_brain if w["brain"] is None else w["brain"] * rds_brain)
        cls = classify(gds, rds_v)
    rd_idx = np.where(cls == "regulation-driven")[0]
    gd_idx = np.where((cls == "gene-driven") | (cls == "gene-driven (relaxed)"))[0]
    or_h, p_h, lo_h, hi_h, k_h, n_h = fisher_ci(rd_idx, har_arr)
    or_c, p_c, lo_c, hi_c, k_c, n_c = fisher_ci(rd_idx, hcondel_arr)
    agree = float(np.mean(cls == cls_base))
    row = dict(variant=name, agreement=agree, **counts(cls),
               HAR_OR=or_h, HAR_p=p_h, HAR_CI=f"[{lo_h:.2f}, {hi_h:.2f}]" if lo_h == lo_h else "NA",
               HAR_k_over_n=f"{k_h}/{n_h}",
               hCONDEL_OR=or_c, hCONDEL_p=p_c,
               hCONDEL_CI=f"[{lo_c:.2f}, {hi_c:.2f}]" if lo_c == lo_c else "NA",
               hCONDEL_k_over_n=f"{k_c}/{n_c}")
    loo_rows.append(row)
    print(f"  {name:12s}: GD={row['gene-driven']}, RD={row['regulation-driven']}, "
          f"HAR OR={or_h:.2f} p={p_h:.1e}, hCONDEL OR={or_c:.2f} p={p_c:.1e}, agree={100*agree:.1f}%")
pd.DataFrame(loo_rows).to_csv(f"{OUT}/loo_full_matrix.csv", index=False)

# ------------------------------------------------------------------
# P1-D: HAR-caMPRA nesting correlation
# ------------------------------------------------------------------
print("\n--- P1-D: HAR proximity vs caMPRA component correlation ---")
phi = np.corrcoef(campra_arr.astype(float), har_arr.astype(float))[0, 1]
chi2_p = stats.chi2_contingency(np.histogram2d(campra_arr.astype(int), har_arr.astype(int), bins=2)[0])[1]
rho, rho_p = stats.spearmanr(campra_arr.astype(float), n_hars)
print(f"  phi(binary caMPRA, binary HAR)={phi:.3f}; Spearman(caMPRA, n_HARs)={rho:.3f} (p={rho_p:.1e})")
print(f"  caMPRA+ genes with HAR: {int((campra_arr & har_arr).sum())}/{int(campra_arr.sum())} "
      f"({100*(campra_arr & har_arr).sum()/max(campra_arr.sum(),1):.1f}%)")

# ------------------------------------------------------------------
# P1-E: unsupervised PCA + k-means control
# ------------------------------------------------------------------
print("\n--- P1-E: unsupervised control (PCA + k-means) ---")
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
X8 = np.column_stack([X_gds, rds_doan, rds_tau, rds_nc, rds_brain])
X8z = StandardScaler().fit_transform(X8)
pca = PCA(n_components=8).fit(X8z)
evr = pca.explained_variance_ratio_
print(f"  PCA explained variance: PC1={evr[0]:.1%}, PC2={evr[1]:.1%}, PC3={evr[2]:.1%}")
Z = PCA(n_components=2).fit_transform(X8z)
km4 = KMeans(n_clusters=4, n_init=20, random_state=42).fit_predict(X8z)
# cross-tab v7 class vs k-means cluster
ct = pd.crosstab(baseline_cls, km4)
ct.to_csv(f"{OUT}/pca_kmeans_crosstab.csv")
# which cluster best matches RD?
rd_best = ct.loc["regulation-driven"].idxmax() if "regulation-driven" in ct.index else None
purity = ct.max(axis=0).sum() / n if False else None
print(f"  k-means(4) cluster sizes: {dict(Counter(km4))}")
print(f"  best cluster for RD: {rd_best} "
      f"(contains {ct.loc['regulation-driven', rd_best]}/{(baseline_cls=='regulation-driven').sum()} RD genes)" if rd_best is not None else "")
# GD/RD separability check in PC space
gd_m = np.isin(baseline_cls, ["gene-driven", "gene-driven (relaxed)"])
rd_m = baseline_cls == "regulation-driven"
u_pcs, p_pcs = stats.mannwhitneyu(Z[rd_m, 0], Z[gd_m, 0])
u_pcs2, p_pcs2 = stats.mannwhitneyu(Z[rd_m, 1], Z[gd_m, 1])
print(f"  RD vs GD in PC1: p={p_pcs:.1e}; PC2: p={p_pcs2:.1e}")

# Bootstrap posteriors (per-gene stability as posterior probability of baseline class)
print("\n--- P1-E: bootstrap class posteriors ---")
stability = stab["stability"].values
print(f"  P(class=v7 class | perturbation): median={np.median(stability):.3f}; "
      f"genes with P>=0.8: {100*np.mean(stability>=0.8):.1f}%")

# ------------------------------------------------------------------
# P0-B: evidence availability quantification
# ------------------------------------------------------------------
print("\n--- P0-B: evidence availability asymmetry ---")
nonsig = ~sig
cov = {
    "coding_evidence": {"available": int(n), "pct": 100.0},
    "caMPRA_active_HAR": {"available": int(campra_arr.sum()), "pct": 100 * campra_arr.sum() / n},
    "any_HAR_pm50kb": {"available": int(har_arr.sum()), "pct": 100 * har_arr.sum() / n},
    "hCONDEL_mapped": {"available": int(hcondel_arr.sum()), "pct": 100 * hcondel_arr.sum() / n},
    "RELAX_attempted": {"available": int((v7['has_relax_active'].astype(str).str.lower().eq('true')).sum()),
                        "pct": 100 * (v7['has_relax_active'].astype(str).str.lower().eq('true')).sum() / n},
}
# RD yield conditioned on regulatory evidence availability
rd_with_campra = int(((baseline_cls == "regulation-driven") & campra_arr & nonsig).sum())
nonsig_with_campra = int((campra_arr & nonsig).sum())
nonsig_without = int((~campra_arr & nonsig).sum())
# how many nonsig genes would pass RDS threshold if caMPRA=1 were granted?
rds_hypo = rds_full + 0.30 * (1 - rds_doan)  # grant caMPRA evidence to all
cls_hypo = classify(gds, rds_hypo)
rd_hypo = int((cls_hypo == "regulation-driven").sum())
print(f"  Coverage: {cov}")
print(f"  Among BUSTED-non-significant genes WITH caMPRA evidence: "
      f"{rd_with_campra}/{nonsig_with_campra} = {100*rd_with_campra/max(nonsig_with_campra,1):.1f}% RD")
print(f"  Counterfactual (grant caMPRA=1 to all genes): RD {int((baseline_cls=='regulation-driven').sum())} -> {rd_hypo}")

# ------------------------------------------------------------------
# P0-A: RELAX selection bias numbers
# ------------------------------------------------------------------
print("\n--- P0-A: RELAX coverage by class ---")
relax_active = v7["has_relax_active"].astype(str).str.strip().str.lower().eq("true").values
for c in ["gene-driven", "gene-driven (relaxed)", "regulation-driven", "dual-driven", "neutral"]:
    m = baseline_cls == c
    if m.sum():
        print(f"  {c:22s}: {int(relax_active[m].sum())}/{int(m.sum())} with active RELAX "
              f"({100*relax_active[m].sum()/m.sum():.1f}%)")

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
print("\n--- Figures ---")
# Fig R1: LOO matrix heatmap of RD/HAR/hCONDEL across variants
fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
variants = [r["variant"] for r in loo_rows]
rd_counts = [r["regulation-driven"] for r in loo_rows]
gd_counts = [r["gene-driven"] + r["gene-driven (relaxed)"] for r in loo_rows]
xpos = np.arange(len(variants)); w = 0.38
axes[0].bar(xpos - w/2, gd_counts, w, color="#2166ac", label="GD (+relaxed)")
axes[0].bar(xpos + w/2, rd_counts, w, color="#d6604d", label="RD")
for x, (a, b) in enumerate(zip(gd_counts, rd_counts)):
    axes[0].text(x - w/2, a + 15, f"{a}", ha="center", fontsize=7)
    axes[0].text(x + w/2, b + 15, f"{b}", ha="center", fontsize=7)
axes[0].set_xticks(xpos); axes[0].set_xticklabels(variants, fontsize=8, rotation=15)
axes[0].set_ylabel("Genes"); axes[0].legend(fontsize=8)
axes[0].set_title("(a) Class counts under leave-one-out variants")

ors_h = [r["HAR_OR"] for r in loo_rows]
ors_c = [r["hCONDEL_OR"] for r in loo_rows]
axes[1].plot(variants, ors_h, "o-", color="#e08214", label="HAR enrichment (RD)")
axes[1].plot(variants, ors_c, "s-", color="#54278f", label="hCONDEL enrichment (RD)")
axes[1].axhline(1.0, color="grey", ls=":", lw=1)
axes[1].set_ylabel("Odds ratio"); axes[1].legend(fontsize=8)
axes[1].set_xticks(xpos); axes[1].set_xticklabels(variants, fontsize=8, rotation=15)
axes[1].set_title("(b) Non-circular RD enrichment under LOO variants")
fig.tight_layout()
fig.savefig(f"{OUT}/figR_loo_matrix.png", dpi=300)
plt.close(fig)

# Fig R2: PCA with classes
fig, ax = plt.subplots(figsize=(6.2, 5.2))
COLORS = {"gene-driven": "#2166ac", "gene-driven (relaxed)": "#92c5de",
          "regulation-driven": "#d6604d", "dual-driven": "#e08214", "neutral": "#c8c8c8"}
for c in ["neutral", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "gene-driven"]:
    m = baseline_cls == c
    ax.scatter(Z[m, 0], Z[m, 1], s=5, alpha=0.5, c=COLORS[c], label=f"{c} (n={m.sum()})", edgecolors="none")
ax.set_xlabel(f"PC1 ({evr[0]:.1%} variance)")
ax.set_ylabel(f"PC2 ({evr[1]:.1%} variance)")
ax.set_title("Unsupervised PCA of 8 classification components")
ax.legend(fontsize=7, loc="best")
fig.tight_layout()
fig.savefig(f"{OUT}/figR_pca.png", dpi=300)
plt.close(fig)
print("  Saved figR_loo_matrix.png, figR_pca.png")

# ------------------------------------------------------------------
# Summary JSON
# ------------------------------------------------------------------
summary = {
    "loo_matrix": loo_rows,
    "har_campra_correlation": {"phi_binary": float(phi), "spearman_nHARs": float(rho),
                               "spearman_p": float(rho_p),
                               "campra_genes_with_HAR_pct": float(100 * (campra_arr & har_arr).sum() / max(campra_arr.sum(), 1))},
    "pca": {"explained_variance": [float(x) for x in evr[:3]],
            "rd_vs_gd_PC1_p": float(p_pcs), "rd_vs_gd_PC2_p": float(p_pcs2),
            "kmeans4_sizes": {int(k): int(v) for k, v in Counter(km4).items()},
            "rd_best_cluster": int(rd_best) if rd_best is not None else None,
            "rd_in_best_cluster": int(ct.loc["regulation-driven", rd_best]) if rd_best is not None else None},
    "bootstrap_posteriors": {"median": float(np.median(stability)), "p_ge_0.8": float(np.mean(stability >= 0.8))},
    "evidence_coverage": cov,
    "rd_yield_with_campra": {"rd": rd_with_campra, "n": nonsig_with_campra,
                             "pct": 100 * rd_with_campra / max(nonsig_with_campra, 1)},
    "counterfactual_campra_universal": {"rd_baseline": int((baseline_cls == "regulation-driven").sum()),
                                        "rd_if_all_campra1": rd_hypo},
    "relax_coverage_by_class": {c: {"active": int(relax_active[baseline_cls == c].sum()),
                                    "total": int((baseline_cls == c).sum())}
                                for c in set(baseline_cls)},
}
with open(f"{OUT}/revision_v7_summary.json", "w") as f:
    json.dump(summary, f, indent=2, default=float)
print(f"\nSaved: {OUT}/revision_v7_summary.json")
print("=" * 70)
print("Revision v7 analyses complete")
print("=" * 70)
