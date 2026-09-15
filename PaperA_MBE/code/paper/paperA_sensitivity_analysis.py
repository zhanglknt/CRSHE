#!/usr/bin/env python3
"""
Paper A (MBE) — P0: Weight & Threshold Sensitivity Analysis
============================================================
Systematic robustness analysis of the GD/RD classification (Phase 7G v7):

1. Threshold grid: 0.05–0.30, classification counts + agreement with baseline
2. Single-weight ±50% perturbation (8 weights individually)
3. Random joint perturbation: 1000 runs, GDS+RDS weights ~ U(0.5, 1.5) multiplier, renormalized
4. Random joint perturbation + threshold ~ U(0.05, 0.30)
5. Per-gene stability score + class-level confusion matrix + figures

Reads component scores from gene_classification_v7.csv (no re-derivation).
Output: results/paper/sensitivity/  (CSV + JSON + 2 figures)

Run in WSL: cd /mnt/d/人类正选择基因项目 && python3 scripts/paperA_sensitivity_analysis.py
"""
import json
import os
import sys
import numpy as np
import pandas as pd
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Path auto-detection: WSL (/mnt/d/...) or Windows (D:/...)
if os.path.exists("/mnt/d"):
    BASE = "/mnt/d/人类正选择基因项目"
else:
    BASE = "D:/人类正选择基因项目"
INPUT_CSV = f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv"
OUT_DIR = f"{BASE}/results/paper/sensitivity"
os.makedirs(OUT_DIR, exist_ok=True)

rng = np.random.default_rng(42)

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
print("=" * 70)
print("Paper A — Weight & Threshold Sensitivity Analysis")
print("=" * 70)

df = pd.read_csv(INPUT_CSV, low_memory=False)
print(f"Loaded: {INPUT_CSV}  ({len(df)} genes)")

GDS_COLS = ["gds_p_pct_resid", "gds_lrt_pct_resid", "gds_relax_pct", "gds_selectome"]
RDS_COLS = ["rds_doan", "rds_tau", "rds_nc", "rds_brain"]
GDS_W = np.array([0.35, 0.30, 0.20, 0.15])
RDS_W = np.array([0.30, 0.25, 0.25, 0.20])
THRESHOLD = 0.15
GDS_LABELS = ["P-value (resid.)", "LRT (resid.)", "RELAX K", "Selectome"]
RDS_LABELS = ["caMPRA", "tau", "nc-conservation", "brain tau"]

X_gds = df[GDS_COLS].apply(pd.to_numeric, errors="coerce").fillna(0).values
X_rds = df[RDS_COLS].apply(pd.to_numeric, errors="coerce").fillna(0).values
sig = df["bh_fdr_sig"].astype(str).str.strip().str.lower().eq("true").values
relaxed = df["relax_relaxed"].astype(str).str.strip().str.lower().eq("true").values
baseline_cls = df["classification_v7"].astype(str).values

# ------------------------------------------------------------------
# Classification function (exact replica of phase7g_v7_classification.py)
# ------------------------------------------------------------------
def classify(gds, rds, sig_arr, relaxed_arr, threshold=0.15):
    out = np.empty(len(gds), dtype=object)
    for i in range(len(gds)):
        g, r = gds[i], rds[i]
        s = sig_arr[i]
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
        if relaxed_arr[i] and c == "gene-driven":
            c = "gene-driven (relaxed)"
        out[i] = c
    return out

def counts(cls_arr):
    c = Counter(cls_arr)
    return {k: int(c.get(k, 0)) for k in
            ["gene-driven", "gene-driven (relaxed)", "regulation-driven", "dual-driven", "neutral"]}

# ------------------------------------------------------------------
# Baseline reproduction check
# ------------------------------------------------------------------
gds_base = X_gds @ GDS_W
rds_base = X_rds @ RDS_W
cls_base = classify(gds_base, rds_base, sig, relaxed, THRESHOLD)
repro_agree = np.mean(cls_base == baseline_cls)
print(f"\nBaseline reproduction agreement with stored classification_v7: {100*repro_agree:.2f}%")
base_counts = counts(cls_base)
print(f"Baseline counts: {base_counts}")
assert repro_agree == 1.0 or repro_agree > 0.999, "Baseline reproduction failed — check component columns!"

# ------------------------------------------------------------------
# 1) Threshold grid
# ------------------------------------------------------------------
print("\n--- 1) Threshold grid ---")
thresholds = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
thr_rows = []
for t in thresholds:
    cls_t = classify(gds_base, rds_base, sig, relaxed, t)
    agree = np.mean(cls_t == cls_base)
    row = {"threshold": t, "agreement_with_baseline": agree, **counts(cls_t)}
    thr_rows.append(row)
    print(f"  t={t:.2f}: GD={row['gene-driven']}, RD={row['regulation-driven']}, "
          f"dual={row['dual-driven']}, neutral={row['neutral']}, agree={100*agree:.1f}%")
thr_df = pd.DataFrame(thr_rows)
thr_df.to_csv(f"{OUT_DIR}/threshold_sensitivity.csv", index=False)

# ------------------------------------------------------------------
# 2) Single-weight ±50% (individual)
# ------------------------------------------------------------------
print("\n--- 2) Single-weight ±50% perturbation ---")
sw_rows = []
for side, cols, w0, labels, X in [("GDS", GDS_COLS, GDS_W, GDS_LABELS, X_gds),
                                   ("RDS", RDS_COLS, RDS_W, RDS_LABELS, X_rds)]:
    for j, lab in enumerate(labels):
        for mult, tag in [(1.5, "+50%"), (0.5, "-50%")]:
            w = w0.copy()
            w[j] = w0[j] * mult
            w = w / w.sum()
            if side == "GDS":
                gds_s, rds_s = X_gds @ w, rds_base
            else:
                gds_s, rds_s = gds_base, X_rds @ w
            cls_s = classify(gds_s, rds_s, sig, relaxed, THRESHOLD)
            agree = np.mean(cls_s == cls_base)
            row = {"score": side, "component": lab, "perturbation": tag,
                   "agreement_with_baseline": agree, **counts(cls_s)}
            sw_rows.append(row)
            print(f"  {side}.{lab} {tag}: GD={row['gene-driven']}, RD={row['regulation-driven']}, "
                  f"agree={100*agree:.1f}%")
sw_df = pd.DataFrame(sw_rows)
sw_df.to_csv(f"{OUT_DIR}/single_weight_sensitivity.csv", index=False)

# ------------------------------------------------------------------
# 3) Joint random perturbation (1000 runs)
# ------------------------------------------------------------------
print("\n--- 3) Joint random weight perturbation (1000 runs) ---")
N_RUNS = 1000
n = len(df)
pert_classes = np.empty((N_RUNS, n), dtype=object)
for run in range(N_RUNS):
    wg = GDS_W * rng.uniform(0.5, 1.5, size=4)
    wg = wg / wg.sum()
    wr = RDS_W * rng.uniform(0.5, 1.5, size=4)
    wr = wr / wr.sum()
    pert_classes[run] = classify(X_gds @ wg, X_rds @ wr, sig, relaxed, THRESHOLD)

agree_runs = np.array([np.mean(pert_classes[k] == cls_base) for k in range(N_RUNS)])
print(f"  Agreement with baseline: mean={100*agree_runs.mean():.1f}%, "
      f"sd={100*agree_runs.std():.1f}%, range=[{100*agree_runs.min():.1f}%, {100*agree_runs.max():.1f}%]")

pert_counts = pd.DataFrame([counts(pert_classes[k]) for k in range(N_RUNS)])
pc_summary = {c: {"mean": float(pert_counts[c].mean()), "sd": float(pert_counts[c].std()),
                  "p5": float(pert_counts[c].quantile(0.05)), "p95": float(pert_counts[c].quantile(0.95))}
              for c in pert_counts.columns}
base_c = counts(cls_base)
print("  Class counts across perturbations (mean±SD [5-95%]) vs baseline:")
for c in ["gene-driven", "regulation-driven", "dual-driven", "neutral"]:
    s = pc_summary[c]
    print(f"    {c:22s}: {s['mean']:7.1f} ± {s['sd']:5.1f} [{s['p5']:.0f}-{s['p95']:.0f}]  (baseline {base_c.get(c, 0)})")

# Per-gene stability
stability = np.mean(pert_classes == cls_base[None, :], axis=0)
stab_df = pd.DataFrame({"gene_id": df["gene_id"].values, "gene_symbol": df.get("gene_symbol", pd.Series([""]*n)).values,
                        "classification_v7": cls_base, "stability": stability})
stab_df.to_csv(f"{OUT_DIR}/per_gene_stability.csv", index=False)
print(f"\n  Per-gene stability: median={np.median(stability):.3f}")
for c in ["gene-driven", "regulation-driven", "dual-driven", "neutral"]:
    m = cls_base == c
    if m.sum() > 0:
        print(f"    {c:22s}: median stability={np.median(stability[m]):.3f} (n={m.sum()})")
stable_core = np.mean(stability >= 0.8)
print(f"  Genes with stability >= 0.80 ('stable core'): {100*stable_core:.1f}%")

# ------------------------------------------------------------------
# 4) Joint perturbation + threshold jitter
# ------------------------------------------------------------------
print("\n--- 4) Joint perturbation + threshold jitter (1000 runs) ---")
thr_jit = rng.uniform(0.05, 0.30, size=N_RUNS)
agree_jit = np.zeros(N_RUNS)
for run in range(N_RUNS):
    wg = GDS_W * rng.uniform(0.5, 1.5, size=4); wg = wg / wg.sum()
    wr = RDS_W * rng.uniform(0.5, 1.5, size=4); wr = wr / wr.sum()
    cj = classify(X_gds @ wg, X_rds @ wr, sig, relaxed, thr_jit[run])
    agree_jit[run] = np.mean(cj == cls_base)
print(f"  Agreement: mean={100*agree_jit.mean():.1f}%, sd={100*agree_jit.std():.1f}%, "
      f"range=[{100*agree_jit.min():.1f}%, {100*agree_jit.max():.1f}%]")

# ------------------------------------------------------------------
# 5) Figures
# ------------------------------------------------------------------
print("\n--- 5) Figures ---")
CLASSES = ["gene-driven", "gene-driven (relaxed)", "regulation-driven", "dual-driven", "neutral"]
COLORS = {"gene-driven": "#2166ac", "gene-driven (relaxed)": "#92c5de",
          "regulation-driven": "#d6604d", "dual-driven": "#f4a582", "neutral": "#bdbdbd"}

# Fig S1: threshold sensitivity
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
bottom = np.zeros(len(thr_df))
for c in CLASSES:
    vals = thr_df[c].values.astype(float)
    ax.bar(thr_df["threshold"], vals, bottom=bottom, color=COLORS[c], label=c, width=0.035)
    bottom += vals
ax.axvline(THRESHOLD, color="black", ls="--", lw=1, label=f"baseline ({THRESHOLD})")
ax.set_xlabel("Classification threshold")
ax.set_ylabel("Number of genes")
ax.set_title("(a) Class composition vs threshold")
ax.legend(fontsize=7, loc="upper right")

ax = axes[1]
ax.plot(thr_df["threshold"], 100*thr_df["agreement_with_baseline"], "o-", color="#1a9850")
ax.axvline(THRESHOLD, color="black", ls="--", lw=1)
ax.set_xlabel("Classification threshold")
ax.set_ylabel("Agreement with baseline (%)")
ax.set_title("(b) Agreement vs threshold")
ax.set_ylim(0, 105)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/figS_threshold_sensitivity.png", dpi=300)
fig.savefig(f"{OUT_DIR}/figS_threshold_sensitivity.pdf", bbox_inches="tight")
plt.close(fig)

# Fig S2: perturbation stability
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
# (a) confusion matrix baseline vs perturbed (mode per gene)
mode_cls = np.array([Counter(pert_classes[:, i]).most_common(1)[0][0] for i in range(n)])
labels = CLASSES
M = np.zeros((len(labels), len(labels)), dtype=float)
for i, ci in enumerate(labels):
    for j, cj in enumerate(labels):
        M[i, j] = np.sum((cls_base == ci) & (mode_cls == cj))
Mnorm = M / M.sum(axis=1, keepdims=True) * 100
ax = axes[0]
im = ax.imshow(Mnorm, cmap="Blues", vmin=0, vmax=100)
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
for i in range(len(labels)):
    for j in range(len(labels)):
        ax.text(j, i, f"{Mnorm[i, j]:.0f}", ha="center", va="center", fontsize=7,
                color="white" if Mnorm[i, j] > 50 else "black")
ax.set_xlabel("Most frequent perturbed class")
ax.set_ylabel("Baseline class (v7)")
ax.set_title("(a) Baseline vs perturbed class (%)")
fig.colorbar(im, ax=ax, shrink=0.7)

# (b) per-gene stability histogram by class
ax = axes[1]
for c in CLASSES:
    m = cls_base == c
    if m.sum() > 0:
        ax.hist(stability[m], bins=20, alpha=0.55, label=f"{c} (n={m.sum()})", color=COLORS[c])
ax.set_xlabel("Per-gene stability (fraction of 1000 runs keeping baseline class)")
ax.set_ylabel("Genes")
ax.set_title("(b) Per-gene classification stability")
ax.legend(fontsize=7)

# (c) agreement distribution across runs
ax = axes[2]
ax.hist(100*agree_runs, bins=30, alpha=0.7, color="#2166ac", label="weights only")
ax.hist(100*agree_jit, bins=30, alpha=0.7, color="#d6604d", label="weights + threshold")
ax.set_xlabel("Agreement with baseline classification (%)")
ax.set_ylabel("Runs")
ax.set_title("(c) Run-level agreement")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/figS_weight_perturbation.png", dpi=300)
fig.savefig(f"{OUT_DIR}/figS_weight_perturbation.pdf", bbox_inches="tight")
plt.close(fig)
print("  Saved: figS_threshold_sensitivity.png, figS_weight_perturbation.png")

# ------------------------------------------------------------------
# 6) Summary JSON
# ------------------------------------------------------------------
summary = {
    "baseline": {
        "gds_weights": dict(zip(GDS_LABELS, GDS_W.tolist())),
        "rds_weights": dict(zip(RDS_LABELS, RDS_W.tolist())),
        "threshold": THRESHOLD,
        "counts": base_counts,
        "reproduction_agreement": float(repro_agree),
    },
    "threshold_grid": thr_rows,
    "single_weight": sw_rows,
    "joint_perturbation_1000": {
        "agreement_mean": float(agree_runs.mean()), "agreement_sd": float(agree_runs.std()),
        "agreement_min": float(agree_runs.min()), "agreement_max": float(agree_runs.max()),
        "class_counts": pc_summary,
    },
    "joint_perturbation_with_threshold_jitter": {
        "agreement_mean": float(agree_jit.mean()), "agreement_sd": float(agree_jit.std()),
        "threshold_range": [0.05, 0.30],
    },
    "per_gene_stability": {
        "median": float(np.median(stability)),
        "stable_core_ge_0.8_fraction": float(stable_core),
        "by_class": {c: float(np.median(stability[cls_base == c])) for c in CLASSES if (cls_base == c).sum() > 0},
    },
}
with open(f"{OUT_DIR}/sensitivity_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"\nSaved: {OUT_DIR}/sensitivity_summary.json")
print("=" * 70)
print("Sensitivity analysis complete")
print("=" * 70)
