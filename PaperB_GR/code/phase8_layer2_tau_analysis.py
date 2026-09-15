#!/usr/bin/env python3
"""
Phase 8 Component 3: Expression Specificity (tau) Analysis

Multi-source tau comparison: GTEx V11 + HPA (3 sources)
Non-circular: uses classification_loo_tau
Reports: Kruskal-Wallis, Mann-Whitney U, Cliff's delta, regression

Run in WSL: python3 phase8_layer2_tau_analysis.py
"""
import csv
import json
import numpy as np
from scipy import stats
from pathlib import Path

import os
BASE = Path(os.environ.get("HSD_BASE") or (
    "/mnt/d/hs_gene_project" if os.path.exists("/mnt/d/hs_gene_project") else "D:/人类正选择基因项目"))  # HSD_BASE-aware
MASTER_CSV = BASE / "results/phase8_tissue_analysis/data/phase8_master_table.csv"
OUTPUT_DIR = BASE / "results/phase8_tissue_analysis"

print("=" * 70)
print("Phase 8 Component 3: Expression Specificity (tau) Analysis")
print("=" * 70)

# ============================================
# Load master table
# ============================================
print("\n--- Loading master table ---")

genes = []
with open(MASTER_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        genes.append(row)

print(f"  Loaded {len(genes)} genes")

# ============================================
# Define tau sources
# ============================================
tau_sources = {
    "GTEx V11 (30 tissues)": "gtex_v11_tau",
    "HPA GTEx (36 tissues)": "hpa_tau_gtex",
    "HPA HPA (40 tissues)": "hpa_tau_hpa",
    "HPA Consensus (51 tissues)": "hpa_tau_consensus",
}

# ============================================
# Define gene groups
# ============================================
print("\n--- Defining gene groups ---")

# LOO-tau classification (non-circular for tau analysis)
gd_loo_tau = [g for g in genes if g.get("classification_loo_tau") == "gene-driven"]
rd_loo_tau = [g for g in genes if g.get("classification_loo_tau") == "regulation-driven"]
dual_loo_tau = [g for g in genes if g.get("classification_loo_tau") == "dual-driven"]
gd_relaxed_loo_tau = [g for g in genes if g.get("classification_loo_tau") == "gene-driven (relaxed)"]
neutral_loo_tau = [g for g in genes if g.get("classification_loo_tau") == "neutral"]

# Full v7 classification (for circular reference comparison)
gd_full = [g for g in genes if g["classification_v7"] == "gene-driven"]
rd_full = [g for g in genes if g["classification_v7"] == "regulation-driven"]
neutral_full = [g for g in genes if g["classification_v7"] == "neutral"]

print(f"  LOO-tau: GD={len(gd_loo_tau)}, RD={len(rd_loo_tau)}, dual={len(dual_loo_tau)}, "
      f"GD(relaxed)={len(gd_relaxed_loo_tau)}, neutral={len(neutral_loo_tau)}")
print(f"  Full v7: GD={len(gd_full)}, RD={len(rd_full)}, neutral={len(neutral_full)}")

# ============================================
# Helper functions
# ============================================

def cliffs_delta(x, y):
    """Compute Cliff's delta effect size."""
    x, y = np.array(x), np.array(y)
    n = len(x) * len(y)
    dominance = 0
    for xi in x:
        dominance += np.sum(xi > y) - np.sum(xi < y)
    return dominance / n

def bh_fdr(pvals):
    pvals = np.array(pvals)
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    fdr = np.zeros(n)
    fdr[sorted_idx[-1]] = pvals[sorted_idx[-1]]
    for i in range(n - 2, -1, -1):
        rank = i + 1
        fdr[sorted_idx[i]] = min(pvals[sorted_idx[i]] * n / rank, fdr[sorted_idx[i + 1]])
    return fdr

# ============================================
# Layer 2: Tau analysis
# ============================================
print("\n--- Layer 2: Tau distribution analysis ---")

results_tau = []
results_tests = []

for source_name, col in tau_sources.items():
    print(f"\n  === {source_name} ({col}) ===")
    
    # Extract tau values for each group
    def get_tau(group):
        vals = []
        for g in group:
            try:
                v = float(g.get(col) or 0)
                if np.isfinite(v):
                    vals.append(v)
            except (ValueError, TypeError):
                pass
        return vals
    
    # LOO-tau groups (non-circular)
    gd_tau_loo = get_tau(gd_loo_tau)
    rd_tau_loo = get_tau(rd_loo_tau)
    dual_tau_loo = get_tau(dual_loo_tau)
    gd_relaxed_tau_loo = get_tau(gd_relaxed_loo_tau)
    neutral_tau_loo = get_tau(neutral_loo_tau)
    
    # Full v7 groups (circular reference)
    gd_tau_circ = get_tau(gd_full)
    rd_tau_circ = get_tau(rd_full)
    neutral_tau_circ = get_tau(neutral_full)
    
    # Kruskal-Wallis across 5 groups (LOO-tau)
    all_groups_loo = [gd_tau_loo, rd_tau_loo, dual_tau_loo, gd_relaxed_tau_loo, neutral_tau_loo]
    all_groups_loo_filtered = [g for g in all_groups_loo if len(g) > 0]
    if len(all_groups_loo_filtered) >= 2:
        kw_stat, kw_p = stats.kruskal(*all_groups_loo_filtered)
    else:
        kw_stat, kw_p = None, None
    
    print(f"    Kruskal-Wallis (5 groups, LOO): H={kw_stat:.2f}, p={kw_p:.2e}")
    
    # Mann-Whitney U: RD vs GD (LOO, non-circular)
    # 2026-09-14: switched one-sided (greater) -> two-sided to match manuscript
    # Methods claim (manuscriptB_v2 reports two-sided values)
    if len(rd_tau_loo) > 0 and len(gd_tau_loo) > 0:
        u_loo, p_loo = stats.mannwhitneyu(rd_tau_loo, gd_tau_loo, alternative="two-sided")
        delta_loo = cliffs_delta(rd_tau_loo, gd_tau_loo)
        print(f"    Mann-Whitney RD(LOO) > GD(LOO): U={u_loo:.0f}, p={p_loo:.2e}")
        print(f"    Cliff's delta: {delta_loo:.4f}")
        print(f"    RD(LOO) tau: median={np.median(rd_tau_loo):.4f}, mean={np.mean(rd_tau_loo):.4f}")
        print(f"    GD(LOO) tau: median={np.median(gd_tau_loo):.4f}, mean={np.mean(gd_tau_loo):.4f}")
    else:
        u_loo, p_loo, delta_loo = None, None, None
    
    # Circular reference (full v7) — two-sided since 2026-09-14 (see above)
    if len(rd_tau_circ) > 0 and len(gd_tau_circ) > 0:
        u_circ, p_circ = stats.mannwhitneyu(rd_tau_circ, gd_tau_circ, alternative="two-sided")
        delta_circ = cliffs_delta(rd_tau_circ, gd_tau_circ)
        print(f"    [CIRCULAR REF] RD(full) > GD(full): U={u_circ:.0f}, p={p_circ:.2e}")
        print(f"    [CIRCULAR REF] Cliff's delta: {delta_circ:.4f}")
    else:
        u_circ, p_circ, delta_circ = None, None, None
    
    # Store distribution stats
    for group_name, group_data, tau_vals in [
        ("gene-driven", gd_loo_tau, gd_tau_loo),
        ("regulation-driven", rd_loo_tau, rd_tau_loo),
        ("dual-driven", dual_loo_tau, dual_tau_loo),
        ("gene-driven(relaxed)", gd_relaxed_loo_tau, gd_relaxed_tau_loo),
        ("neutral", neutral_loo_tau, neutral_tau_loo),
    ]:
        if tau_vals:
            results_tau.append({
                "tau_source": source_name,
                "classification": group_name,
                "n": len(tau_vals),
                "mean": np.mean(tau_vals),
                "median": np.median(tau_vals),
                "std": np.std(tau_vals),
                "q25": np.percentile(tau_vals, 25),
                "q75": np.percentile(tau_vals, 75),
            })
    
    # Store test results
    results_tests.append({
        "tau_source": source_name,
        "kw_stat": kw_stat,
        "kw_p": kw_p,
        "mw_u_loo": u_loo,
        "mw_p_loo": p_loo,
        "cliffs_delta_loo": delta_loo,
        "rd_median_loo": np.median(rd_tau_loo) if rd_tau_loo else None,
        "gd_median_loo": np.median(gd_tau_loo) if gd_tau_loo else None,
        "mw_u_circular": u_circ,
        "mw_p_circular": p_circ,
        "cliffs_delta_circular": delta_circ,
        "rd_median_circular": np.median(rd_tau_circ) if rd_tau_circ else None,
        "gd_median_circular": np.median(gd_tau_circ) if gd_tau_circ else None,
        "non_circular": True,
    })

# ============================================
# Regression: tau ~ classification + log10(CDS length)
# ============================================
print(f"\n--- Regression: tau ~ classification + log10(CDS) ---")

results_regression = []

for source_name, col in tau_sources.items():
    # Use LOO-tau classification
    # Encode: GD=1, RD=2, neutral=0
    x_class = []
    y_tau = []
    x_cds = []
    
    for g in genes:
        cls = g.get("classification_loo_tau", "neutral")
        try:
            tau_val = float(g.get(col) or 0)
            cds_val = float(g.get("cds_length") or 0)
            if np.isfinite(tau_val) and cds_val > 0:
                y_tau.append(tau_val)
                x_cds.append(np.log10(cds_val))
                if cls == "gene-driven":
                    x_class.append(1)
                elif cls == "regulation-driven":
                    x_class.append(2)
                else:
                    x_class.append(0)
        except (ValueError, TypeError):
            pass
    
    y_tau = np.array(y_tau)
    x_class = np.array(x_class)
    x_cds = np.array(x_cds)
    
    if len(y_tau) > 100:
        # Multiple regression: tau ~ class + log10(CDS)
        X = np.column_stack([x_class, x_cds, np.ones(len(y_tau))])
        try:
            beta, residuals, rank, sv = np.linalg.lstsq(X, y_tau, rcond=None)
            y_pred = X @ beta
            ss_res = np.sum((y_tau - y_pred) ** 2)
            ss_tot = np.sum((y_tau - np.mean(y_tau)) ** 2)
            r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            
            # t-test for class coefficient
            n = len(y_tau)
            k = 3
            mse = ss_res / (n - k) if n > k else 0
            se = np.sqrt(np.diag(np.linalg.inv(X.T @ X) * mse)) if mse > 0 else [0, 0, 0]
            t_stat = beta[0] / se[0] if se[0] > 0 else 0
            p_val = 2 * stats.t.sf(np.abs(t_stat), df=n-k)
            
            print(f"  {source_name}: R²={r_squared:.4f}, beta_class={beta[0]:.4f}, "
                  f"t={t_stat:.2f}, p={p_val:.2e}")
            
            results_regression.append({
                "tau_source": source_name,
                "n": len(y_tau),
                "r_squared": r_squared,
                "beta_class": beta[0],
                "beta_cds": beta[1],
                "beta_intercept": beta[2],
                "t_stat_class": t_stat,
                "p_val_class": p_val,
            })
        except Exception as e:
            print(f"  {source_name}: regression error: {e}")

# ============================================
# Per-tissue expression specificity
# ============================================
print(f"\n--- Per-tissue expression specificity ---")

# For each GTEx tissue, compute specificity = TPM(tissue) / max(TPM(all tissues))
gtex_tissue_cols = [h for h in genes[0].keys() if h.startswith("gtex_") and h != "gtex_v11_tau"]

per_tissue_results = []

for tissue_col in gtex_tissue_cols:
    tissue_name = tissue_col.replace("gtex_", "")
    
    # For each gene, compute tissue specificity
    tissue_specs_gd = []
    tissue_specs_rd = []
    tissue_specs_neutral = []
    
    # Use LOO-tau classification
    for g in genes:
        cls = g.get("classification_loo_tau", "neutral")
        
        # Get all tissue TPMs for this gene
        all_tpms = []
        for tc in gtex_tissue_cols:
            try:
                v = float(g.get(tc) or 0)
                all_tpms.append(v)
            except:
                all_tpms.append(0)
        
        max_tpm = max(all_tpms) if all_tpms else 0
        tissue_tpm = float(g.get(tissue_col) or 0)
        spec = tissue_tpm / max_tpm if max_tpm > 0 else 0
        
        if cls == "gene-driven":
            tissue_specs_gd.append(spec)
        elif cls == "regulation-driven":
            tissue_specs_rd.append(spec)
        elif cls == "neutral":
            tissue_specs_neutral.append(spec)
    
    # Mann-Whitney: RD vs GD for tissue specificity — two-sided since 2026-09-14
    # (was one-sided "greater"; manuscript Methods declares two-sided)
    # Metric definition (for Methods): spec_g(t) = TPM_g(t) / max_t' TPM_g(t'),
    # i.e. each tissue's share of the gene's maximum TPM across the 30 GTEx
    # tissues (0 when the gene is unexpressed everywhere); LOO-tau classes.
    if len(tissue_specs_rd) > 5 and len(tissue_specs_gd) > 5:
        u, p = stats.mannwhitneyu(tissue_specs_rd, tissue_specs_gd, alternative="two-sided")
        delta = cliffs_delta(tissue_specs_rd, tissue_specs_gd)
        
        per_tissue_results.append({
            "tissue": tissue_name,
            "gd_mean_spec": np.mean(tissue_specs_gd),
            "rd_mean_spec": np.mean(tissue_specs_rd),
            "neutral_mean_spec": np.mean(tissue_specs_neutral),
            "mw_u": u,
            "mw_p": p,
            "cliffs_delta": delta,
            "rd_gt_gd": p < 0.05,
        })

# BH FDR for per-tissue tests
if per_tissue_results:
    pvals = [r["mw_p"] for r in per_tissue_results]
    fdrs = bh_fdr(pvals)
    for i, r in enumerate(per_tissue_results):
        r["fdr"] = fdrs[i]
        r["fdr_sig"] = fdrs[i] < 0.05

# Print significant per-tissue results
print(f"\n  Tissues where RD > GD in tissue specificity (FDR < 0.05):")
sig_tissues = [r for r in per_tissue_results if r.get("fdr_sig", False)]
for r in sorted(sig_tissues, key=lambda x: x["fdr"]):
    print(f"    {r['tissue']:25s} RD_spec={r['rd_mean_spec']:.4f} "
          f"GD_spec={r['gd_mean_spec']:.4f} "
          f"delta={r['cliffs_delta']:.4f} "
          f"FDR={r['fdr']:.2e}")

if not sig_tissues:
    print(f"    (none FDR-significant)")
    # Show top 5 nominal
    for r in sorted(per_tissue_results, key=lambda x: x["mw_p"])[:5]:
        print(f"    [nominal] {r['tissue']:25s} RD={r['rd_mean_spec']:.4f} "
              f"GD={r['gd_mean_spec']:.4f} p={r['mw_p']:.2e}")

# ============================================
# Write outputs
# ============================================
print(f"\n--- Writing outputs ---")

# Tau distributions
tau_dist_csv = OUTPUT_DIR / "layer2_tau_distributions.csv"
with open(tau_dist_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["tau_source", "classification", "n", "mean", "median", "std", "q25", "q75"])
    writer.writeheader()
    for r in results_tau:
        writer.writerow(r)
print(f"  {tau_dist_csv}")

# Statistical tests
tests_csv = OUTPUT_DIR / "layer2_statistical_tests.csv"
with open(tests_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(results_tests[0].keys()))
    writer.writeheader()
    for r in results_tests:
        writer.writerow(r)
print(f"  {tests_csv}")

# Regression results
reg_csv = OUTPUT_DIR / "layer2_regression_results.csv"
if results_regression:
    with open(reg_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results_regression[0].keys()))
        writer.writeheader()
        for r in results_regression:
            writer.writerow(r)
    print(f"  {reg_csv}")

# Per-tissue specificity
per_tissue_csv = OUTPUT_DIR / "layer2_per_tissue_specificity.csv"
if per_tissue_results:
    with open(per_tissue_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(per_tissue_results[0].keys()))
        writer.writeheader()
        for r in per_tissue_results:
            writer.writerow(r)
    print(f"  {per_tissue_csv}")

print("\n" + "=" * 70)
print("Layer 2 Tau Analysis Complete!")
print("=" * 70)
