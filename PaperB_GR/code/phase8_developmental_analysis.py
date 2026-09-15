#!/usr/bin/env python3
"""
Phase 8 Component 6: Full Developmental Trajectory Analysis

BrainSpan 31 developmental stages — full trajectory, not just 4 periods.
- Per-stage GD/RD enrichment
- Developmental tau across stages
- Prenatal vs postnatal comparison
- Region × stage matrix

Uses classification_loo_brain for non-circular brain analysis.

Run in WSL: python3 phase8_developmental_analysis.py
"""
import csv
import json
import os
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path

BASE = Path(os.environ.get("HSD_BASE") or (
    "/mnt/d/hs_gene_project" if os.path.exists("/mnt/d/hs_gene_project") else "D:/人类正选择基因项目"))
MASTER_CSV = BASE / "results/phase8_tissue_analysis/data/phase8_master_table.csv"
BRAINSPAN_EXPR = BASE / "data/brainspan/expression_matrix.csv"
BRAINSPAN_COLS = BASE / "data/brainspan/columns_metadata.csv"
BRAINSPAN_ROWS = BASE / "data/brainspan/rows_metadata.csv"
OUTPUT_DIR = BASE / "results/phase8_tissue_analysis"

print("=" * 70)
print("Phase 8 Component 6: Full Developmental Trajectory Analysis")
print("=" * 70)

# ============================================
# Helper functions
# ============================================
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

def fisher_test(group_set, bg_set, feature_set):
    a = len(group_set & feature_set)
    b = len(group_set - feature_set)
    c = len(feature_set - group_set)
    d = len(bg_set - group_set - feature_set)
    if a + b == 0 or c + d == 0:
        return 1.0, 1.0, 0, 0
    table = [[a, b], [c, d]]
    odds, p = stats.fisher_exact(table, alternative="greater")
    return odds, p, a, len(group_set)

def compute_tau(expr_values):
    """Compute tissue specificity index tau."""
    vals = np.array(expr_values, dtype=float)
    if np.all(vals == 0) or np.max(vals) == 0:
        return 0.0
    if len(vals) < 2:
        return 0.0
    max_val = np.max(vals)
    if max_val == 0:
        return 0.0
    normalized = vals / max_val
    return 1.0 - np.mean(normalized)

def parse_age_to_weeks(age_str):
    """Parse BrainSpan age string to gestational weeks for sorting."""
    if not age_str:
        return 9999
    age_str = str(age_str).strip('"').strip()
    if "pcw" in age_str:
        try:
            return float(age_str.replace("pcw", "").strip())
        except:
            return 9999
    elif "mos" in age_str or "months" in age_str:
        try:
            months = float(age_str.replace("mos", "").replace("months", "").strip())
            return 40 + months * 4.345  # Convert to gestational weeks
        except:
            return 9999
    elif "yrs" in age_str or "years" in age_str:
        try:
            years = float(age_str.replace("yrs", "").replace("years", "").strip())
            return 40 + years * 52  # Convert to gestational weeks
        except:
            return 9999
    return 9999

def get_dev_period(age_str):
    """Classify age into developmental period."""
    weeks = parse_age_to_weeks(age_str)
    if weeks <= 17:
        return "early_prenatal"
    elif weeks <= 38:
        return "late_prenatal"
    elif weeks <= 40 + 5 * 52:  # ~5 years
        return "early_postnatal"
    else:
        return "late_postnatal"

# ============================================
# Load master table
# ============================================
print("\n--- Loading master table ---")
genes = {}
with open(MASTER_CSV) as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        genes[row["gene_id"]] = row

# LOO-brain classification (non-circular)
gd_genes = set(g for g, v in genes.items() if v.get("classification_loo_brain") == "gene-driven")
rd_genes = set(g for g, v in genes.items() if v.get("classification_loo_brain") == "regulation-driven")
neutral_genes = set(g for g, v in genes.items() if v.get("classification_loo_brain") == "neutral")
all_genes = set(genes.keys())

print(f"  GD (LOO-brain): {len(gd_genes)}")
print(f"  RD (LOO-brain): {len(rd_genes)}")
print(f"  Neutral (LOO-brain): {len(neutral_genes)}")

# ============================================
# Load BrainSpan metadata
# ============================================
print("\n--- Loading BrainSpan metadata ---")
columns = []
with open(BRAINSPAN_COLS) as f:
    reader = csv.DictReader(f)
    for row in reader:
        columns.append(row)
print(f"  {len(columns)} samples")

# Parse sample info: age, region, dev period
sample_info = []
for col in columns:
    age = col.get("age", "").strip('"')
    region = col.get("structure_acronym", col.get("structure_name", "")).strip('"')
    period = get_dev_period(age)
    weeks = parse_age_to_weeks(age)
    sample_info.append({"age": age, "region": region, "period": period, "weeks": weeks})

# Get unique stages sorted by developmental order
unique_ages = sorted(set(s["age"] for s in sample_info if s["age"]),
                     key=parse_age_to_weeks)
print(f"  Unique developmental stages: {len(unique_ages)}")
for a in unique_ages:
    n = sum(1 for s in sample_info if s["age"] == a)
    print(f"    {a}: {n} samples, ~{parse_age_to_weeks(a):.0f} weeks")

# Unique brain regions
unique_regions = sorted(set(s["region"] for s in sample_info if s["region"]))
print(f"  Brain regions: {len(unique_regions)}")

# ============================================
# Load gene mapping
# ============================================
print("\n--- Loading BrainSpan gene mapping ---")
row_map = {}
with open(BRAINSPAN_ROWS) as f:
    reader = csv.DictReader(f)
    for i, row in enumerate(reader):
        eid = row.get("ensembl_gene_id", row.get("gene_id", ""))
        if eid:
            eid = eid.split(".")[0]
            row_map[i] = eid

row_map_filtered = {k: v for k, v in row_map.items() if v in genes}
print(f"  {len(row_map)} genes mapped, {len(row_map_filtered)} overlap with our set")

# ============================================
# Load expression matrix (only for our genes)
# ============================================
print("\n--- Loading BrainSpan expression matrix (183MB) ---")
print("  Reading only rows matching our 4,974 gene set...")

gene_expr = {}  # ensembl_id -> [rpkm values across all samples]
with open(BRAINSPAN_EXPR) as f:
    reader = csv.reader(f)
    header = next(reader)
    n_samples = len(header) - 1
    for i, row in enumerate(reader):
        if i in row_map_filtered:
            eid = row_map_filtered[i]
            try:
                expr_vals = [float(v) if v else 0.0 for v in row[1:]]
                gene_expr[eid] = expr_vals
            except ValueError:
                pass

print(f"  Expression data loaded for {len(gene_expr)} genes, {n_samples} samples")

# ============================================
# Analysis 1: Per-Stage Enrichment (31 stages)
# ============================================
print(f"\n{'='*50}")
print("Analysis 1: Per-Stage Enrichment (31 developmental stages)")
print(f"{'='*50}")

stage_results = []

for age in unique_ages:
    # Get sample indices for this stage
    stage_indices = [i for i, s in enumerate(sample_info) if s["age"] == age]
    
    if len(stage_indices) < 3:
        stage_results.append({
            "stage": age, "n_samples": len(stage_indices),
            "gd_OR": 1.0, "gd_p": 1.0, "gd_fdr": 1.0,
            "rd_OR": 1.0, "rd_p": 1.0, "rd_fdr": 1.0,
            "n_genes_expressed": 0, "n_high_expr": 0,
        })
        continue
    
    # Compute mean RPKM per gene for this stage
    gene_rpkm = {}
    for eid, expr_vals in gene_expr.items():
        stage_vals = [expr_vals[i] for i in stage_indices if i < len(expr_vals)]
        if stage_vals:
            gene_rpkm[eid] = float(np.mean(stage_vals))
    
    if len(gene_rpkm) < 100:
        continue
    
    # High-expression genes (top 25%)
    rpkm_vals = np.array(list(gene_rpkm.values()))
    non_zero = rpkm_vals[rpkm_vals > 0]
    if len(non_zero) < 50:
        continue
    threshold = np.percentile(non_zero, 75)
    high_expr = set(g for g, v in gene_rpkm.items() if v >= threshold and v > 0)
    
    # Fisher tests
    gd_odds, gd_p, gd_n, _ = fisher_test(gd_genes, all_genes, high_expr)
    rd_odds, rd_p, rd_n, _ = fisher_test(rd_genes, all_genes, high_expr)
    
    stage_results.append({
        "stage": age,
        "n_samples": len(stage_indices),
        "n_genes_expressed": len(gene_rpkm),
        "n_high_expr": len(high_expr),
        "gd_OR": gd_odds,
        "gd_p": gd_p,
        "gd_n_high": gd_n,
        "rd_OR": rd_odds,
        "rd_p": rd_p,
        "rd_n_high": rd_n,
    })

# BH FDR
if stage_results:
    gd_pvals = [r["gd_p"] for r in stage_results]
    rd_pvals = [r["rd_p"] for r in stage_results]
    gd_fdrs = bh_fdr(gd_pvals)
    rd_fdrs = bh_fdr(rd_pvals)
    for i, r in enumerate(stage_results):
        r["gd_fdr"] = gd_fdrs[i]
        r["rd_fdr"] = rd_fdrs[i]

print(f"\n  Per-stage enrichment ({len(stage_results)} stages):")
print(f"  {'Stage':15s} {'N':5s} {'GD OR':8s} {'GD FDR':10s} {'RD OR':8s} {'RD FDR':10s}")
for r in stage_results:
    gd_s = "*" if r.get("gd_fdr", 1) < 0.05 else " "
    rd_s = "*" if r.get("rd_fdr", 1) < 0.05 else " "
    print(f"  {r['stage']:15s} {r['n_samples']:5d} "
          f"{r['gd_OR']:8.2f}{gd_s} {r.get('gd_fdr', 1):10.2e} "
          f"{r['rd_OR']:8.2f}{rd_s} {r.get('rd_fdr', 1):10.2e}")

# Write stage results
stage_csv = OUTPUT_DIR / "developmental_full_trajectory.csv"
with open(stage_csv, "w", newline="") as f:
    if stage_results:
        writer = csv.DictWriter(f, fieldnames=list(stage_results[0].keys()))
        writer.writeheader()
        for r in stage_results:
            writer.writerow(r)
print(f"\n  Output: {stage_csv}")

# ============================================
# Analysis 2: Developmental Tau
# ============================================
print(f"\n{'='*50}")
print("Analysis 2: Developmental Tau (cross-stage specificity)")
print(f"{'='*50}")

# For each gene, compute mean expression per stage
gene_stage_expr = {}  # gene_id -> [mean RPKM per stage]

for eid, expr_vals in gene_expr.items():
    stage_means = []
    for age in unique_ages:
        stage_indices = [i for i, s in enumerate(sample_info) if s["age"] == age]
        if stage_indices:
            vals = [expr_vals[i] for i in stage_indices if i < len(expr_vals)]
            stage_means.append(float(np.mean(vals)) if vals else 0.0)
        else:
            stage_means.append(0.0)
    gene_stage_expr[eid] = stage_means

# Compute developmental tau
gene_dev_tau = {}
for eid, stage_vals in gene_stage_expr.items():
    gene_dev_tau[eid] = compute_tau(stage_vals)

# Compare tau between groups
gd_tau = [gene_dev_tau[g] for g in gd_genes if g in gene_dev_tau]
rd_tau = [gene_dev_tau[g] for g in rd_genes if g in gene_dev_tau]
neutral_tau = [gene_dev_tau[g] for g in neutral_genes if g in gene_dev_tau]

print(f"\n  Developmental tau (cross {len(unique_ages)} stages):")
print(f"    GD: median={np.median(gd_tau):.4f}, mean={np.mean(gd_tau):.4f} (n={len(gd_tau)})")
print(f"    RD: median={np.median(rd_tau):.4f}, mean={np.mean(rd_tau):.4f} (n={len(rd_tau)})")
print(f"    Neutral: median={np.median(neutral_tau):.4f}, mean={np.mean(neutral_tau):.4f} (n={len(neutral_tau)})")

# Mann-Whitney U tests
if gd_tau and rd_tau:
    u, p = stats.mannwhitneyu(rd_tau, gd_tau, alternative="greater")
    print(f"    RD > GD: U={u:.0f}, p={p:.2e}")
    u2, p2 = stats.mannwhitneyu(gd_tau, rd_tau, alternative="greater")
    print(f"    GD > RD: U={u2:.0f}, p={p2:.2e}")
    
    # Cliff's delta
    n_total = len(gd_tau) * len(rd_tau)
    n_greater = sum(1 for r in rd_tau for g in gd_tau if r > g)
    n_less = sum(1 for r in rd_tau for g in gd_tau if r < g)
    cliff_d = (n_greater - n_less) / n_total if n_total > 0 else 0
    print(f"    Cliff's delta (RD-GD): {cliff_d:.4f}")

# Kruskal-Wallis across all groups
if gd_tau and rd_tau and neutral_tau:
    h, p_kw = stats.kruskal(gd_tau, rd_tau, neutral_tau)
    print(f"    Kruskal-Wallis: H={h:.2f}, p={p_kw:.2e}")

# Write developmental tau
dev_tau_csv = OUTPUT_DIR / "developmental_tau_comparison.csv"
with open(dev_tau_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gene_id", "gene_symbol", "classification_loo_brain", "dev_tau", "n_stages"])
    for gid in sorted(gene_dev_tau.keys()):
        sym = genes[gid].get("gene_symbol", "")
        cls = genes[gid].get("classification_loo_brain", "")
        writer.writerow([gid, sym, cls, gene_dev_tau[gid], len(unique_ages)])
print(f"\n  Output: {dev_tau_csv}")

# ============================================
# Analysis 3: Prenatal vs Postnatal
# ============================================
print(f"\n{'='*50}")
print("Analysis 3: Prenatal vs Postnatal Expression")
print(f"{'='*50}")

prenatal_indices = [i for i, s in enumerate(sample_info) if s["period"] in ("early_prenatal", "late_prenatal")]
postnatal_indices = [i for i, s in enumerate(sample_info) if s["period"] in ("early_postnatal", "late_postnatal")]

print(f"  Prenatal samples: {len(prenatal_indices)}")
print(f"  Postnatal samples: {len(postnatal_indices)}")

gene_prenatal = {}
gene_postnatal = {}
for eid, expr_vals in gene_expr.items():
    if prenatal_indices:
        vals = [expr_vals[i] for i in prenatal_indices if i < len(expr_vals)]
        gene_prenatal[eid] = float(np.mean(vals)) if vals else 0
    if postnatal_indices:
        vals = [expr_vals[i] for i in postnatal_indices if i < len(expr_vals)]
        gene_postnatal[eid] = float(np.mean(vals)) if vals else 0

# Prenatal/postnatal ratio
# NOTE: fixed 2026-09-14 — the previous version filled ratio=1.0 when both means
# were <0.1 (~887 genes, uninformative fill) and inf when prenatal>=0.1 &
# postnatal<0.1 (78 genes), silently distorting downstream statistics.
# ratio is now the pure prenatal_mean/postnatal_mean; NaN only when both are 0.
# See scripts/pB_fig4d_adjudicate.py and results/paperB/_verify/fig4d_adjudication.txt.
gene_ratio = {}
for eid in gene_prenatal:
    pre = gene_prenatal[eid]
    post = gene_postnatal.get(eid, 0)
    if post > 0:
        gene_ratio[eid] = pre / post
    elif pre > 0:
        gene_ratio[eid] = float('inf')
    else:
        gene_ratio[eid] = float('nan')

gd_ratios = [gene_ratio[g] for g in gd_genes if g in gene_ratio and np.isfinite(gene_ratio[g])]
rd_ratios = [gene_ratio[g] for g in rd_genes if g in gene_ratio and np.isfinite(gene_ratio[g])]

if gd_ratios and rd_ratios:
    u_pre, p_pre = stats.mannwhitneyu(rd_ratios, gd_ratios, alternative="greater")
    print(f"  Prenatal/Postnatal ratio RD > GD: U={u_pre:.0f}, p={p_pre:.2e}")
    print(f"    GD median ratio: {np.median(gd_ratios):.4f} (n={len(gd_ratios)})")
    print(f"    RD median ratio: {np.median(rd_ratios):.4f} (n={len(rd_ratios)})")
    
    u_gd, p_gd = stats.mannwhitneyu(gd_ratios, rd_ratios, alternative="greater")
    print(f"  Prenatal/Postnatal ratio GD > RD: U={u_gd:.0f}, p={p_gd:.2e}")

# Write prenatal/postnatal
prepost_csv = OUTPUT_DIR / "prenatal_vs_postnatal.csv"
with open(prepost_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gene_id", "gene_symbol", "classification_loo_brain", 
                      "prenatal_mean", "postnatal_mean", "ratio"])
    for eid in sorted(gene_prenatal.keys()):
        sym = genes[eid].get("gene_symbol", "")
        cls = genes[eid].get("classification_loo_brain", "")
        pre = gene_prenatal[eid]
        post = gene_postnatal.get(eid, 0)
        r = gene_ratio[eid]
        writer.writerow([eid, sym, cls, pre, post, r])
print(f"\n  Output: {prepost_csv}")

# ============================================
# Analysis 4: Region × Stage Matrix
# ============================================
print(f"\n{'='*50}")
print("Analysis 4: Region × Stage Enrichment Matrix")
print(f"{'='*50}")

# For each region × stage combination: compute GD/RD enrichment
region_stage_results = []

for region in unique_regions:
    for age in unique_ages:
        rs_indices = [i for i, s in enumerate(sample_info) 
                       if s["region"] == region and s["age"] == age]
        
        if len(rs_indices) < 2:
            continue
        
        # Compute mean RPKM per gene
        gene_rpkm = {}
        for eid, expr_vals in gene_expr.items():
            rs_vals = [expr_vals[i] for i in rs_indices if i < len(expr_vals)]
            if rs_vals:
                gene_rpkm[eid] = float(np.mean(rs_vals))
        
        if len(gene_rpkm) < 50:
            continue
        
        rpkm_vals = np.array(list(gene_rpkm.values()))
        non_zero = rpkm_vals[rpkm_vals > 0]
        if len(non_zero) < 20:
            continue
        threshold = np.percentile(non_zero, 75)
        high_expr = set(g for g, v in gene_rpkm.items() if v >= threshold and v > 0)
        
        gd_odds, gd_p, _, _ = fisher_test(gd_genes, all_genes, high_expr)
        rd_odds, rd_p, _, _ = fisher_test(rd_genes, all_genes, high_expr)
        
        region_stage_results.append({
            "region": region,
            "stage": age,
            "n_samples": len(rs_indices),
            "n_genes": len(gene_rpkm),
            "gd_OR": gd_odds,
            "gd_p": gd_p,
            "rd_OR": rd_odds,
            "rd_p": rd_p,
        })

print(f"  Region × Stage combinations: {len(region_stage_results)}")

# Write region × stage matrix
rs_csv = OUTPUT_DIR / "region_x_stage_heatmap.csv"
with open(rs_csv, "w", newline="") as f:
    if region_stage_results:
        writer = csv.DictWriter(f, fieldnames=list(region_stage_results[0].keys()))
        writer.writeheader()
        for r in region_stage_results:
            writer.writerow(r)
print(f"  Output: {rs_csv}")

# ============================================
# Analysis 5: Developmental Period Summary
# ============================================
print(f"\n{'='*50}")
print("Analysis 5: Developmental Period Summary")
print(f"{'='*50}")

period_results = []
for period in ["early_prenatal", "late_prenatal", "early_postnatal", "late_postnatal"]:
    period_indices = [i for i, s in enumerate(sample_info) if s["period"] == period]
    
    if len(period_indices) < 3:
        continue
    
    gene_rpkm = {}
    for eid, expr_vals in gene_expr.items():
        vals = [expr_vals[i] for i in period_indices if i < len(expr_vals)]
        if vals:
            gene_rpkm[eid] = float(np.mean(vals))
    
    if len(gene_rpkm) < 100:
        continue
    
    rpkm_vals = np.array(list(gene_rpkm.values()))
    non_zero = rpkm_vals[rpkm_vals > 0]
    threshold = np.percentile(non_zero, 75)
    high_expr = set(g for g, v in gene_rpkm.items() if v >= threshold and v > 0)
    
    gd_odds, gd_p, gd_n, _ = fisher_test(gd_genes, all_genes, high_expr)
    rd_odds, rd_p, rd_n, _ = fisher_test(rd_genes, all_genes, high_expr)
    
    period_results.append({
        "period": period,
        "n_samples": len(period_indices),
        "gd_OR": gd_odds,
        "gd_p": gd_p,
        "gd_n_high": gd_n,
        "rd_OR": rd_odds,
        "rd_p": rd_p,
        "rd_n_high": rd_n,
    })

print(f"\n  Developmental period results:")
for r in period_results:
    print(f"    {r['period']:20s} N={r['n_samples']:4d} "
          f"GD:OR={r['gd_OR']:.2f},p={r['gd_p']:.2e} "
          f"RD:OR={r['rd_OR']:.2f},p={r['rd_p']:.2e}")

dev_period_csv = OUTPUT_DIR / "developmental_period_enrichment.csv"
with open(dev_period_csv, "w", newline="") as f:
    if period_results:
        writer = csv.DictWriter(f, fieldnames=list(period_results[0].keys()))
        writer.writeheader()
        for r in period_results:
            writer.writerow(r)
print(f"\n  Output: {dev_period_csv}")

# ============================================
# Summary
# ============================================
summary = {
    "n_stages": len(unique_ages),
    "n_regions": len(unique_regions),
    "n_genes_matched": len(gene_expr),
    "n_region_stage_combos": len(region_stage_results),
    "gd_dev_tau_median": float(np.median(gd_tau)) if gd_tau else None,
    "rd_dev_tau_median": float(np.median(rd_tau)) if rd_tau else None,
    "dev_tau_rd_vs_gd_p": float(p) if gd_tau and rd_tau else None,
    "prenatal_postnatal_rd_vs_gd_p": float(p_pre) if gd_ratios and rd_ratios else None,
    "gd_fdr_sig_stages": [r["stage"] for r in stage_results if r.get("gd_fdr", 1) < 0.05],
    "rd_fdr_sig_stages": [r["stage"] for r in stage_results if r.get("rd_fdr", 1) < 0.05],
}

summary_file = OUTPUT_DIR / "developmental_analysis_summary.json"
with open(summary_file, "w") as f:
    json.dump(summary, f, indent=2)
print(f"\n  Summary: {summary_file}")

print("\n" + "=" * 70)
print("Full Developmental Trajectory Analysis Complete!")
print("=" * 70)
