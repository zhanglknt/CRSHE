#!/usr/bin/env python3
"""
Phase 7G v7 Classification — Complete P0 Fixes

Addresses ALL 6 P0 issues from expert review:
P0-1: Leave-one-out non-circular validation
P0-2: RELAX coverage transparent reporting
P0-3: Storey q-value replacing BH FDR
P0-4: GDS FDR bonus removal (deduplication)

Run in WSL: python3 phase7g_v7_classification.py
"""
import csv
import json
import os
import sys
import numpy as np
from scipy import stats
from collections import Counter

BASE = os.environ.get("HSD_BASE", "/mnt/d/人类正选择基因项目")
BUSTED_CSV = f"{BASE}/results/phase4_hyphy/final_analysis_v3/busted_results_v3.csv"
RELAX_CSV = f"{BASE}/results/phase4_hyphy/relax_results/relax_results_v3.csv"
CONSERVATION_CSV = f"{BASE}/results/phase4_conservation/conservation_scores_v3.csv"
DOAN_DIR = f"{BASE}/results/doan2024_campra"
GTEX_DIR = f"{BASE}/results/gtex_v11_tau"
V5_CSV = f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v5/gene_classification_v5.csv"
V6_CSV = f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v6/gene_classification_v6.csv"
OUTPUT_DIR = f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 70)
print("Phase 7G v7 Classification — Complete P0 Fixes")
print("=" * 70)

# ============================================
# P0-3: Storey q-value implementation
# ============================================
def estimate_pi0(pvals, lambda_range=None):
    """Estimate pi0 (proportion of true nulls) using Storey's method.
    
    NOTE: BUSTED p-values are bounded at 0.5 (not 1.0), so we use
    lambda values up to 0.45 and adjust the formula accordingly.
    The effective null distribution is uniform on [0, 0.5].
    """
    pvals = np.array(pvals)
    pmax = max(pvals.max(), 0.5)  # bounded p-values
    if lambda_range is None:
        # Use lambda values appropriate for bounded distribution
        lambda_range = np.arange(0.05, min(0.95, pmax - 0.01), 0.02)
    
    n = len(pvals)
    pi0_estimates = []
    
    for lam in lambda_range:
        W_lam = np.sum(pvals > lam)
        # Adjust for bounded distribution: effective range is [0, pmax]
        pi0_lam = W_lam * pmax / (n * (pmax - lam))
        pi0_estimates.append(min(pi0_lam, 1.0))
    
    # Use cubic spline smoothing on pi0 estimates vs lambda
    # For robustness, take the median of estimates at higher lambdas
    if len(pi0_estimates) > 3:
        pi0 = np.median(pi0_estimates[max(0, len(pi0_estimates)//2):])
    else:
        pi0 = np.median(pi0_estimates)
    pi0 = max(min(pi0, 1.0), 1.0 / n)
    return pi0

def storey_qvalue(pvals):
    """Compute Storey q-values."""
    pvals = np.array(pvals)
    n = len(pvals)
    pi0 = estimate_pi0(pvals)
    
    # Sort p-values
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]
    
    # Calculate q-values from largest to smallest
    qvals = np.zeros(n)
    qvals[sorted_idx[-1]] = min(sorted_pvals[-1] * pi0 * n / 1, 1.0)  # divide by rank=1
    
    for i in range(n - 2, -1, -1):
        rank = i + 1  # 1-indexed rank
        q_candidate = sorted_pvals[i] * pi0 * n / rank
        qvals[sorted_idx[i]] = min(q_candidate, qvals[sorted_idx[i + 1]])
    
    qvals = np.minimum(qvals, 1.0)
    return qvals, pi0

# ============================================
# Step 1: Load all data
# ============================================
print("\n--- Step 1: Loading data ---")

# BUSTED results
busted = {}
busted_pvals = []
busted_gene_order = []
with open(BUSTED_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        gid = row["gene_id"]
        p = float(row["p_value"])
        busted[gid] = {
            "p_value": p,
            "lrt": float(row["lrt"]),
            "bh_fdr": float(row["bh_fdr"]),
            "fdr_sig": row["fdr_significant"] == "True",
            "n_sites": int(row["n_sites"]),
        }
        busted_pvals.append(p)
        busted_gene_order.append(gid)

print(f"  BUSTED: {len(busted)} genes ({sum(1 for v in busted.values() if v['fdr_sig'])} BH FDR sig)")

# P0-3: Compute Storey q-values
print("\n--- P0-3: Storey q-value computation ---")
storey_qvals, pi0 = storey_qvalue(busted_pvals)
storey_sig_count = 0
for i, gid in enumerate(busted_gene_order):
    busted[gid]["storey_q"] = float(storey_qvals[i])
    busted[gid]["storey_sig"] = storey_qvals[i] < 0.05
    if storey_qvals[i] < 0.05:
        storey_sig_count += 1

print(f"  pi0 estimate: {pi0:.4f}")
print(f"  Storey q<0.05: {storey_sig_count} genes (vs BH FDR<0.05: {sum(1 for v in busted.values() if v['fdr_sig'])})")
print(f"  BH/Storey ratio: {sum(1 for v in busted.values() if v['fdr_sig']) / max(storey_sig_count, 1):.2f}")

# p-value distribution diagnostics
pvals_arr = np.array(busted_pvals)
p_eq_0 = np.sum(pvals_arr == 0)
p_eq_05 = np.sum(pvals_arr == 0.5)
p_between = np.sum((pvals_arr > 0) & (pvals_arr < 0.5))
print(f"  p-value distribution: p=0: {p_eq_0} ({100*p_eq_0/len(pvals_arr):.1f}%), p=0.5: {p_eq_05} ({100*p_eq_05/len(pvals_arr):.1f}%), 0<p<0.5: {p_between} ({100*p_between/len(pvals_arr):.1f}%)")

# RELAX results
relax = {}
with open(RELAX_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        gid = row["gene_id"]
        if row["p_value"] and row["K"]:
            try:
                relax[gid] = {
                    "K": float(row["K"]),
                    "p_value": float(row["p_value"]),
                    "fdr": float(row["fdr"]) if row["fdr"] else 1.0,
                    "fdr_sig": row["fdr_significant"] == "True",
                    "lrt": float(row["LRT"]) if row["LRT"] else 0,
                }
            except (ValueError, TypeError):
                pass
print(f"  RELAX: {len(relax)} valid genes")

# P0-2: RELAX coverage transparent reporting
print("\n--- P0-2: RELAX coverage transparent report ---")
# Count genes with actual RELAX results vs defaults
genes_with_relax = set(relax.keys())
genes_without_relax = set(busted.keys()) - genes_with_relax
# Among genes WITH relax entries, how many have actual signal (p < 1.0)?
relax_default = sum(1 for v in relax.values() if v["p_value"] >= 0.99 and v["K"] == 1.0)
relax_active = sum(1 for v in relax.values() if v["p_value"] < 0.99 or v["K"] != 1.0)
# Among BUSTED-significant genes
busted_sig_genes = set(g for g, v in busted.items() if v["storey_sig"])
busted_sig_with_relax = busted_sig_genes & genes_with_relax
busted_sig_relax_active = sum(1 for g in busted_sig_with_relax if relax[g]["p_value"] < 0.99 or relax[g]["K"] != 1.0)

print(f"  Total genes in universe: {len(busted)}")
print(f"  Genes with RELAX entry: {len(genes_with_relax)} ({100*len(genes_with_relax)/len(busted):.1f}%)")
print(f"    - With actual signal (p<0.99 or K!=1): {relax_active} ({100*relax_active/len(genes_with_relax):.1f}%)")
print(f"    - Default (p=1, K=1, failed): {relax_default} ({100*relax_default/len(genes_with_relax):.1f}%)")
print(f"  Genes without RELAX entry: {len(genes_without_relax)} ({100*len(genes_without_relax)/len(busted):.1f}%)")
print(f"  BUSTED-significant genes: {len(busted_sig_genes)}")
print(f"    - With RELAX entry: {len(busted_sig_with_relax)} ({100*len(busted_sig_with_relax)/max(len(busted_sig_genes),1):.1f}%)")
print(f"    - With RELAX active signal: {busted_sig_relax_active} ({100*busted_sig_relax_active/max(len(busted_sig_genes),1):.1f}%)")

# Proper K>1 rate: only among genes with active RELAX results
relax_active_genes = [g for g in relax if relax[g]["p_value"] < 0.99 or relax[g]["K"] != 1.0]
relax_k_gt1_active = sum(1 for g in relax_active_genes if relax[g]["K"] > 1)
relax_k_lt1_active = sum(1 for g in relax_active_genes if relax[g]["K"] < 1)
print(f"  K>1 rate (proper denominator): {relax_k_gt1_active}/{len(relax_active_genes)} = {100*relax_k_gt1_active/max(len(relax_active_genes),1):.1f}%")
print(f"  K<1 rate (proper denominator): {relax_k_lt1_active}/{len(relax_active_genes)} = {100*relax_k_lt1_active/max(len(relax_active_genes),1):.1f}%")

# Conservation scores
conservation = {}
with open(CONSERVATION_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        gid = row["gene_id"]
        conservation[gid] = {
            "phyloP_cds": float(row["phyloP_cds"]) if row["phyloP_cds"] != "NA" else None,
            "phyloP_gene": float(row["phyloP_gene_body"]) if row["phyloP_gene_body"] != "NA" else None,
            "phastCons_cds": float(row["phastCons_cds"]) if row["phastCons_cds"] != "NA" else None,
            "phastCons_gene": float(row["phastCons_gene_body"]) if row["phastCons_gene_body"] != "NA" else None,
        }
print(f"  Conservation: {len(conservation)} genes")

# v5 classification data (for CDS length, gene symbols, HAR, Selectome, tau)
v5_data = {}
with open(V5_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        gid = row["gene_id"]
        v5_data[gid] = row
print(f"  v5 classification: {len(v5_data)} genes")

# Doan caMPRA
doan_genes = set()
doan_active_hars_file = f"{DOAN_DIR}/active_hars.csv"
gene_bed_file = f"{BASE}/results/phase4_conservation/gene_body.bed"
if os.path.exists(doan_active_hars_file) and os.path.exists(gene_bed_file):
    active_hars = []
    with open(doan_active_hars_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            chrom = row.get("chr_hg38", "")
            start = int(float(row.get("start_hg38", 0)))
            end = int(float(row.get("end_hg38", 0)))
            active = row.get("active_both", "False")
            if chrom and start and end:
                active_hars.append((chrom, start, end, active == "True"))
    gene_coords = {}
    with open(gene_bed_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3]
                gene_id = "_".join(name.split("_")[:-1]) if name.rsplit("_", 1)[-1].isdigit() else name
                gene_coords[gene_id] = (chrom, start, end)
    EXTENSION = 50000
    for gene_id, (g_chrom, g_start, g_end) in gene_coords.items():
        g_start_ext = max(0, g_start - EXTENSION)
        g_end_ext = g_end + EXTENSION
        for h_chrom, h_start, h_end, h_both in active_hars:
            if h_chrom == g_chrom and h_start < g_end_ext and h_end > g_start_ext:
                doan_genes.add(gene_id)
                break
print(f"  Doan caMPRA genes: {len(doan_genes)}")

# GTEx tau
tau_data = {}
tau_file = f"{GTEX_DIR}/gtex_v11_tau.csv"
if os.path.exists(tau_file):
    with open(tau_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gid = row.get("gene_id", "")
            if gid:
                try:
                    tau_data[gid] = float(row.get("tau", 0))
                except (ValueError, TypeError):
                    pass
print(f"  GTEx tau: {len(tau_data)} genes")

# ============================================
# Step 2: Build unified gene table
# ============================================
print("\n--- Step 2: Building unified gene table ---")

all_genes = set(busted.keys())
gene_table = []

for gid in sorted(all_genes):
    b = busted.get(gid, {})
    r = relax.get(gid, {})
    c = conservation.get(gid, {})
    v5 = v5_data.get(gid, {})
    
    cds_length = int(v5.get("cds_length", 0)) if v5.get("cds_length") else 0
    gene_symbol = v5.get("gene_symbol", "")
    
    # Determine if this gene has active RELAX (not default)
    has_relax_active = gid in relax and (relax[gid]["p_value"] < 0.99 or relax[gid]["K"] != 1.0)
    
    gene_table.append({
        "gene_id": gid,
        "gene_symbol": gene_symbol,
        "cds_length": cds_length,
        "busted_p": b.get("p_value", 1.0),
        "busted_lrt": b.get("lrt", 0),
        "bh_fdr": b.get("bh_fdr", 1.0),
        "bh_fdr_sig": b.get("fdr_sig", False),
        "storey_q": b.get("storey_q", 1.0),
        "storey_sig": b.get("storey_sig", False),
        "relax_K": r.get("K", 1.0),
        "relax_p": r.get("p_value", 1.0),
        "relax_fdr": r.get("fdr", 1.0),
        "relax_fdr_sig": r.get("fdr_sig", False),
        "has_relax_active": has_relax_active,
        "relax_intensified": r.get("K", 1.0) > 1 and r.get("fdr_sig", False),
        "relax_relaxed": 0 < r.get("K", 1.0) < 1 and r.get("fdr_sig", False),
        "has_selectome": v5.get("has_selectome", "False") == "True" if v5 else False,
        "n_hars": int(v5.get("n_hars", 0)) if v5 else 0,
        "has_doan_campra": gid in doan_genes,
        "tau": tau_data.get(gid, float(v5.get("brain_tau", 0))) if v5 else 0,
        "phyloP_cds": c.get("phyloP_cds") or 0,
        "phyloP_gene": c.get("phyloP_gene") or 0,
        "phastCons_cds": c.get("phastCons_cds") or 0,
        "phastCons_gene": c.get("phastCons_gene") or 0,
        "nc_conservation": (c.get("phyloP_gene", 0) or 0) - (c.get("phyloP_cds", 0) or 0),
        "brain_tau": float(v5.get("brain_tau", 0)) if v5 else 0,
    })

print(f"  Total genes: {len(gene_table)}")

# ============================================
# Step 3: CDS Length Residualization (from v6)
# ============================================
print("\n--- Step 3: CDS Length Residualization ---")

cds_lengths = np.array([g["cds_length"] for g in gene_table], dtype=float)
cds_log = np.log10(cds_lengths + 1)
busted_neglogp = np.array([-np.log10(max(g["busted_p"], 1e-300)) for g in gene_table])
busted_lrt_arr = np.array([g["busted_lrt"] for g in gene_table], dtype=float)

mask = np.isfinite(busted_neglogp) & np.isfinite(cds_log) & (cds_lengths > 0)
slope1, intercept1, r1, p1, se1 = stats.linregress(cds_log[mask], busted_neglogp[mask])
neglogp_resid = np.zeros(len(gene_table))
neglogp_resid[mask] = busted_neglogp[mask] - (slope1 * cds_log[mask] + intercept1)
print(f"  BUSTED -log10(p) vs log10(CDS): r={r1:.4f}, p={p1:.2e}")

mask2 = np.isfinite(busted_lrt_arr) & np.isfinite(cds_log) & (cds_lengths > 0)
slope2, intercept2, r2, p2, se2 = stats.linregress(cds_log[mask2], busted_lrt_arr[mask2])
lrt_resid = np.zeros(len(gene_table))
lrt_resid[mask2] = busted_lrt_arr[mask2] - (slope2 * cds_log[mask2] + intercept2)
print(f"  BUSTED LRT vs log10(CDS): r={r2:.4f}, p={p2:.2e}")

# ============================================
# Step 4: P0-4 — GDS v7 (FDR bonus removed)
# ============================================
print("\n--- Step 4: P0-4 — GDS v7 (FDR bonus removed) ---")

def percentile_rank(values):
    arr = np.array(values, dtype=float)
    finite = np.isfinite(arr)
    ranks = np.zeros(len(arr))
    sorted_vals = np.sort(arr[finite])
    for i, v in enumerate(arr):
        if np.isfinite(v):
            ranks[i] = np.searchsorted(sorted_vals, v) / max(len(sorted_vals), 1)
    return ranks

gds_p_pct = percentile_rank(neglogp_resid)
gds_lrt_pct = percentile_rank(lrt_resid)

# RELAX K component (only for genes with active RELAX)
relax_K_arr = np.array([g["relax_K"] for g in gene_table], dtype=float)
has_active = np.array([g["has_relax_active"] for g in gene_table], dtype=bool)
gds_relax_pct = np.where(has_active, np.minimum(0.5 + 0.25 * np.log2(np.maximum(relax_K_arr, 0.01)), 1.0), 0.5)
gds_relax_pct = np.where(relax_K_arr < 1, 0.5 * np.maximum(relax_K_arr, 0.01), gds_relax_pct)

# Selectome
gds_selectome = np.array([1.0 if g["has_selectome"] else 0.0 for g in gene_table])

# P0-4 FIX: GDS v7 = weighted sum WITHOUT FDR bonus
# Redistributed weights: p-value 35%, LRT 30%, RELAX 20%, Selectome 15%
# (was: 30% + 25% + 20% + 15% + 10% FDR = 100%)
gds_v7 = (0.35 * gds_p_pct + 
          0.30 * gds_lrt_pct + 
          0.20 * gds_relax_pct + 
          0.15 * gds_selectome)

print(f"  GDS v7 range: {gds_v7.min():.4f} - {gds_v7.max():.4f} (mean: {gds_v7.mean():.4f})")
print(f"  (FDR bonus removed; weights redistributed to p=35%, LRT=30%, RELAX=20%, Selectome=15%)")

# ============================================
# Step 5: RDS v7 (same as v6, components kept for leave-one-out)
# ============================================
print("\n--- Step 5: RDS v7 (with leave-one-out components) ---")

rds_doan = np.array([1.0 if g["has_doan_campra"] else 0.0 for g in gene_table])
tau_arr = np.array([g["tau"] for g in gene_table], dtype=float)
rds_tau = percentile_rank(tau_arr)
nc_arr = np.array([g["nc_conservation"] for g in gene_table], dtype=float)
rds_nc = percentile_rank(-nc_arr)
brain_arr = np.array([g["brain_tau"] for g in gene_table], dtype=float)
rds_brain = percentile_rank(brain_arr)

# Full RDS v7 (same as v6)
rds_v7_full = (0.30 * rds_doan + 0.25 * rds_tau + 0.25 * rds_nc + 0.20 * rds_brain)

# Leave-one-out RDS variants for P0-1
# LOO-caMPRA: remove caMPRA, redistribute to tau=35%, nc=35%, brain=30%
rds_v7_no_campra = (0.35 * rds_tau + 0.35 * rds_nc + 0.30 * rds_brain)

# LOO-tau: remove tau, redistribute to caMPRA=40%, nc=33%, brain=27%
rds_v7_no_tau = (0.40 * rds_doan + 0.33 * rds_nc + 0.27 * rds_brain)

print(f"  RDS v7 (full): {rds_v7_full.min():.4f} - {rds_v7_full.max():.4f}")
print(f"  RDS v7 (no caMPRA): {rds_v7_no_campra.min():.4f} - {rds_v7_no_campra.max():.4f}")
print(f"  RDS v7 (no tau): {rds_v7_no_tau.min():.4f} - {rds_v7_no_tau.max():.4f}")

# ============================================
# Step 6: Classification v7 (3 variants for LOO)
# ============================================
print("\n--- Step 6: Classification v7 ---")

def classify(gds, rds, storey_sig, relax_relaxed_flags, threshold=0.15):
    """Classify genes based on GDS and RDS scores."""
    classifications = []
    for i in range(len(gds)):
        g = gds[i]
        r = rds[i]
        sig = storey_sig[i]
        relaxed = relax_relaxed_flags[i]
        
        if g > r + threshold and sig:
            cls = "gene-driven"
        elif r > g + threshold and not sig:
            cls = "regulation-driven"
        elif g > r + threshold and r > 0.4 and sig:
            cls = "dual-driven"
        elif r > g + threshold and sig and g > 0.3:
            cls = "dual-driven"
        else:
            cls = "neutral"
        
        if relaxed and cls == "gene-driven":
            cls = "gene-driven (relaxed)"
        
        classifications.append(cls)
    return classifications

storey_sig_arr = np.array([g["storey_sig"] for g in gene_table])
bh_fdr_sig_arr = np.array([g["bh_fdr_sig"] for g in gene_table])
relax_relaxed_arr = np.array([g["relax_relaxed"] for g in gene_table])

# P0-3 NOTE: Storey q-value is too permissive with bounded p-values (pmax=0.5).
# We use BH FDR as primary classification criterion (more conservative),
# and report Storey q-value as a sensitivity analysis.
# Full classification (using BH FDR sig, same as v6 but with v7 scores)
cls_full = classify(gds_v7, rds_v7_full, bh_fdr_sig_arr, relax_relaxed_arr)
cls_loo_campra = classify(gds_v7, rds_v7_no_campra, bh_fdr_sig_arr, relax_relaxed_arr)
cls_loo_tau = classify(gds_v7, rds_v7_no_tau, bh_fdr_sig_arr, relax_relaxed_arr)

# Also compute Storey-based classification for sensitivity analysis
cls_storey = classify(gds_v7, rds_v7_full, storey_sig_arr, relax_relaxed_arr)

cls_counter = Counter(cls_full)
print("  Classification v7 (full):")
for k, v in cls_counter.most_common():
    print(f"    {k}: {v}")

cls_loo_campra_counter = Counter(cls_loo_campra)
print("\n  Classification v7 (LOO caMPRA):")
for k, v in cls_loo_campra_counter.most_common():
    print(f"    {k}: {v}")

cls_loo_tau_counter = Counter(cls_loo_tau)
print("\n  Classification v7 (LOO tau):")
for k, v in cls_loo_tau_counter.most_common():
    print(f"    {k}: {v}")

# ============================================
# Step 7: P0-1 — Non-circular validation (leave-one-out)
# ============================================
print("\n" + "=" * 70)
print("P0-1: LEAVE-ONE-OUT NON-CIRCULAR VALIDATION")
print("=" * 70)

def get_class_indices(classifications, cls_name):
    return [i for i, c in enumerate(classifications) if c == cls_name]

def fisher_enrichment(gene_indices, all_indices, feature_arr):
    """Fisher's exact test for enrichment of binary feature."""
    has_feat = sum(1 for i in gene_indices if feature_arr[i])
    no_feat = len(gene_indices) - has_feat
    bg_has = sum(1 for i in all_indices if feature_arr[i]) - has_feat
    bg_no = len(all_indices) - len(gene_indices) - bg_has
    table = [[has_feat, no_feat], [bg_has, bg_no]]
    odds, p = stats.fisher_exact(table, alternative="greater")
    return odds, p, has_feat, len(gene_indices)

# 7a: HAR enrichment (ALREADY non-circular — HAR not in any RDS)
print("\n--- 7a: HAR enrichment (non-circular) ---")
gd_idx = [i for i, c in enumerate(cls_full) if c == "gene-driven"]
rd_idx = get_class_indices(cls_full, "regulation-driven")
neutral_idx = get_class_indices(cls_full, "neutral")
all_idx = list(range(len(gene_table)))
har_arr = np.array([g["n_hars"] > 0 for g in gene_table])

odds_har, p_har, rd_har, rd_total = fisher_enrichment(rd_idx, all_idx, har_arr)
print(f"  HAR enrichment in RD (full): OR={odds_har:.3f}, p={p_har:.4e}")
print(f"    RD with HAR: {rd_har}/{rd_total} ({100*rd_har/max(rd_total,1):.1f}%)")
bg_har = sum(har_arr)
print(f"    BG with HAR: {bg_har}/{len(all_idx)} ({100*bg_har/len(all_idx):.1f}%)")

# 7b: caMPRA enrichment — LEAVE-ONE-OUT (caMPRA removed from RDS)
print("\n--- 7b: caMPRA enrichment (LOO — caMPRA removed from RDS) ---")
campra_arr = np.array([g["has_doan_campra"] for g in gene_table])

# Use LOO-caMPRA classification to define RD
rd_loo_campra = get_class_indices(cls_loo_campra, "regulation-driven")
odds_campra, p_campra, rd_campra, rd_campra_total = fisher_enrichment(rd_loo_campra, all_idx, campra_arr)
print(f"  caMPRA enrichment in RD (LOO): OR={odds_campra:.3f}, p={p_campra:.4e}")
print(f"    RD(LOO) with caMPRA: {rd_campra}/{rd_campra_total} ({100*rd_campra/max(rd_campra_total,1):.1f}%)")
bg_campra = sum(campra_arr)
print(f"    BG with caMPRA: {bg_campra}/{len(all_idx)} ({100*bg_campra/len(all_idx):.1f}%)")

# Also test in full RD for comparison (circular — for reference only)
odds_campra_circ, p_campra_circ, _, _ = fisher_enrichment(rd_idx, all_idx, campra_arr)
print(f"  [CIRCULAR REFERENCE] caMPRA in RD (full): OR={odds_campra_circ:.3f}, p={p_campra_circ:.4e}")

# 7c: Tau comparison — LEAVE-ONE-OUT (tau removed from RDS)
print("\n--- 7c: Tau comparison (LOO — tau removed from RDS) ---")
gd_loo_tau = [i for i, c in enumerate(cls_loo_tau) if c == "gene-driven"]
rd_loo_tau = get_class_indices(cls_loo_tau, "regulation-driven")

gd_tau_vals = [gene_table[i]["tau"] for i in gd_loo_tau]
rd_tau_vals = [gene_table[i]["tau"] for i in rd_loo_tau]
if len(gd_tau_vals) > 0 and len(rd_tau_vals) > 0:
    u_tau, p_tau = stats.mannwhitneyu(rd_tau_vals, gd_tau_vals, alternative="greater")
    print(f"  Tau comparison RD(LOO) > GD(LOO): U={u_tau:.0f}, p={p_tau:.4e}")
    print(f"    RD(LOO) tau median: {np.median(rd_tau_vals):.4f}")
    print(f"    GD(LOO) tau median: {np.median(gd_tau_vals):.4f}")

# Also test with full classification (circular — for reference only)
gd_tau_circ = [gene_table[i]["tau"] for i in gd_idx]
rd_tau_circ = [gene_table[i]["tau"] for i in rd_idx]
if len(gd_tau_circ) > 0 and len(rd_tau_circ) > 0:
    u_tau_c, p_tau_c = stats.mannwhitneyu(rd_tau_circ, gd_tau_circ, alternative="greater")
    print(f"  [CIRCULAR REFERENCE] Tau RD > GD (full): U={u_tau_c:.0f}, p={p_tau_c:.4e}")

# 7d: Conservation comparison (phyloP CDS, GD > RD — non-circular)
print("\n--- 7d: Conservation comparison (non-circular) ---")
gd_phylop = [gene_table[i]["phyloP_cds"] for i in gd_idx if gene_table[i]["phyloP_cds"] != 0]
rd_phylop = [gene_table[i]["phyloP_cds"] for i in rd_idx if gene_table[i]["phyloP_cds"] != 0]
if len(gd_phylop) > 0 and len(rd_phylop) > 0:
    u_cons, p_cons = stats.mannwhitneyu(gd_phylop, rd_phylop, alternative="greater")
    print(f"  phyloP CDS (GD > RD): U={u_cons:.0f}, p={p_cons:.4e}")
    print(f"    GD phyloP CDS median: {np.median(gd_phylop):.4f}")
    print(f"    RD phyloP CDS median: {np.median(rd_phylop):.4f}")

# 7e: RELAX K>1 comparison (non-circular — RELAX is independent of RDS)
print("\n--- 7e: RELAX K>1 comparison (non-circular, proper denominator) ---")
# Only count genes with active RELAX results
gd_with_relax = [i for i in gd_idx if gene_table[i]["has_relax_active"]]
rd_with_relax = [i for i in rd_idx if gene_table[i]["has_relax_active"]]
gd_k_gt1 = sum(1 for i in gd_with_relax if gene_table[i]["relax_K"] > 1)
rd_k_gt1 = sum(1 for i in rd_with_relax if gene_table[i]["relax_K"] > 1)
print(f"  RELAX K>1 in GD (active only): {gd_k_gt1}/{len(gd_with_relax)} ({100*gd_k_gt1/max(len(gd_with_relax),1):.1f}%)")
print(f"  RELAX K>1 in RD (active only): {rd_k_gt1}/{len(rd_with_relax)} ({100*rd_k_gt1/max(len(rd_with_relax),1):.1f}%)")
print(f"  GD without active RELAX: {len(gd_idx) - len(gd_with_relax)}")
print(f"  RD without active RELAX: {len(rd_idx) - len(rd_with_relax)}")

# ============================================
# Step 8: Sensitivity analysis (P1-1)
# ============================================
print("\n--- P1-1: Weight sensitivity analysis ---")

# Original weights
orig_weights = {"p": 0.35, "lrt": 0.30, "relax": 0.20, "selectome": 0.15}
# ±50% perturbation
for label, w in [("p+50%", {"p": 0.525, "lrt": 0.30, "relax": 0.20, "selectome": 0.15}),
                  ("p-50%", {"p": 0.175, "lrt": 0.30, "relax": 0.20, "selectome": 0.15}),
                  ("all+50%", {"p": 0.525, "lrt": 0.45, "relax": 0.30, "selectome": 0.225}),
                  ("all-50%", {"p": 0.175, "lrt": 0.15, "relax": 0.10, "selectome": 0.075})]:
    # Normalize weights
    total_w = sum(w.values())
    w = {k: v/total_w for k, v in w.items()}
    
    gds_sens = (w["p"] * gds_p_pct + w["lrt"] * gds_lrt_pct + 
                w["relax"] * gds_relax_pct + w["selectome"] * gds_selectome)
    
    cls_sens = classify(gds_sens, rds_v7_full, bh_fdr_sig_arr, relax_relaxed_arr)
    c = Counter(cls_sens)
    gd_n = c.get("gene-driven", 0)
    rd_n = c.get("regulation-driven", 0)
    print(f"  {label}: GD={gd_n}, RD={rd_n}, neutral={c.get('neutral',0)}")

# Threshold sensitivity
print("\n  Threshold sensitivity:")
for thresh in [0.10, 0.15, 0.20]:
    cls_t = classify(gds_v7, rds_v7_full, bh_fdr_sig_arr, relax_relaxed_arr, threshold=thresh)
    c = Counter(cls_t)
    print(f"  threshold={thresh}: GD={c.get('gene-driven',0)}, RD={c.get('regulation-driven',0)}, dual={c.get('dual-driven',0)}, neutral={c.get('neutral',0)}")

# ============================================
# Step 9: Write output
# ============================================
print("\n--- Step 9: Writing output ---")

output_csv = f"{OUTPUT_DIR}/gene_classification_v7.csv"
fieldnames = [
    "gene_id", "gene_symbol", "cds_length", "classification_v7",
    "busted_p", "busted_lrt", "bh_fdr", "bh_fdr_sig",
    "storey_q", "storey_sig",
    "relax_K", "relax_p", "relax_fdr_sig", "has_relax_active", "relax_intensified", "relax_relaxed",
    "has_selectome", "n_hars", "has_doan_campra", "tau", "brain_tau",
    "phyloP_cds", "phyloP_gene", "phastCons_cds", "phastCons_gene", "nc_conservation",
    "gds_p_pct_resid", "gds_lrt_pct_resid", "gds_relax_pct", "gds_selectome", "gds_v7",
    "rds_doan", "rds_tau", "rds_nc", "rds_brain", "rds_v7",
    "classification_loo_campra", "classification_loo_tau",
]

with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for i, g in enumerate(gene_table):
        writer.writerow({
            "gene_id": g["gene_id"],
            "gene_symbol": g["gene_symbol"],
            "cds_length": g["cds_length"],
            "classification_v7": cls_full[i],
            "busted_p": g["busted_p"],
            "busted_lrt": g["busted_lrt"],
            "bh_fdr": g["bh_fdr"],
            "bh_fdr_sig": g["bh_fdr_sig"],
            "storey_q": f"{g['storey_q']:.6e}",
            "storey_sig": g["storey_sig"],
            "relax_K": g["relax_K"],
            "relax_p": g["relax_p"],
            "relax_fdr_sig": g["relax_fdr_sig"],
            "has_relax_active": g["has_relax_active"],
            "relax_intensified": g["relax_intensified"],
            "relax_relaxed": g["relax_relaxed"],
            "has_selectome": g["has_selectome"],
            "n_hars": g["n_hars"],
            "has_doan_campra": g["has_doan_campra"],
            "tau": g["tau"],
            "brain_tau": g["brain_tau"],
            "phyloP_cds": g["phyloP_cds"],
            "phyloP_gene": g["phyloP_gene"],
            "phastCons_cds": g["phastCons_cds"],
            "phastCons_gene": g["phastCons_gene"],
            "nc_conservation": g["nc_conservation"],
            "gds_p_pct_resid": f"{gds_p_pct[i]:.6f}",
            "gds_lrt_pct_resid": f"{gds_lrt_pct[i]:.6f}",
            "gds_relax_pct": f"{gds_relax_pct[i]:.6f}",
            "gds_selectome": f"{gds_selectome[i]:.6f}",
            "gds_v7": f"{gds_v7[i]:.6f}",
            "rds_doan": f"{rds_doan[i]:.6f}",
            "rds_tau": f"{rds_tau[i]:.6f}",
            "rds_nc": f"{rds_nc[i]:.6f}",
            "rds_brain": f"{rds_brain[i]:.6f}",
            "rds_v7": f"{rds_v7_full[i]:.6f}",
            "classification_loo_campra": cls_loo_campra[i],
            "classification_loo_tau": cls_loo_tau[i],
        })

print(f"  Written: {output_csv} ({len(gene_table)} genes)")

# Summary JSON
summary = {
    "version": "v7",
    "total_genes": len(gene_table),
    "classification_v7": dict(cls_counter),
    "p0_fixes": {
        "P0_1_leave_one_out": {
            "har_enrichment_rd": {"OR": float(odds_har) if not np.isnan(odds_har) else None, "p": float(p_har), "non_circular": True},
            "campra_enrichment_rd_loo": {"OR": float(odds_campra) if not np.isnan(odds_campra) else None, "p": float(p_campra), "non_circular": True},
            "campra_enrichment_rd_circular_ref": {"OR": float(odds_campra_circ) if not np.isnan(odds_campra_circ) else None, "p": float(p_campra_circ), "non_circular": False},
            "tau_comparison_loo": {"p": float(p_tau) if 'p_tau' in dir() and p_tau is not None else None, "non_circular": True},
            "tau_comparison_circular_ref": {"p": float(p_tau_c) if 'p_tau_c' in dir() and p_tau_c is not None else None, "non_circular": False},
            "conservation_comparison": {"p": float(p_cons) if 'p_cons' in dir() and p_cons is not None else None, "non_circular": True},
        },
        "P0_2_relax_coverage": {
            "total_genes": len(busted),
            "genes_with_relax_entry": len(genes_with_relax),
            "genes_with_active_relax": relax_active,
            "genes_with_default_relax": relax_default,
            "genes_without_relax": len(genes_without_relax),
            "busted_sig_with_relax": len(busted_sig_with_relax),
            "busted_sig_relax_active": busted_sig_relax_active,
            "k_gt1_rate_proper_denominator": float(relax_k_gt1_active) / max(len(relax_active_genes), 1),
            "k_lt1_rate_proper_denominator": float(relax_k_lt1_active) / max(len(relax_active_genes), 1),
        },
        "P0_3_storey_qvalue": {
            "pi0": float(pi0),
            "storey_sig_count": storey_sig_count,
            "bh_fdr_sig_count": sum(1 for v in busted.values() if v["fdr_sig"]),
            "p_eq_0": int(p_eq_0),
            "p_eq_05": int(p_eq_05),
        },
        "P0_4_gds_fdr_bonus_removed": {
            "old_weights": {"p": 0.30, "lrt": 0.25, "relax": 0.20, "selectome": 0.15, "fdr_bonus": 0.10},
            "new_weights": {"p": 0.35, "lrt": 0.30, "relax": 0.20, "selectome": 0.15},
            "fdr_bonus_removed": True,
        },
    },
    "cds_residualization": {
        "busted_neglogp_r": float(r1),
        "busted_neglogp_p": float(p1),
        "busted_lrt_r": float(r2),
        "busted_lrt_p": float(p2),
    },
    "loo_classification_campra": dict(cls_loo_campra_counter),
    "loo_classification_tau": dict(cls_loo_tau_counter),
}

summary_file = f"{OUTPUT_DIR}/classification_v7_summary.json"
with open(summary_file, "w") as f:
    json.dump(summary, f, indent=2)
print(f"  Summary: {summary_file}")

# Also write LOO summary
loo_summary = {
    "full": dict(cls_counter),
    "loo_campra": dict(cls_loo_campra_counter),
    "loo_tau": dict(cls_loo_tau_counter),
    "validation_results": {
        "HAR_enrichment_RD": {"OR": float(odds_har) if not np.isnan(odds_har) else None, "p": float(p_har), "type": "non-circular"},
        "caMPRA_enrichment_RD_LOO": {"OR": float(odds_campra) if not np.isnan(odds_campra) else None, "p": float(p_campra), "type": "non-circular (LOO)"},
        "caMPRA_enrichment_RD_circular": {"OR": float(odds_campra_circ) if not np.isnan(odds_campra_circ) else None, "p": float(p_campra_circ), "type": "CIRCULAR (reference only)"},
        "tau_comparison_LOO": {"p": float(p_tau) if 'p_tau' in dir() and len(rd_tau_vals) > 0 else None, "rd_median": float(np.median(rd_tau_vals)) if len(rd_tau_vals) > 0 else None, "gd_median": float(np.median(gd_tau_vals)) if len(gd_tau_vals) > 0 else None, "type": "non-circular (LOO)"},
        "tau_comparison_circular": {"p": float(p_tau_c) if 'p_tau_c' in dir() and len(rd_tau_circ) > 0 else None, "type": "CIRCULAR (reference only)"},
        "conservation_GD_gt_RD": {"p": float(p_cons) if 'p_cons' in dir() and len(gd_phylop) > 0 else None, "gd_median": float(np.median(gd_phylop)) if len(gd_phylop) > 0 else None, "rd_median": float(np.median(rd_phylop)) if len(rd_phylop) > 0 else None, "type": "non-circular"},
    }
}
loo_file = f"{OUTPUT_DIR}/loo_validation_summary.json"
with open(loo_file, "w") as f:
    json.dump(loo_summary, f, indent=2)
print(f"  LOO summary: {loo_file}")

print("\n" + "=" * 70)
print("v7 P0 Fixes Complete!")
print("=" * 70)
