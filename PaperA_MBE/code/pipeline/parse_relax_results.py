#!/usr/bin/env python3
"""Parse existing RELAX JSON results and identify remaining genes to run."""
import json
import os
import sys
from pathlib import Path

FASTA_DIR = "/home/linux/hyphy_work_v3/fasta"
FDR_LIST = "/home/linux/relax_work_v3/fdr_genes.txt"
OUTPUT_CSV = "/mnt/d/人类正选择基因项目/results/phase4_hyphy/relax_results/relax_results_v3.csv"
REMAINING_LIST = "/home/linux/relax_work_v3/remaining_genes_v2.txt"

# Read FDR gene list
with open(FDR_LIST) as f:
    fdr_genes = [line.strip() for line in f if line.strip()]
print(f"FDR significant genes: {len(fdr_genes)}")

# Scan all RELAX JSON files
valid_results = []
empty_jsons = []
missing_genes = []

for gene in fdr_genes:
    json_path = os.path.join(FASTA_DIR, f"{gene}.fasta.RELAX.json")
    if not os.path.exists(json_path):
        missing_genes.append(gene)
        continue
    
    size = os.path.getsize(json_path)
    if size < 100:
        empty_jsons.append(gene)
        continue
    
    try:
        with open(json_path) as f:
            data = json.load(f)
        
        test_results = data.get("test results", {})
        lrt = test_results.get("LRT", None)
        p_value = test_results.get("p-value", None)
        
        # Extract K - stored as "relaxation or intensification parameter" in test results
        k_value = test_results.get("relaxation or intensification parameter", None)
        
        # Extract omega rates if available
        fits = data.get("fits", {})
        
        valid_results.append({
            "gene_id": gene,
            "K": k_value,
            "LRT": lrt,
            "p_value": p_value,
            "status": "complete"
        })
    except json.JSONDecodeError:
        empty_jsons.append(gene)
    except Exception as e:
        print(f"Error parsing {gene}: {e}", file=sys.stderr)
        empty_jsons.append(gene)

# Compute BH FDR
import numpy as np
from scipy import stats

p_values = np.array([r["p_value"] for r in valid_results if r["p_value"] is not None])
n = len(p_values)
if n > 0:
    sorted_indices = np.argsort(p_values)
    bh_fdr = np.zeros(n)
    for i, idx in enumerate(sorted_indices):
        rank = i + 1
        bh_fdr[idx] = p_values[idx] * n / rank
    # Enforce monotonicity
    for i in range(n - 2, -1, -1):
        bh_fdr[sorted_indices[i]] = min(bh_fdr[sorted_indices[i]], bh_fdr[sorted_indices[i + 1]])
    bh_fdr = np.minimum(bh_fdr, 1.0)
    
    fdr_idx = 0
    for r in valid_results:
        if r["p_value"] is not None:
            r["bh_fdr"] = bh_fdr[fdr_idx]
            r["fdr_significant"] = bool(r["bh_fdr"] < 0.05)
            r["k_category"] = "intensified" if r["K"] is not None and r["K"] > 1 else ("relaxed" if r["K"] is not None and r["K"] < 1 else "unknown")
            fdr_idx += 1
        else:
            r["bh_fdr"] = None
            r["fdr_significant"] = False
            r["k_category"] = "unknown"

# Write CSV
os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
with open(OUTPUT_CSV, "w") as f:
    f.write("gene_id,K,LRT,p_value,bh_fdr,fdr_significant,k_category,status\n")
    for r in valid_results:
        k_str = f"{r['K']:.6f}" if r['K'] is not None else ""
        lrt_str = f"{r['LRT']:.6f}" if r['LRT'] is not None else ""
        p_str = f"{r['p_value']:.6e}" if r['p_value'] is not None else ""
        fdr_str = f"{r['bh_fdr']:.6e}" if r.get('bh_fdr') is not None else ""
        f.write(f"{r['gene_id']},{k_str},{lrt_str},{p_str},{fdr_str},{r.get('fdr_significant',False)},{r.get('k_category','unknown')},{r['status']}\n")

# Write remaining genes (need to be run)
remaining = missing_genes + empty_jsons
with open(REMAINING_LIST, "w") as f:
    for g in remaining:
        f.write(g + "\n")

# Summary
k_values = [r["K"] for r in valid_results if r["K"] is not None]
k_greater_1 = sum(1 for k in k_values if k > 1)
k_less_1 = sum(1 for k in k_values if k < 1)
fdr_sig = sum(1 for r in valid_results if r.get("fdr_significant", False))

print(f"\n=== RELAX Results Summary ===")
print(f"Valid JSON results: {len(valid_results)}")
print(f"Empty/failed JSONs: {len(empty_jsons)}")
print(f"Missing (not run): {len(missing_genes)}")
print(f"Remaining to run: {len(remaining)}")
print(f"\nK > 1 (intensified selection): {k_greater_1}")
print(f"K < 1 (relaxed constraint): {k_less_1}")
print(f"K = 1 (no change): {len(k_values) - k_greater_1 - k_less_1}")
print(f"FDR < 0.05: {fdr_sig}")
if k_values:
    print(f"\nK range: {min(k_values):.4f} - {max(k_values):.4f}")
    print(f"K median: {np.median(k_values):.4f}")
print(f"\nCSV written: {OUTPUT_CSV}")
print(f"Remaining gene list: {REMAINING_LIST} ({len(remaining)} genes)")
