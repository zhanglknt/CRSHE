#!/usr/bin/env python3
"""
Phase 8: Extract BUSTED v1 Omega values from JSON files

Parses 4,832 BUSTED v1 JSON files to extract omega values from the
Unconstrained model → Test (human foreground branch) rate distribution.

Extracts:
- omega_1, omega_2, omega_3 (three rate classes)
- proportion_1, proportion_2, proportion_3
- omega_max (max of the three)
- omega_weighted (sum of omega_i * proportion_i)
- omega_positive_class (omega > 1, if any)
- proportion_positive (proportion of positive selection class)
- has_positive_selection (any omega > 1)

Run in WSL: python3 phase8_extract_omega_v1.py
"""
import json
import csv
import os
import sys
import time
from pathlib import Path

BASE = "/mnt/d/hs_gene_project"
BUSTED_DIR = Path(f"{BASE}/results/phase4_hyphy/busted_results_v1_backup")
OUTPUT_DIR = Path(f"{BASE}/results/phase8_tissue_analysis")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_CSV = OUTPUT_DIR / "busted_v1_omega.csv"

print("=" * 70)
print("Phase 8: BUSTED v1 Omega Extraction")
print("=" * 70)

# ============================================
# Parse each JSON and extract omega
# ============================================

json_files = sorted(BUSTED_DIR.glob("*.json"))
print(f"\nFound {len(json_files)} BUSTED v1 JSON files")

results = []
errors = []
t0 = time.time()

for i, jf in enumerate(json_files):
    gene_id = jf.stem  # ENSG00000xxxxx
    
    try:
        with open(jf) as f:
            d = json.load(f)
        
        # Navigate to Unconstrained model → Rate Distributions → Test
        fits = d.get("fits", {})
        unconstrained = fits.get("Unconstrained model", {})
        rate_dist = unconstrained.get("Rate Distributions", {})
        test_rates = rate_dist.get("Test", {})
        
        if not test_rates:
            # Try alternative key
            for key in rate_dist:
                if "test" in key.lower():
                    test_rates = rate_dist[key]
                    break
        
        if not test_rates:
            errors.append((gene_id, "No Test rate distribution found"))
            continue
        
        # Extract omega and proportion for each class
        omegas = []
        proportions = []
        
        for class_key in sorted(test_rates.keys()):
            cls_data = test_rates[class_key]
            if isinstance(cls_data, dict):
                omega = float(cls_data.get("omega", 0))
                prop = float(cls_data.get("proportion", 0))
                omegas.append(omega)
                proportions.append(prop)
        
        if len(omegas) < 2:
            errors.append((gene_id, f"Only {len(omegas)} rate classes"))
            continue
        
        # Pad to 3 classes if needed
        while len(omegas) < 3:
            omegas.append(0.0)
            proportions.append(0.0)
        
        # Compute derived metrics
        omega_max = max(omegas)
        omega_weighted = sum(o * p for o, p in zip(omegas, proportions))
        
        # Find positive selection class (omega > 1)
        positive_omegas = [(o, p) for o, p in zip(omegas, proportions) if o > 1.0]
        if positive_omegas:
            omega_positive = max(positive_omegas, key=lambda x: x[0])[0]
            proportion_positive = sum(p for _, p in positive_omegas)
            has_pos_sel = True
        else:
            omega_positive = 0.0
            proportion_positive = 0.0
            has_pos_sel = False
        
        # Also extract background omega for comparison
        bg_rates = rate_dist.get("Background", {})
        bg_omegas = []
        bg_proportions = []
        for class_key in sorted(bg_rates.keys()) if isinstance(bg_rates, dict) else []:
            cls_data = bg_rates[class_key]
            if isinstance(cls_data, dict):
                bg_omegas.append(float(cls_data.get("omega", 0)))
                bg_proportions.append(float(cls_data.get("proportion", 0)))
        
        bg_omega_max = max(bg_omegas) if bg_omegas else 0
        bg_omega_weighted = sum(o * p for o, p in zip(bg_omegas, bg_proportions)) if bg_omegas else 0
        
        results.append({
            "gene_id": gene_id,
            "omega_1": omegas[0],
            "omega_2": omegas[1],
            "omega_3": omegas[2] if len(omegas) > 2 else 0,
            "proportion_1": proportions[0],
            "proportion_2": proportions[1],
            "proportion_3": proportions[2] if len(proportions) > 2 else 0,
            "omega_max": omega_max,
            "omega_weighted": omega_weighted,
            "omega_positive": omega_positive,
            "proportion_positive": proportion_positive,
            "has_positive_selection": has_pos_sel,
            "bg_omega_max": bg_omega_max,
            "bg_omega_weighted": bg_omega_weighted,
        })
        
    except json.JSONDecodeError as e:
        errors.append((gene_id, f"JSON decode error: {e}"))
    except Exception as e:
        errors.append((gene_id, f"Error: {e}"))
    
    # Progress
    if (i + 1) % 1000 == 0:
        elapsed = time.time() - t0
        rate = (i + 1) / elapsed
        eta = (len(json_files) - i - 1) / rate
        print(f"  Processed {i+1}/{len(json_files)} ({100*(i+1)/len(json_files):.1f}%) "
              f"| {rate:.0f} genes/s | ETA: {eta:.0f}s")

elapsed = time.time() - t0
print(f"\nDone: {len(results)} extracted, {len(errors)} errors in {elapsed:.1f}s")

# ============================================
# Write output CSV
# ============================================
print(f"\n--- Writing output ---")

fieldnames = [
    "gene_id",
    "omega_1", "omega_2", "omega_3",
    "proportion_1", "proportion_2", "proportion_3",
    "omega_max", "omega_weighted",
    "omega_positive", "proportion_positive",
    "has_positive_selection",
    "bg_omega_max", "bg_omega_weighted",
]

with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for r in results:
        writer.writerow(r)

print(f"  Output: {OUTPUT_CSV} ({len(results)} genes)")

# ============================================
# Summary statistics
# ============================================
print(f"\n--- Summary Statistics ---")

import numpy as np

omega_max_arr = np.array([r["omega_max"] for r in results])
omega_weighted_arr = np.array([r["omega_weighted"] for r in results])
prop_positive_arr = np.array([r["proportion_positive"] for r in results])
has_pos_sel_arr = np.array([r["has_positive_selection"] for r in results])

print(f"  Genes with omega_max > 1: {np.sum(omega_max_arr > 1)} ({100*np.sum(omega_max_arr > 1)/len(results):.1f}%)")
print(f"  Genes with omega_max > 5: {np.sum(omega_max_arr > 5)} ({100*np.sum(omega_max_arr > 5)/len(results):.1f}%)")
print(f"  Genes with omega_max > 10: {np.sum(omega_max_arr > 10)} ({100*np.sum(omega_max_arr > 10)/len(results):.1f}%)")
print(f"  omega_max median: {np.median(omega_max_arr):.4f}")
print(f"  omega_max mean: {np.mean(omega_max_arr):.4f}")
print(f"  omega_weighted median: {np.median(omega_weighted_arr):.4f}")
print(f"  proportion_positive median (among genes with pos sel): "
      f"{np.median(prop_positive_arr[has_pos_sel_arr]):.4f}")

# Top 10 by omega_max
print(f"\n  Top 10 by omega_max:")
sorted_results = sorted(results, key=lambda x: x["omega_max"], reverse=True)
for r in sorted_results[:10]:
    print(f"    {r['gene_id']}: omega_max={r['omega_max']:.2f}, "
          f"omega_weighted={r['omega_weighted']:.4f}, "
          f"prop_pos={r['proportion_positive']:.4f}")

# Error summary
if errors:
    print(f"\n  Errors ({len(errors)}):")
    for gid, err in errors[:5]:
        print(f"    {gid}: {err}")
    if len(errors) > 5:
        print(f"    ... and {len(errors)-5} more")

print("\n" + "=" * 70)
print("BUSTED v1 Omega Extraction Complete!")
print("=" * 70)
