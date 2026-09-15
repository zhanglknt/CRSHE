#!/usr/bin/env python3
"""
Phase 8 Pre-step: Add classification_loo_brain variant

LOO-brain: Remove brain_tau (0.20 weight) from RDS
Redistribute: caMPRA=0.375, tau=0.3125, nc=0.3125
(Original: caMPRA=0.30, tau=0.25, nc=0.25 → normalize to sum=1.0)

Run in WSL: python3 phase8_loo_brain.py
"""
import csv
import json
import os
import numpy as np
from scipy import stats
from collections import Counter

BASE = os.environ.get("HSD_BASE", "/mnt/d/hs_gene_project")
CLASSIFICATION_CSV = f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv"
SUMMARY_JSON = f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7/classification_v7_summary.json"

print("=" * 70)
print("Phase 8 Pre-step: LOO-Brain Classification Variant")
print("=" * 70)

# ============================================
# Step 1: Load existing v7 classification
# ============================================
print("\n--- Step 1: Loading v7 classification ---")

gene_data = []
with open(CLASSIFICATION_CSV) as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        gene_data.append(row)

print(f"  Loaded {len(gene_data)} genes")
print(f"  Existing columns: {len(fieldnames)}")

# Extract arrays needed for classification
gds_v7 = np.array([float(g["gds_v7"]) for g in gene_data])
rds_doan = np.array([float(g["rds_doan"]) for g in gene_data])
rds_tau = np.array([float(g["rds_tau"]) for g in gene_data])
rds_nc = np.array([float(g["rds_nc"]) for g in gene_data])
rds_brain = np.array([float(g["rds_brain"]) for g in gene_data])
rds_v7_full = np.array([float(g["rds_v7"]) for g in gene_data])
bh_fdr_sig = np.array([g["bh_fdr_sig"] == "True" for g in gene_data])
relax_relaxed = np.array([g["relax_relaxed"] == "True" for g in gene_data])

# Verify RDS components sum correctly
rds_check = (0.30 * rds_doan + 0.25 * rds_tau + 0.25 * rds_nc + 0.20 * rds_brain)
diff = np.max(np.abs(rds_check - rds_v7_full))
print(f"  RDS verification: max diff = {diff:.10f} (should be ~0)")

# ============================================
# Step 2: Compute LOO-brain RDS
# ============================================
print("\n--- Step 2: Computing LOO-brain RDS ---")

# Remove brain_tau (0.20 weight), redistribute remaining:
# caMPRA: 0.30/0.80 = 0.375
# tau:    0.25/0.80 = 0.3125
# nc:     0.25/0.80 = 0.3125
rds_v7_no_brain = (0.375 * rds_doan + 0.3125 * rds_tau + 0.3125 * rds_nc)

print(f"  RDS v7 (no brain): {rds_v7_no_brain.min():.4f} - {rds_v7_no_brain.max():.4f}")
print(f"  RDS v7 (full):     {rds_v7_full.min():.4f} - {rds_v7_full.max():.4f}")
print(f"  Correlation (full vs no-brain): r={np.corrcoef(rds_v7_full, rds_v7_no_brain)[0,1]:.4f}")

# ============================================
# Step 3: Classify with LOO-brain
# ============================================
print("\n--- Step 3: Classification with LOO-brain ---")

def classify(gds, rds, storey_sig, relax_relaxed_flags, threshold=0.15):
    """Classify genes based on GDS and RDS scores (same as v7)."""
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

cls_loo_brain = classify(gds_v7, rds_v7_no_brain, bh_fdr_sig, relax_relaxed)

cls_loo_brain_counter = Counter(cls_loo_brain)
print("\n  Classification v7 (LOO brain):")
for k, v in cls_loo_brain_counter.most_common():
    print(f"    {k}: {v}")

# Compare with full classification
cls_full = [g["classification_v7"] for g in gene_data]
cls_full_counter = Counter(cls_full)
print("\n  Classification v7 (full, for reference):")
for k, v in cls_full_counter.most_common():
    print(f"    {k}: {v}")

# Count changes
changes = sum(1 for i in range(len(gene_data)) if cls_loo_brain[i] != cls_full[i])
print(f"\n  Genes that changed classification: {changes} ({100*changes/len(gene_data):.1f}%)")

# ============================================
# Step 4: Non-circular brain_tau validation
# ============================================
print("\n--- Step 4: Non-circular brain_tau validation ---")

gd_loo_brain = [i for i, c in enumerate(cls_loo_brain) if c == "gene-driven"]
rd_loo_brain = [i for i, c in enumerate(cls_loo_brain) if c == "regulation-driven"]

gd_brain_tau = [float(gene_data[i]["brain_tau"]) for i in gd_loo_brain]
rd_brain_tau = [float(gene_data[i]["brain_tau"]) for i in rd_loo_brain]

if len(gd_brain_tau) > 0 and len(rd_brain_tau) > 0:
    u_brain, p_brain = stats.mannwhitneyu(rd_brain_tau, gd_brain_tau, alternative="greater")
    print(f"  Brain tau comparison RD(LOO) > GD(LOO): U={u_brain:.0f}, p={p_brain:.4e}")
    print(f"    RD(LOO) brain_tau median: {np.median(rd_brain_tau):.4f} (n={len(rd_brain_tau)})")
    print(f"    GD(LOO) brain_tau median: {np.median(gd_brain_tau):.4f} (n={len(gd_brain_tau)})")
else:
    u_brain, p_brain = None, None
    print(f"  Insufficient data: GD={len(gd_brain_tau)}, RD={len(rd_brain_tau)}")

# Circular reference (for comparison)
gd_full = [i for i, c in enumerate(cls_full) if c == "gene-driven"]
rd_full = [i for i, c in enumerate(cls_full) if c == "regulation-driven"]
gd_brain_circ = [float(gene_data[i]["brain_tau"]) for i in gd_full]
rd_brain_circ = [float(gene_data[i]["brain_tau"]) for i in rd_full]
if len(gd_brain_circ) > 0 and len(rd_brain_circ) > 0:
    u_brain_c, p_brain_c = stats.mannwhitneyu(rd_brain_circ, gd_brain_circ, alternative="greater")
    print(f"  [CIRCULAR REF] Brain tau RD > GD (full): U={u_brain_c:.0f}, p={p_brain_c:.4e}")

# ============================================
# Step 5: Write updated CSV with LOO-brain column
# ============================================
print("\n--- Step 5: Writing updated CSV ---")

# Add new columns
new_fields = list(fieldnames)
if "rds_v7_no_brain" not in new_fields:
    new_fields.append("rds_v7_no_brain")
if "classification_loo_brain" not in new_fields:
    new_fields.append("classification_loo_brain")

output_csv = CLASSIFICATION_CSV  # Overwrite in place
with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=new_fields)
    writer.writeheader()
    for i, g in enumerate(gene_data):
        row = dict(g)
        row["rds_v7_no_brain"] = f"{rds_v7_no_brain[i]:.6f}"
        row["classification_loo_brain"] = cls_loo_brain[i]
        writer.writerow(row)

print(f"  Updated: {output_csv} ({len(gene_data)} genes, {len(new_fields)} columns)")
print(f"  New columns: rds_v7_no_brain, classification_loo_brain")

# ============================================
# Step 6: Update summary JSON
# ============================================
print("\n--- Step 6: Updating summary JSON ---")

if os.path.exists(SUMMARY_JSON):
    with open(SUMMARY_JSON) as f:
        summary = json.load(f)
else:
    summary = {}

summary["loo_classification_brain"] = dict(cls_loo_brain_counter)
summary["p0_fixes"]["P0_1_leave_one_out"]["brain_tau_comparison_loo"] = {
    "p": float(p_brain) if p_brain is not None else None,
    "rd_median": float(np.median(rd_brain_tau)) if len(rd_brain_tau) > 0 else None,
    "gd_median": float(np.median(gd_brain_tau)) if len(gd_brain_tau) > 0 else None,
    "rd_n": len(rd_brain_tau),
    "gd_n": len(gd_brain_tau),
    "non_circular": True,
}
summary["p0_fixes"]["P0_1_leave_one_out"]["brain_tau_comparison_circular_ref"] = {
    "p": float(p_brain_c) if p_brain_c is not None else None,
    "non_circular": False,
}

with open(SUMMARY_JSON, "w") as f:
    json.dump(summary, f, indent=2)
print(f"  Updated: {SUMMARY_JSON}")

# Also update LOO summary
loo_file = f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7/loo_validation_summary.json"
if os.path.exists(loo_file):
    with open(loo_file) as f:
        loo_summary = json.load(f)
else:
    loo_summary = {}

loo_summary["loo_brain"] = dict(cls_loo_brain_counter)
loo_summary["validation_results"]["brain_tau_comparison_LOO"] = {
    "p": float(p_brain) if p_brain is not None else None,
    "rd_median": float(np.median(rd_brain_tau)) if len(rd_brain_tau) > 0 else None,
    "gd_median": float(np.median(gd_brain_tau)) if len(gd_brain_tau) > 0 else None,
    "type": "non-circular (LOO)",
}
loo_summary["validation_results"]["brain_tau_comparison_circular"] = {
    "p": float(p_brain_c) if p_brain_c is not None else None,
    "type": "CIRCULAR (reference only)",
}

with open(loo_file, "w") as f:
    json.dump(loo_summary, f, indent=2)
print(f"  Updated: {loo_file}")

print("\n" + "=" * 70)
print("LOO-Brain Classification Complete!")
print("=" * 70)
