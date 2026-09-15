#!/usr/bin/env python3
"""
Phase 8 extension: HPA whole-body single-cell-type enrichment (154 cell types, v25.1)
=====================================================================================
Question: is the regulation-driven (RD) high-expression signal specific to neuronal
cell types when compared against ALL other organs' cell types (hepatocytes,
cardiomyocytes, immune cells, ...)? The brain single-nuclei analysis (34 types)
was brain-only; this closes the whole-body gap.

Data: HPA v25.1 rna_single_cell_type.tsv (long format: Gene, Gene name, Cell type, nCPM)
      + rna_single_cell_type_cell_types.tsv (cell type -> group -> class metadata)

Three classification variants are tested to control circularity:
  - v7 (full): reference, partially circular (RDS includes brain-tau and GTEx-tau)
  - LOO-brain: non-circular for neuronal/brain cell types
  - LOO-tau: non-circular for expression-specificity-driven signal

GD excludes "gene-driven (relaxed)"; dual-driven excluded from both (Table 3 convention).

Outputs (results/phase8_tissue_analysis/):
  hpa_wholebody_celltype_enrichment.csv   (154 cell types x 3 variants)
  hpa_wholebody_cellclass_enrichment.csv  (15 cell classes x 3 variants)
  hpa_wholebody_summary.json
  figures/fig9_hpa_wholebody_cellclass.png/.pdf
"""
import csv
import json
import os
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(os.environ.get("HSD_BASE") or (
    "/mnt/d/人类正选择基因项目" if os.path.exists("/mnt/d") else "D:/人类正选择基因项目"))

MASTER_CSV = BASE / "results/phase8_tissue_analysis/data/phase8_master_table.csv"
HPA_EXPR = BASE / "data/hpa_single_cell/rna_single_cell_type.tsv"
HPA_META = BASE / "data/hpa_single_cell/rna_single_cell_type_cell_types.tsv"
OUT_DIR = BASE / "results/phase8_tissue_analysis"
FIG_DIR = OUT_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

print("=" * 70)
print("Phase 8 extension: HPA whole-body single-cell-type enrichment")
print("=" * 70)


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    fdr = np.empty(n)
    fdr[order[-1]] = pvals[order[-1]]
    for i in range(n - 2, -1, -1):
        fdr[order[i]] = min(pvals[order[i]] * n / (i + 1), fdr[order[i + 1]])
    return fdr


def fisher_greater(group_set, bg_set, feature_set):
    a = len(group_set & feature_set)
    b = len(group_set - feature_set)
    c = len(feature_set - group_set)
    d = len(bg_set - group_set - feature_set)
    if a + b == 0 or c + d == 0:
        return np.nan, 1.0, 0
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return odds, p, a


# ---------------- load classifications ----------------
print("\n--- Loading master table ---")
variants = {
    "v7": "classification_v7",
    "loo_brain": "classification_loo_brain",
    "loo_tau": "classification_loo_tau",
    "loo_campra": "classification_loo_campra",
    "loo_double": "classification_loo_double",  # computed below (not in master)
}
with open(MASTER_CSV) as f:
    rows = list(csv.DictReader(f))

# --- double-LOO variant (remove GTEx-tau AND brain-tau from RDS) ---
# RDS_double = (0.30*rds_doan + 0.25*rds_nc) / 0.55; GDS unchanged (gds_v7).
# Same classify() rule as phase7g_v7_classification.py (threshold 0.15, BH-FDR sig gate).
def _classify(g, r, sig, relaxed, thr=0.15):
    if g > r + thr and sig:
        cls = "gene-driven"
    elif r > g + thr and not sig:
        cls = "regulation-driven"
    elif g > r + thr and r > 0.4 and sig:
        cls = "dual-driven"
    elif r > g + thr and sig and g > 0.3:
        cls = "dual-driven"
    else:
        cls = "neutral"
    if relaxed and cls == "gene-driven":
        cls = "gene-driven (relaxed)"
    return cls

if "classification_loo_double" not in rows[0]:
    for r in rows:
        g = float(r["gds_v7"])
        rd = (0.30 * float(r["rds_doan"]) + 0.25 * float(r["rds_nc"])) / 0.55
        sig = str(r["bh_fdr_sig"]) in ("True", "1", "true")
        relaxed = str(r["relax_relaxed"]) in ("True", "1", "true")
        r["classification_loo_double"] = _classify(g, rd, sig, relaxed)
    gd_v7 = {r["gene_id"] for r in rows if r["classification_v7"] == "gene-driven"}
    gd_db = {r["gene_id"] for r in rows if r["classification_loo_double"] == "gene-driven"}
    print(f"  [double-LOO] GD n={len(gd_db)} vs v7 GD n={len(gd_v7)}; "
          f"lost={len(gd_v7 - gd_db)}, gained={len(gd_db - gd_v7)} "
          f"(GD set identical to v7: {gd_v7 == gd_db})")

genesets = {}  # variant -> {"gd": set, "rd": set}
all_genes = set()
for vname, col in variants.items():
    gd = {r["gene_id"] for r in rows if r.get(col) == "gene-driven"}
    rd = {r["gene_id"] for r in rows if r.get(col) == "regulation-driven"}
    genesets[vname] = {"gd": gd, "rd": rd}
    print(f"  {vname:10s}: GD={len(gd):5d}  RD={len(rd):5d}")
all_genes = {r["gene_id"] for r in rows}
print(f"  Universe: {len(all_genes)}")

# ---------------- load cell type metadata ----------------
print("\n--- Loading cell type metadata ---")
ct_group, ct_class = {}, {}
with open(HPA_META) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        ct_group[r["Cell type"]] = r["Cell type group"]
        ct_class[r["Cell type"]] = r["Cell type class"]
print(f"  {len(ct_class)} cell types, {len(set(ct_class.values()))} classes")

# ---------------- parse expression (long -> wide) ----------------
print("\n--- Parsing 142MB expression TSV ---")
gene_ct = defaultdict(dict)  # gene -> {cell_type: nCPM}
with open(HPA_EXPR, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        gid = row["Gene"].split(".")[0]
        if gid not in all_genes:
            continue
        try:
            gene_ct[gid][row["Cell type"]] = float(row["nCPM"] or 0)
        except ValueError:
            pass
matched = len(gene_ct)
all_cts = sorted({ct for v in gene_ct.values() for ct in v})
print(f"  Matched genes: {matched} / {len(all_genes)}; cell types: {len(all_cts)}")

# ---------------- per-cell-type enrichment ----------------
print("\n--- Per-cell-type enrichment (154 types x 3 variants) ---")


def high_expr_set(expr_vals):
    vals = np.asarray(list(expr_vals.values()))
    nz = vals[vals > 0]
    if len(nz) < 50:
        return None
    thr = np.percentile(nz, 75)
    he = {g for g, v in expr_vals.items() if v >= thr and v > 0}
    return he if len(he) >= 10 else None


results = []
for ct in all_cts:
    expr = {g: v.get(ct, 0.0) for g, v in gene_ct.items()}
    if len(expr) < 100:
        continue
    he = high_expr_set(expr)
    if he is None:
        continue
    bg = set(expr)
    rec = {"cell_type": ct,
           "cell_type_group": ct_group.get(ct, "NA"),
           "cell_type_class": ct_class.get(ct, "NA"),
           "n_genes": len(expr), "n_high_expr": len(he)}
    for vname in variants:
        gd_odds, gd_p, gd_n = fisher_greater(genesets[vname]["gd"] & bg, bg, he)
        rd_odds, rd_p, rd_n = fisher_greater(genesets[vname]["rd"] & bg, bg, he)
        rec[f"{vname}_gd_OR"] = gd_odds
        rec[f"{vname}_gd_p"] = gd_p
        rec[f"{vname}_gd_n"] = gd_n
        rec[f"{vname}_rd_OR"] = rd_odds
        rec[f"{vname}_rd_p"] = rd_p
        rec[f"{vname}_rd_n"] = rd_n
    results.append(rec)

for vname in variants:
    for tgt in ("gd", "rd"):
        fdrs = bh_fdr([r[f"{vname}_{tgt}_p"] for r in results])
        for i, r in enumerate(results):
            r[f"{vname}_{tgt}_fdr"] = fdrs[i]

out_csv = OUT_DIR / "hpa_wholebody_celltype_enrichment.csv"
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader()
    w.writerows(results)
print(f"  Output: {out_csv} ({len(results)} cell types)")

# ---------------- per-class enrichment ----------------
print("\n--- Per-class enrichment (15 classes x 3 variants) ---")
class_results = []
for cc in sorted(set(ct_class.values())):
    cts = [ct for ct in all_cts if ct_class.get(ct) == cc]
    if not cts:
        continue
    expr = {g: float(np.mean([v.get(ct, 0.0) for ct in cts])) for g, v in gene_ct.items()}
    he = high_expr_set(expr)
    if he is None:
        continue
    bg = set(expr)
    rec = {"cell_type_class": cc, "n_types": len(cts),
           "n_genes": len(expr), "n_high_expr": len(he)}
    for vname in variants:
        gd_odds, gd_p, gd_n = fisher_greater(genesets[vname]["gd"] & bg, bg, he)
        rd_odds, rd_p, rd_n = fisher_greater(genesets[vname]["rd"] & bg, bg, he)
        rec[f"{vname}_gd_OR"] = gd_odds
        rec[f"{vname}_gd_p"] = gd_p
        rec[f"{vname}_rd_OR"] = rd_odds
        rec[f"{vname}_rd_p"] = rd_p
        rec[f"{vname}_rd_n"] = rd_n
    class_results.append(rec)

for vname in variants:
    for tgt in ("gd", "rd"):
        fdrs = bh_fdr([r[f"{vname}_{tgt}_p"] for r in class_results])
        for i, r in enumerate(class_results):
            r[f"{vname}_{tgt}_fdr"] = fdrs[i]

out_csv2 = OUT_DIR / "hpa_wholebody_cellclass_enrichment.csv"
with open(out_csv2, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(class_results[0].keys()))
    w.writeheader()
    w.writerows(class_results)
print(f"  Output: {out_csv2} ({len(class_results)} classes)")

# ---------------- key comparison: neuronal vs everything else ----------------
print("\n--- Neuronal vs non-neuronal comparison ---")
summary = {"n_celltypes": len(results), "n_classes": len(class_results),
           "n_matched_genes": matched}

for vname in variants:
    neu = [r[f"{vname}_rd_OR"] for r in results if r["cell_type_class"] == "neuronal cells"]
    non = [r[f"{vname}_rd_OR"] for r in results if r["cell_type_class"] != "neuronal cells"]
    gli = [r[f"{vname}_rd_OR"] for r in results if r["cell_type_class"] == "glial cells"]
    u, p = stats.mannwhitneyu(neu, non, alternative="greater")
    ranked = sorted(results, key=lambda r: -(r[f"{vname}_rd_OR"] or 0))
    top10 = [(r["cell_type"], r["cell_type_class"], round(r[f"{vname}_rd_OR"], 2)) for r in ranked[:10]]
    n_rd_sig = sum(1 for r in results if r[f"{vname}_rd_fdr"] < 0.05)
    n_gd_sig = sum(1 for r in results if r[f"{vname}_gd_fdr"] < 0.05)
    neu_sig = sum(1 for r in results
                  if r["cell_type_class"] == "neuronal cells" and r[f"{vname}_rd_fdr"] < 0.05)
    summary[vname] = {
        "neuronal_rd_OR_median": float(np.median(neu)),
        "nonneuronal_rd_OR_median": float(np.median(non)),
        "glial_rd_OR_median": float(np.median(gli)),
        "mwu_neuronal_gt_non_p": float(p),
        "n_types_rd_fdr_sig": n_rd_sig,
        "n_types_gd_fdr_sig": n_gd_sig,
        "n_neuronal_types": len(neu),
        "n_neuronal_rd_sig": neu_sig,
        "top10_by_rd_OR": top10,
    }
    print(f"  [{vname}] neuronal median RD OR={np.median(neu):.2f} vs "
          f"non-neuronal={np.median(non):.2f} (MWU p={p:.2e}); "
          f"RD FDR<0.05 in {n_rd_sig}/{len(results)} types "
          f"(neuronal {neu_sig}/{len(neu)})")
    print(f"    Top-5 RD OR: {[(t[0], t[2]) for t in top10[:5]]}")

with open(OUT_DIR / "hpa_wholebody_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"  Output: {OUT_DIR / 'hpa_wholebody_summary.json'}")

# ---------------- figure ----------------
print("\n--- Figure: fig9_hpa_wholebody_cellclass ---")
vname = "loo_brain"  # primary non-circular variant for the neuronal question
cr = sorted(class_results, key=lambda r: r[f"{vname}_rd_OR"])
labels = [f"{r['cell_type_class']} ({r['n_types']})" for r in cr]
rd_or = [r[f"{vname}_rd_OR"] for r in cr]
gd_or = [r[f"{vname}_gd_OR"] for r in cr]
rd_sig = [r[f"{vname}_rd_fdr"] < 0.05 for r in cr]
colors = ["#d6604d" if r["cell_type_class"] == "neuronal cells"
          else "#f4a582" if s else "#fddbc7" for r, s in zip(cr, rd_sig)]

fig, ax = plt.subplots(figsize=(9, 6.5))
y = np.arange(len(cr))
ax.barh(y, rd_or, color=colors, edgecolor="white", height=0.7)
ax.scatter(gd_or, y, color="#2166ac", s=28, zorder=3)
ax.axvline(1.0, color="grey", ls="--", lw=1)
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel("Odds ratio of high-expression enrichment (top-quartile nCPM)")
ax.set_title("HPA whole-body single-cell types (154 types, 15 classes): RD vs GD enrichment\n"
             "(LOO-brain classification)")
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
ax.legend(handles=[
    Patch(fc="#d6604d", label="RD OR - neuronal (FDR<0.05)"),
    Patch(fc="#f4a582", label="RD OR - FDR<0.05"),
    Patch(fc="#fddbc7", label="RD OR - n.s."),
    Line2D([0], [0], marker="o", color="none", mfc="#2166ac", ms=7, label="GD OR (LOO-brain)"),
], fontsize=9, loc="lower right")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig9_hpa_wholebody_cellclass.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG_DIR / "fig9_hpa_wholebody_cellclass.pdf", bbox_inches="tight")
plt.close(fig)
print(f"  saved fig9_hpa_wholebody_cellclass.png/.pdf")

print("\n" + "=" * 70)
print("Whole-body single-cell analysis complete!")
print("=" * 70)
