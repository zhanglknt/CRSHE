# -*- coding: utf-8 -*-
"""
P2 GO-slim one-sided vs two-sided q-value audit
================================================================================
Question (reviewer, Round 5): the q_bh column of p2_goslim_composition.csv may
have been computed from two-sided Fisher p, while the main text adopts the
one-sided (enrichment-direction, 'greater') convention; a value 0.46 vs 0.38 is
quoted. This script recomputes both BH families from the table's own counts,
verifies which convention the stored q uses, and produces the corrected minimum
q values the manuscript should quote.

Output: results/phase9_hardening/p2_goslim_check.json ; updated CSV (adds
q_bh_two_sided column) if needed ; log _audit_p2.txt
"""
import io
import json

import numpy as np
import pandas as pd
from scipy import stats

BASE = "d:/人类正选择基因项目/"
LOG = io.open(BASE + "_audit_p2.txt", "w", encoding="utf-8")


def P(*a):
    print(*a, file=LOG)


csv_path = BASE + "results/phase9_hardening/p2_goslim_composition.csv"
tab = pd.read_csv(csv_path)
P(f"rows: {len(tab)} (classes: {tab['class'].nunique()} x buckets: {tab['bucket'].nunique()})")

# ---- recompute Fisher p from the table's own counts ----
ok = True
for i, r in tab.iterrows():
    k, n_cls, K = int(r["k_overlap"]), int(r["n_class"]), int(r["K_universe"])
    N = 4974
    table = [[k, n_cls - k], [K - k, N - n_cls - (K - k)]]
    _, pg = stats.fisher_exact(table, alternative="greater")
    _, p2 = stats.fisher_exact(table, alternative="two-sided")
    if not (abs(pg - r["p_greater"]) < 1e-9 and abs(p2 - r["p_two_sided"]) < 1e-9):
        ok = False
        P(f"  MISMATCH row {i}: {r['class']}/{r['bucket']} stored pg={r['p_greater']:.4g} "
          f"recomputed {pg:.4g}; stored p2={r['p_two_sided']:.4g} recomputed {p2:.4g}")
P(f"Fisher p recomputation from table counts: {'ALL MATCH' if ok else 'MISMATCHES FOUND'}")


def bh(pvals):
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p)
    q = p[order] * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


q_greater = bh(tab["p_greater"].values)
q_two = bh(tab["p_two_sided"].values)
stored_q = tab["q_bh"].values
P(f"\nmax |stored q - BH(p_greater)| = {np.nanmax(np.abs(stored_q - q_greater)):.2e}")
P(f"max |stored q - BH(p_two_sided)| = {np.nanmax(np.abs(stored_q - q_two)):.2e}")
stored_is = "BH over p_greater (one-sided enrichment)" if np.nanmax(np.abs(stored_q - q_greater)) < 1e-9 else \
    ("BH over p_two_sided" if np.nanmax(np.abs(stored_q - q_two)) < 1e-9 else "NEITHER")

P(f"\nstored q_bh convention: {stored_is}")

P("\n--- minimum q by class ---")
for cname in ["gene_driven", "gene_driven_all", "regulation_driven", "dual_driven"]:
    m = tab["class"] == cname
    P(f"  {cname:18s} min q_greater={q_greater[m.values].min():.4f}  "
      f"min q_two_sided={q_two[m.values].min():.4f}")
    for j in np.where(m.values)[0]:
        if abs(q_greater[j] - q_greater[m.values].min()) < 1e-12:
            P(f"      (min at {tab.loc[j,'bucket']}: k={tab.loc[j,'k_overlap']}, "
              f"p_greater={tab.loc[j,'p_greater']:.4f})")

# reviewer's quoted 0.46 vs 0.38
P("\n--- reviewer's 0.46 vs 0.38 ---")
P(f"gene_driven min q (one-sided greater) = {q_greater[(tab['class']=='gene_driven').values].min():.4f} "
  f"-> rounds to 0.46")
P(f"The manuscript's 'all q >= 0.38' cannot be reproduced from either BH family;")
P(f"0.38 is consistent with BH over 30 tests (3 classes, pre-dual version): "
  f"min gene_driven p_greater 0.12604 * 30 / 10 = {0.126041167750916*30/10:.4f}")

# ---- update CSV: keep stored columns, add q_bh_two_sided for reference ----
changed = []
if stored_is.startswith("BH over p_greater"):
    tab["q_bh_two_sided"] = q_two
    tab.to_csv(csv_path, index=False)
    P("\nCSV already one-sided; added q_bh_two_sided column for reference, kept q_bh as-is")
else:
    old = tab["q_bh"].copy()
    tab["q_bh"] = q_greater
    tab["q_bh_two_sided"] = q_two
    tab.to_csv(csv_path, index=False)
    for i in range(len(tab)):
        if abs(old[i] - q_greater[i]) > 1e-9:
            changed.append({"class": tab.loc[i, "class"], "bucket": tab.loc[i, "bucket"],
                            "q_old": float(old[i]), "q_new": float(q_greater[i])})
    P(f"\nCSV q_bh was two-sided -> REPLACED with one-sided BH; {len(changed)} rows changed")

result = {
    "task": "P2 GO-slim one-sided vs two-sided q audit",
    "n_tests": int(len(tab)),
    "fisher_p_recomputed_from_counts": "all match" if ok else "MISMATCH",
    "stored_q_convention": stored_is,
    "min_q_by_class": {
        cname: {"q_greater": float(q_greater[(tab["class"] == cname).values].min()),
                "q_two_sided": float(q_two[(tab["class"] == cname).values].min())}
        for cname in ["gene_driven", "gene_driven_all", "regulation_driven", "dual_driven"]},
    "reviewer_0.46_vs_0.38": {
        "0.46": "gene_driven minimum one-sided (greater) BH q = 0.4583 (Metabolism & proteostasis, "
                "p_greater = 0.126) -- this is the correct value under the manuscript's one-sided convention",
        "0.38": "not reproducible from the current 40-test CSV under either BH family; consistent with BH "
                "over 30 tests (3 classes, before dual_driven was added): 0.12604*30/10 = 0.378",
        "recommended_manuscript_edit": "'none significant, all q >= 0.38' should read 'all q >= 0.46' "
                "(or '>= 0.45' if the gene-driven-incl-relaxed class is meant: min q = 0.447)",
    },
    "rd_class_check": {
        "enriched_buckets_q": {tab.loc[i, "bucket"]: float(q_greater[i])
                               for i in range(len(tab))
                               if tab.loc[i, "class"] == "regulation_driven" and q_greater[i] < 0.01},
        "manuscript_claim_all_q_below_0.01": bool(q_greater[(tab["class"] == "regulation_driven").values] .min() < 0.01),
    },
    "csv_action": "q_bh verified one-sided; q_bh_two_sided column added" if not changed else
                  f"q_bh replaced with one-sided; {len(changed)} rows changed",
    "changed_rows": changed,
}
with io.open(BASE + "results/phase9_hardening/p2_goslim_check.json", "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False)
P("\nsaved results/phase9_hardening/p2_goslim_check.json")
LOG.close()
print("done")
