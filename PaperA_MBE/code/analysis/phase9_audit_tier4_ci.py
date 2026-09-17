# -*- coding: utf-8 -*-
"""
Tier 4: rate differences vs genome-wide BUSTED rate with 95% CI (Newcombe)
================================================================================
For HAR-proximate (476), hCONDEL-proximate (183), and combined (630) subsets:
BUSTED FDR-significant rate vs the genome-wide rate (1,690/4,974 = 34.0%).
CI method: Newcombe hybrid score interval (Wilson intervals for each proportion,
combined as (L1 - U0, U1 - L0)) — appropriate for small-to-moderate n and does
not require the normal approximation of the difference. Wald CI reported for
comparison.

Output: updated tier4_empirical_disjointness.{csv,json} ; log _audit_t4.txt
"""
import io
import json
import math

import numpy as np
import pandas as pd
from scipy import stats

BASE = "d:/人类正选择基因项目/"
LOG = io.open(BASE + "_audit_t4.txt", "w", encoding="utf-8")


def P(*a):
    print(*a, file=LOG)


def wilson(k, n, z=1.959963984540054):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def newcombe(k1, n1, k0, n0):
    l1, u1 = wilson(k1, n1)
    l0, u0 = wilson(k0, n0)
    return l1 - u0, u1 - l0


def wald(k1, n1, k0, n0):
    p1, p0 = k1 / n1, k0 / n0
    se = math.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    return p1 - p0 - 1.959963984540054 * se, p1 - p0 + 1.959963984540054 * se


csv_path = BASE + "results/phase9_hardening/tier4_empirical_disjointness.csv"
tab = pd.read_csv(csv_path)

# genome-wide reference from the universe file (independent check of 34.0%)
gc = pd.read_csv(BASE + "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv")
k0, n0 = int(gc["bh_fdr_sig"].sum()), len(gc)
P(f"genome-wide BUSTED FDR<0.05: {k0}/{n0} = {100*k0/n0:.2f}%")

rows = []
for _, r in tab.iterrows():
    k1, n1 = int(r["busted_sig"]), int(r["n"])
    assert abs(k1 / n1 - r["busted_rate"]) < 5e-4, r
    ln, hn = newcombe(k1, n1, k0, n0)
    lw, hw = wald(k1, n1, k0, n0)
    _, pfish = stats.fisher_exact([[k1, n1 - k1], [k0, n0 - k0]])
    d = k1 / n1 - k0 / n0
    rows.append({
        "subset": r["subset"], "n": n1, "busted_sig": k1,
        "subset_rate": round(k1 / n1, 4),
        "genome_rate": round(k0 / n0, 4),
        "rate_diff": round(d, 4),
        "rate_diff_ci95_newcombe_low": round(ln, 4),
        "rate_diff_ci95_newcombe_high": round(hn, 4),
        "rate_diff_ci95_wald_low": round(lw, 4),
        "rate_diff_ci95_wald_high": round(hw, 4),
        "fisher_p_vs_genome": pfish,
    })
    P(f"{r['subset']}: {k1}/{n1} = {100*k1/n1:.2f}% vs genome {100*k0/n0:.2f}% "
      f"-> diff {100*d:+.2f} pp, Newcombe 95% CI ({100*ln:+.2f}, {100*hn:+.2f}) "
      f"[Wald ({100*lw:+.2f}, {100*hw:+.2f})], Fisher p={pfish:.3f}")

diffs = pd.DataFrame(rows)
diffs.to_csv(BASE + "results/phase9_hardening/tier4_rate_differences.csv", index=False)

# append to the existing CSV as extra columns
for col in ["rate_diff", "rate_diff_ci95_newcombe_low", "rate_diff_ci95_newcombe_high",
            "rate_diff_ci95_wald_low", "rate_diff_ci95_wald_high"]:
    tab[col] = [r[col] for r in rows]
tab.to_csv(csv_path, index=False)

# update JSON
jpath = BASE + "results/phase9_hardening/tier4_empirical_disjointness.json"
with io.open(jpath, encoding="utf-8") as f:
    j = json.load(f)
j["rate_difference_vs_genome"] = {
    "genome_wide": {"k": k0, "n": n0, "rate": round(k0 / n0, 4)},
    "method": "Newcombe hybrid score interval (Wilson per proportion, combined L1-U0 / U1-L0); "
              "Wald shown for comparison",
    "subsets": {r["subset"]: {kk: (vv if not isinstance(vv, float) else round(vv, 6))
                              for kk, vv in r.items()} for r in rows},
}
with io.open(jpath, "w", encoding="utf-8") as f:
    json.dump(j, f, indent=2, ensure_ascii=False)

P("\nupdated tier4_empirical_disjointness.csv/.json + tier4_rate_differences.csv")
LOG.close()
print("done")
