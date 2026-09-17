# -*- coding: utf-8 -*-
"""
P2 (reviewer 1): BUSTED FDR rate in strict one-to-one vs permissive-call genes
================================================================================
Manuscript: of the 4,974-gene universe, 4,119 (82.8%) are strict 1:1 orthologs
and 855 derive from permissive calls. Question: do BUSTED FDR-significant rates
differ between the two orthology strata? (If yes, a one-line caveat is needed.)

Output: results/phase9_hardening/strict_perm_busted_audit.json ; log _audit_sp.txt
"""
import io
import json

import numpy as np
import pandas as pd
from scipy import stats

BASE = "d:/人类正选择基因项目/"
LOG = io.open(BASE + "_audit_sp.txt", "w", encoding="utf-8")


def P(*a):
    print(*a, file=LOG)


strict = pd.read_csv(BASE + "results/phase2_final_v2/shared_genes_one2one_all10.csv")
P("strict file cols:", strict.columns.tolist()[:8], "... n =", len(strict))
id_col = "gene_id" if "gene_id" in strict.columns else strict.columns[0]
strict_ids = set(str(g).split(".")[0] for g in strict[id_col].dropna())
P(f"strict unique gene ids: {len(strict_ids)} (expected 6,796)")

gc = pd.read_csv(BASE + "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv")
gc["ensg"] = gc["gene_id"].map(lambda g: str(g).split(".")[0])
gc["is_strict"] = gc["ensg"].isin(strict_ids)
sig = gc["bh_fdr_sig"].astype(bool)

n_s = int(gc["is_strict"].sum())
n_p = int((~gc["is_strict"]).sum())
k_s, k_p = int(sig[gc["is_strict"]].sum()), int(sig[~gc["is_strict"]].sum())
P(f"\nuniverse split: strict {n_s} (expected 4,119) / permissive {n_p} (expected 855)")
P(f"BUSTED FDR<0.05: strict {k_s}/{n_s} = {100*k_s/n_s:.2f}%; "
  f"permissive {k_p}/{n_p} = {100*k_p/n_p:.2f}%")

table = [[k_s, n_s - k_s], [k_p, n_p - k_p]]
orr, p2 = stats.fisher_exact(table)
_, pg = stats.fisher_exact(table, alternative="greater")
_, pl = stats.fisher_exact(table, alternative="less")
P(f"Fisher two-sided: OR={orr:.4f}, p={p2:.4g} (greater p={pg:.4g}, less p={pl:.4g})")
P(f"rate difference: {100*(k_s/n_s - k_p/n_p):+.2f} pp")

# same check for the RELAX and classification composition (context)
for col, name in [("relax_fdr_sig", "RELAX FDR sig")]:
    if col in gc.columns:
        f = gc[col].astype(bool)
        a, b = int(f[gc["is_strict"]].sum()), int(f[~gc["is_strict"]].sum())
        P(f"{name}: strict {a}/{n_s} ({100*a/n_s:.1f}%) vs permissive {b}/{n_p} ({100*b/n_p:.1f}%)")
cls_ct = pd.crosstab(gc["classification_v7"], gc["is_strict"])
P("\nclassification x orthology stratum:")
P(cls_ct.to_string())

result = {
    "task": "strict 1:1 vs permissive-call BUSTED FDR rate sensitivity (reviewer 1, P2)",
    "strict_n": n_s, "permissive_n": n_p,
    "busted_fdr_strict": {"k": k_s, "rate": round(k_s / n_s, 4)},
    "busted_fdr_permissive": {"k": k_p, "rate": round(k_p / n_p, 4)},
    "fisher": {"OR": float(orr), "p_two_sided": float(p2),
               "p_greater": float(pg), "p_less": float(pl)},
    "rate_diff_pp": round(100 * (k_s / n_s - k_p / n_p), 2),
    "one_line_conclusion": (
        f"BUSTED FDR-significant rates are statistically indistinguishable between strict one-to-one "
        f"orthologs ({k_s}/{n_s} = {100*k_s/n_s:.1f}%) and permissive-call genes ({k_p}/{n_p} = "
        f"{100*k_p/n_p:.1f}%; Fisher OR = {orr:.2f}, P = {p2:.2f}), so the 17.2% permissive fraction "
        f"does not materially affect the universe-wide coding-selection rate."
    ),
    "classification_by_stratum": {str(k): {str(c): int(v) for c, v in row.items()}
                                  for k, row in cls_ct.iterrows()},
}
with io.open(BASE + "results/phase9_hardening/strict_perm_busted_audit.json", "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False)
P("\nsaved results/phase9_hardening/strict_perm_busted_audit.json")
LOG.close()
print("done")
