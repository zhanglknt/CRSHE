# -*- coding: utf-8 -*-
"""
P0-D audit (final): Selectome "OR = 2.94, P = 0.030" claim in manuscript v9 L78
================================================================================
Findings so far (see _audit_sel*.txt):
  1. The stored has_selectome flag (16 genes in universe) comes from
     results/phase4_hyphy/selectome_analysis/selectome_primate_genes_summary_v3.tsv
     (41 genes), whose family->gene mapping in extract_selectome_detailed_v3.py
     joins Selectome SQL internal tree ids (small integers) with ENSGT-derived
     family numbers -- two DIFFERENT id spaces. Numeric collision => arbitrary
     gene assignment (pps ∩ summary_v3 = 3 genes; none of the 8 manuscript
     Tier-1 genes is in the NHX-parsed primate list).
  2. The manuscript OR=2.94/P=0.030 reproduces from the stored flag as
     [[8,1258],[8,3700]] one-sided-greater Fisher.
This script recomputes everything from the raw Selectome v6 data (NHX-parsed
primate-ancestral positive selection, plus the REST any-branch list as
sensitivity) and writes selectome_audit.json.

Outputs: results/phase9_hardening/selectome_audit.json ; log _audit_sel_final.txt
"""
import io
import json

import numpy as np
import pandas as pd
from scipy import stats

BASE = "d:/人类正选择基因项目/"
LOG = io.open(BASE + "_audit_sel_final.txt", "w", encoding="utf-8")


def P(*a):
    print(*a, file=LOG)


def fisher_block(table, name):
    t = np.asarray(table, dtype=int)  # int for scipy.stats.contingency
    orr, p2 = stats.fisher_exact(t, alternative="two-sided")
    _, pg = stats.fisher_exact(t, alternative="greater")
    lo = hi = None
    try:
        r = stats.contingency.odds_ratio(t, kind="conditional")
        ci = r.confidence_interval(0.95)
        lo = 0.0 if not np.isfinite(ci.low) else float(ci.low)
        hi = float("inf") if not np.isfinite(ci.high) else float(ci.high)
    except Exception as e:  # noqa: BLE001
        P(f"  (conditional CI failed: {e})")
    out = {"name": name,
           "table": [[int(t[0, 0]), int(t[0, 1])], [int(t[1, 0]), int(t[1, 1])]],
           "OR": float(orr) if np.isfinite(orr) else None,
           "p_two_sided": float(p2), "p_greater": float(pg),
           "OR_ci95_low": lo, "OR_ci95_high": hi}
    P(f"  {name}: table={out['table']}")
    P(f"    OR={out['OR']}, p_two={p2:.6g}, p_greater={pg:.6g}, CI95=({lo}, {hi})")
    return out


# ---------------------------------------------------------------- load
gc = pd.read_csv(BASE + "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv")
n_uni = len(gc)
sig = gc["bh_fdr_sig"].astype(bool)
gmap = pd.read_csv(BASE + "data/gene_id_to_symbol_gencode47.csv")
id2sym = {str(g).split(".")[0]: s for g, s in zip(gmap["gene_id"], gmap["gene_symbol"])}
gc["sym"] = gc["gene_id"].map(lambda g: id2sym.get(str(g).split(".")[0]))
P(f"universe: {n_uni}; BUSTED FDR<0.05: {int(sig.sum())} ({100 * sig.mean():.2f}%)")

pps = pd.read_csv(BASE + "data/selectome/selectome_primate_positive_selection.tsv", sep="\t")
hps = pd.read_csv(BASE + "data/selectome/selectome_human_positive_selection.tsv", sep="\t")
s3 = pd.read_csv(BASE + "results/phase4_hyphy/selectome_analysis/selectome_primate_genes_summary_v3.tsv", sep="\t")
pps_ids = set(str(g).split(".")[0] for g in pps["gene_id"].dropna())
hps_ids = set(str(g).split(".")[0] for g in hps["gene_id"].dropna())
s3_ids = set(str(g).split(".")[0] for g in s3["gene_id"].dropna())
stored_ids = set(gc.loc[gc["has_selectome"].astype(bool), "gene_id"].map(lambda g: str(g).split(".")[0]))
P(f"\npps (NHX primate-ancestral, raw v6): {len(pps_ids)} unique genes")
P(f"hps (REST any-branch selection): {len(hps_ids)} unique genes")
P(f"summary_v3 (SQL dump, suspect mapping): {len(s3_ids)} unique genes")
P(f"pps ⊆ hps: {pps_ids <= hps_ids}")
P(f"pps ∩ summary_v3 = {len(pps_ids & s3_ids)}; stored ∩ pps = {len(stored_ids & pps_ids)}")

gdi = gc["classification_v7"].isin(["gene-driven", "gene-driven (relaxed)"])


def run_flag(flag, label):
    """Full analysis for a given Selectome flag: class breakdown, dual genes, Fisher."""
    P(f"\n================ {label} ================")
    n_in = int(flag.sum())
    by = {}
    for c in ["gene-driven", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "neutral"]:
        m = gc["classification_v7"] == c
        by[c] = [int(flag[m].sum()), int(m.sum())]
        P(f"  {c:24s} {by[c][0]}/{by[c][1]}")
    both = gc[flag & sig]
    P(f"  BUSTED∩Selectome (dual-method): {len(both)} genes: "
      f"{sorted(str(s) for s in both['sym'])}")
    k1, n1 = by["gene-driven"][0] + by["gene-driven (relaxed)"][0], \
        by["gene-driven"][1] + by["gene-driven (relaxed)"][1]
    k0 = n_in - k1
    n0 = n_uni - n1
    t = fisher_block([[k1, n1 - k1], [k0, n0 - k0]],
                     f"{label}: GD incl. relaxed vs rest")
    return {"n_in_universe": n_in, "by_class": by,
            "dual_method_genes": sorted(str(s) for s in both["sym"]), "fisher": t}


# stored flag (manuscript basis)
res_stored = run_flag(gc["has_selectome"].astype(bool), "STORED has_selectome (summary_v3, suspect)")
# corrected: raw NHX primate-ancestral list
flag_pps = gc["gene_id"].map(lambda g: str(g).split(".")[0] in pps_ids)
res_pps = run_flag(flag_pps, "CORRECTED pps (NHX primate-ancestral, raw Selectome v6)")
# sensitivity: REST any-branch list
flag_hps = gc["gene_id"].map(lambda g: str(g).split(".")[0] in hps_ids)
res_hps = run_flag(flag_hps, "SENSITIVITY hps (REST any-branch selection)")

# reviewer hand-calc tables for the record
P("\n================ reviewer hand-calculations ================")
hand = {
    "eight_vs_zero_rest": fisher_block([[8, 1258], [0, 3708]],
                                       "8 dual genes vs 0 in rest (reviewer calc 1)"),
    "eight_vs_rd_only": fisher_block([[8, 1258], [0, 293]],
                                     "8 in GD+GDr vs 0/293 RD (reviewer calc 2)"),
}

# ---------------------------------------------------------------- verdict
fA = res_stored["fisher"]
P("\n================ VERDICT ================")
P(f"1. OR=2.94/P=0.030 IS reproducible from the stored 16-gene flag: "
  f"[[8,1258],[8,3700]] OR={fA['OR']:.4f}, p_greater={fA['p_greater']:.4f} "
  f"(one-sided; the manuscript's P=0.030), p_two_sided={fA['p_two_sided']:.4f}.")
P("   It is NOT a copy-paste of the hCONDEL OR=2.94 (different 2x2), but a numerical coincidence.")
P("2. HOWEVER the stored flag itself is invalid: it derives from a family->gene map that")
P("   joins Selectome SQL internal ids with ENSGT family numbers (different id spaces);")
P("   none of the 8 manuscript Tier-1 genes appears in the raw NHX-parsed primate list.")
P(f"3. Corrected (pps, raw v6 NHX): {res_pps['n_in_universe']} genes in universe, "
  f"{len(res_pps['dual_method_genes'])} dual-method; Fisher GD vs rest "
  f"OR={res_pps['fisher']['OR']:.3f}, p_greater={res_pps['fisher']['p_greater']:.3f}, "
  f"p_two={res_pps['fisher']['p_two_sided']:.3f} -> NOT significant.")

result = {
    "claim_audited": "manuscript_english_v9.md L78 & Fig 5d: 'Eight genes detected by both BUSTED and "
                     "Selectome v6 ... Fisher OR = 2.94, P = 0.030'",
    "universe_n": int(n_uni),
    "busted_fdr_sig": int(sig.sum()),
    "selectome_sources": {
        "pps_nhx_primate_ancestral": {"file": "data/selectome/selectome_primate_positive_selection.tsv",
                                      "unique_genes": len(pps_ids)},
        "hps_rest_any_branch": {"file": "data/selectome/selectome_human_positive_selection.tsv",
                                "unique_genes": len(hps_ids)},
        "summary_v3_sql": {"file": "results/phase4_hyphy/selectome_analysis/selectome_primate_genes_summary_v3.tsv",
                           "unique_genes": len(s3_ids),
                           "status": "SUSPECT: family->gene mapping joins Selectome SQL internal tree ids "
                                     "(1-6 digit integers) with ENSGT-derived family numbers -- different id "
                                     "spaces; only 3/94 genes shared with the NHX-parsed list"},
    },
    "stored_flag_analysis": res_stored,
    "corrected_pps_analysis": res_pps,
    "sensitivity_hps_any_branch": res_hps,
    "reviewer_hand_calculations": hand,
    "verdict": {
        "or_2.94_p_0.030_reproduced": True,
        "reproduced_from": "stored has_selectome (16-gene summary_v3 flag); 2x2 = [[8,1258],[8,3708-8=3700]]",
        "reported_P_is": "one-sided greater (p=0.0299); two-sided p=0.0389",
        "hcondel_OR_identity": "numerical coincidence, not copy-paste (hCONDEL test is a different 2x2)",
        "but_stored_flag_invalid": True,
        "invalidity_evidence": [
            "extract_selectome_detailed_v3.py maps SQL branch-table 'id' (Selectome internal tree ids, "
            "small integers) to ENSGT-derived family numbers from the NHX/REST files -- different id spaces",
            "pps(NHX, 94 genes) ∩ summary_v3(41 genes) = 3 genes",
            "0 of the 8 manuscript Tier-1 genes (E4F1,TAPT1,GLRX,RBM22,C2orf74,LDB1,YPEL5,TXNDC16) appear "
            "in the raw NHX primate list; only 2 (E4F1,TXNDC16) appear even in the 1,332-gene any-branch "
            "REST list",
        ],
        "corrected_numbers": {
            "selectome_primate_genes_in_universe": res_pps["n_in_universe"],
            "dual_method_genes_busted_and_selectome": res_pps["dual_method_genes"],
            "fisher_gd_vs_rest": {"OR": res_pps["fisher"]["OR"],
                                  "p_one_sided_greater": res_pps["fisher"]["p_greater"],
                                  "p_two_sided": res_pps["fisher"]["p_two_sided"],
                                  "ci95": [res_pps["fisher"]["OR_ci95_low"],
                                           res_pps["fisher"]["OR_ci95_high"]]},
        },
        "manuscript_impact": "The Selectome enrichment claim (OR=2.94, P=0.030) and the 'eight Tier-1 "
                             "highest-confidence genes' list do not survive recomputation from the raw "
                             "Selectome v6 data. With the corrected primate-ancestral list the enrichment "
                             "is null (OR~1.2, p~0.2). Recommend either removing the Selectome enrichment "
                             "sentence + Fig 5d Fisher annotation + Table S3-1 list, or rebuilding the "
                             "Selectome cross-reference from the NHX/REST data with verified ENSG matching "
                             "and reporting the corrected (null) result.",
    },
}
with io.open(BASE + "results/phase9_hardening/selectome_audit.json", "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False, default=float)
P("\nsaved results/phase9_hardening/selectome_audit.json")
LOG.close()
print("done")
