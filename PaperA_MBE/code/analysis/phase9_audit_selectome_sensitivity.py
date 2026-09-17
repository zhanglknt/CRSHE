# -*- coding: utf-8 -*-
"""
P0-D follow-up: sensitivity reclassification with the corrected Selectome flag
================================================================================
The stored has_selectome (16 genes, from the invalid summary_v3 mapping) feeds
the GDS Selectome component (weight 0.15). Replace it with the corrected NHX
primate-ancestral flag (24 genes in universe), keep every other component and
the v9 classification rule identical, and quantify the impact.

GDS v7 = 0.35*gds_p_pct_resid + 0.30*gds_lrt_pct_resid + 0.20*gds_relax_pct
          + 0.15*selectome_flag
RDS v7 = 0.30*doan + 0.25*rds_tau + 0.25*rds_nc + 0.20*rds_brain   (unchanged)
Rules (phase7g_v7_classification.py classify(), threshold 0.15):
  GD:  gds > rds + 0.15 & bh_fdr_sig
  RD:  rds > gds + 0.15 & ~bh_fdr_sig
  dual: rds > gds + 0.15 & bh_fdr_sig & gds > 0.3
  neutral otherwise; gene-driven & relax_relaxed -> gene-driven (relaxed)

Outputs: results/phase9_hardening/selectome_sensitivity_reclass.json
         results/phase9_hardening/selectome_sensitivity_flips.csv
         log: verification/_audit_sel_sens.txt
"""
import io
import json

import numpy as np
import pandas as pd
from scipy import stats

BASE = "d:/人类正选择基因项目/"
LOG = io.open(BASE + "results/phase9_hardening/verification/_audit_sel_sens.txt", "w", encoding="utf-8")


def P(*a):
    print(*a, file=LOG)


gc = pd.read_csv(BASE + "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv")
strip = lambda s: str(s).split(".")[0]
gc["ensg"] = gc["gene_id"].map(strip)
gmap = pd.read_csv(BASE + "data/gene_id_to_symbol_gencode47.csv")
id2sym = {strip(g): s for g, s in zip(gmap["gene_id"], gmap["gene_symbol"])}
gc["sym"] = gc["ensg"].map(id2sym)

# ---- corrected flag from raw NHX primate-ancestral list ----
pps = pd.read_csv(BASE + "data/selectome/selectome_primate_positive_selection.tsv", sep="\t")
pps_ids = set(strip(g) for g in pps["gene_id"].dropna())
stored = gc["has_selectome"].astype(bool).values
corrected = gc["ensg"].isin(pps_ids).values
P(f"stored flag: {int(stored.sum())} genes; corrected NHX flag: {int(corrected.sum())} genes")
old_only = set(gc.loc[stored & ~corrected, "sym"])
new_only = set(gc.loc[corrected & ~stored, "sym"])
both = set(gc.loc[stored & corrected, "sym"])
P(f"overlap: {len(both)} ({sorted(both)}); old-only: {len(old_only)}; new-only: {len(new_only)}")
P(f"symmetric difference: {len(old_only) + len(new_only)} genes")
P(f"old-only: {sorted(s for s in old_only if isinstance(s, str))}")
P(f"new-only: {sorted(s for s in new_only if isinstance(s, str))}")

# ---- verify GDS formula against stored gds_v7 ----
gds_stored = gc["gds_v7"].values
gds_recomputed = (0.35 * gc["gds_p_pct_resid"] + 0.30 * gc["gds_lrt_pct_resid"]
                  + 0.20 * gc["gds_relax_pct"] + 0.15 * gc["gds_selectome"]).values
P(f"\nGDS formula check: max|stored - recomputed| = {np.max(np.abs(gds_stored - gds_recomputed)):.2e}")

# ---- scores with corrected flag ----
rds = gc["rds_v7"].values  # RDS unchanged
gds_new = (0.35 * gc["gds_p_pct_resid"] + 0.30 * gc["gds_lrt_pct_resid"]
           + 0.20 * gc["gds_relax_pct"] + 0.15 * corrected.astype(float)).values
sig = gc["bh_fdr_sig"].astype(bool).values
relaxed = gc["relax_relaxed"].astype(bool).values
THRESH = 0.15


def classify(gds):
    gd = (gds > rds + THRESH) & sig
    rd = (rds > gds + THRESH) & (~sig)
    dual = (rds > gds + THRESH) & sig & (gds > 0.3)
    cls = np.where(gd, "gene-driven",
                   np.where(rd, "regulation-driven",
                            np.where(dual, "dual-driven", "neutral")))
    cls = np.where(relaxed & (cls == "gene-driven"), "gene-driven (relaxed)", cls)
    return cls


cls_old = classify(gds_stored)
P(f"replication of stored classification: {(cls_old == gc['classification_v7']).mean():.4f}")
cls_new = classify(gds_new)

# ---- flips ----
flip_mask = cls_old != cls_new
flips = pd.DataFrame({
    "gene_id": gc.loc[flip_mask, "gene_id"],
    "gene_symbol": gc.loc[flip_mask, "sym"],
    "old_class": cls_old[flip_mask],
    "new_class": cls_new[flip_mask],
    "stored_selectome_flag": stored[flip_mask],
    "corrected_selectome_flag": corrected[flip_mask],
    "gds_old": gds_stored[flip_mask].round(4),
    "gds_new": gds_new[flip_mask].round(4),
    "rds_v7": rds[flip_mask].round(4),
    "bh_fdr_sig": sig[flip_mask],
})
flips.to_csv(BASE + "results/phase9_hardening/selectome_sensitivity_flips.csv", index=False)
P(f"\nclassification flips: {int(flip_mask.sum())} genes")
if flip_mask.sum():
    P(flips.to_string(index=False))

cnt_old = pd.Series(cls_old).value_counts().to_dict()
cnt_new = pd.Series(cls_new).value_counts().to_dict()
P("\nclass counts old -> new:")
for c in ["gene-driven", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "neutral"]:
    P(f"  {c:24s} {cnt_old.get(c, 0):5d} -> {cnt_new.get(c, 0):5d}")

# ---- class-level robustness checks (only meaningful if classes changed) ----
robust = {}
if flip_mask.sum() > 0:
    # (1) hCONDEL RD enrichment (manuscript: RD OR = 2.94, CI 1.92-4.51, p = 6.7e-6)
    hm = pd.read_csv(BASE + "results/phase8_tissue_analysis/hcondel_gene_mapping.csv")
    hc_ids = set(strip(g) for g in hm["gene_id"])
    hc = gc["ensg"].isin(hc_ids).values
    for tag, cls in [("old", cls_old), ("new", cls_new)]:
        rd = cls == "regulation-driven"
        t = [[int((rd & hc).sum()), int((rd & ~hc).sum())],
             [int((~rd & hc).sum()), int((~rd & ~hc).sum())]]
        orr, p = stats.fisher_exact(t)
        r = stats.contingency.odds_ratio(np.asarray(t), kind="conditional")
        ci = r.confidence_interval(0.95)
        robust[f"hcondel_RD_enrichment_{tag}"] = {
            "table": t, "OR": float(orr), "p": float(p),
            "ci95": [float(ci.low), float(ci.high)]}
        P(f"\nhCONDEL RD enrichment ({tag}): table={t}, OR={orr:.3f}, p={p:.3g}")
    # (2)/(3) SynGO 11 terms & GD zero FDR terms: flip count by class as proxy
    for tag, cls in [("old", cls_old), ("new", cls_new)]:
        rd = cls == "regulation-driven"
        gd_all = np.isin(cls, ["gene-driven", "gene-driven (relaxed)"])
        robust[f"class_sizes_{tag}"] = {"RD": int(rd.sum()), "GD_all": int(gd_all.sum())}
    rd_flips = int(((cls_old == "regulation-driven") ^ (cls_new == "regulation-driven")).sum())
    gd_flips = int((np.isin(cls_old, ["gene-driven", "gene-driven (relaxed)"])
                    ^ np.isin(cls_new, ["gene-driven", "gene-driven (relaxed)"])).sum())
    robust["RD_class_flips"] = rd_flips
    robust["GD_class_flips"] = gd_flips
    P(f"\nRD class flips: {rd_flips}; GD(incl.rel) class flips: {gd_flips}")

# tier1 universe-enrichment margins (robustness of 'GD zero FDR terms' /
# 'RD SynGO 11 terms' against the class change) -- read existing tables
q_margins = {}
for f in ["gene_driven_GO_BP", "gene_driven_SynGO", "gene_driven_all_GO_BP",
          "gene_driven_all_SynGO", "regulation_driven_GO_BP", "regulation_driven_SynGO"]:
    d = pd.read_csv(BASE + f"results/phase9_hardening/tier1_universe_enrichment_{f}.csv")
    q_margins[f] = {"n_terms": int(len(d)),
                    "min_p_bh": float(d["p_bh"].min()),
                    "n_p_bh_lt_0.05": int((d["p_bh"] < 0.05).sum()),
                    "n_p_lt_0.05": int((d["p"] < 0.05).sum())}
    P(f"tier1 {f}: min p_bh={d['p_bh'].min():.4f}, FDR<0.05 terms={int((d['p_bh'] < 0.05).sum())}")
robust["tier1_q_margins"] = q_margins

result = {
    "task": "Selectome-flag sensitivity reclassification (P0-D follow-up)",
    "stored_flag_n": int(stored.sum()),
    "corrected_nhx_flag_n": int(corrected.sum()),
    "flag_overlap": len(both),
    "flag_symmetric_difference": len(old_only) + len(new_only),
    "flag_old_only_genes": sorted(str(s) for s in old_only),
    "flag_new_only_genes": sorted(str(s) for s in new_only),
    "gds_formula_max_abs_error": float(np.max(np.abs(gds_stored - gds_recomputed))),
    "classification_replication_agreement": float((cls_old == gc["classification_v7"]).mean()),
    "n_flips": int(flip_mask.sum()),
    "flips_csv": "results/phase9_hardening/selectome_sensitivity_flips.csv",
    "class_counts_old": {k: int(v) for k, v in cnt_old.items()},
    "class_counts_new": {k: int(v) for k, v in cnt_new.items()},
    "robustness_checks": robust,
}
with io.open(BASE + "results/phase9_hardening/selectome_sensitivity_reclass.json", "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False, default=float)
P("\nsaved selectome_sensitivity_reclass.json + selectome_sensitivity_flips.csv")
LOG.close()
print("done")
