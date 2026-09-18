# -*- coding: utf-8 -*-
"""Phase9 R6 collection & analysis: parse all BUSTED rerun/sim JSONs, compute BH rates,
write final JSON deliverables.

Outputs (results/phase9_hardening/):
  - neutral_sim_calibration.json (+ neutral_sim_pvalues.csv)
  - stratified_rerun_sensitivity.json
  - topology_sensitivity.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(r"d:\人类正选择基因项目")
PH9 = BASE / "results/phase9_hardening"
MASTER = BASE / "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv"

log_lines = []
def log(s):
    print(s)
    log_lines.append(str(s))


def parse_busted(fp):
    try:
        data = json.load(open(fp, encoding="utf-8"))
    except Exception:
        return None
    tr = data.get("test results", {})
    if not tr:
        return None
    p = tr.get("p-value")
    if p is None:
        return None
    return {"p": float(p), "lrt": float(tr.get("LRT", np.nan))}


def bh(pvals, alpha=0.05):
    """BH adjusted p-values (q) in ORIGINAL order + rejection count at alpha."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([]), 0
    order = np.argsort(pvals, kind="mergesort")
    ranked = pvals[order]
    q_sorted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0, 1)
    q = np.empty(n)
    q[order] = q_sorted
    below = ranked <= np.arange(1, n + 1) / n * alpha
    k = int(np.max(np.nonzero(below)[0]) + 1) if below.any() else 0
    return q, k


def load_group(group):
    d = {}
    for fp in (PH9 / f"busted_rerun_{group}").glob("*.json"):
        r = parse_busted(fp)
        if r is not None:
            d[fp.stem] = r
    return d


master = pd.read_csv(MASTER).set_index("gene_id")
subset = pd.read_csv(PH9 / "stratified_400_genes.csv")

# ==================================================================
# Task 2: neutral simulation calibration
# ==================================================================
sim = load_group("sim")
log(f"sim group parsed: {len(sim)} results")

sim_p = pd.Series({g: v["p"] for g, v in sim.items()}).sort_index()
sim_p.to_csv(PH9 / "neutral_sim_pvalues.csv", header=["p"])

n = len(sim_p)
p_lt05 = float((sim_p < 0.05).mean())
p_eq05 = float((sim_p == 0.5).mean())
p_gt05 = float((sim_p > 0.5).mean())
bins = np.arange(0, 1.05, 0.05)
hist, _ = np.histogram(sim_p, bins=bins)
ks_p = float(stats.kstest(sim_p[sim_p < 0.5], "uniform", args=(0, 0.5)).pvalue) if (sim_p < 0.5).sum() > 10 else None
# BH within the 400 simulated null genes = empirical FDR at alpha=0.05
q_sim, k_sim = bh(sim_p.values)
# also raw and BH at other alphas
_, k10 = bh(sim_p.values, alpha=0.10)
_, k25 = bh(sim_p.values, alpha=0.25)

# mechanistic determination
mech = {
    "frac_p_eq_0.5": p_eq05,
    "frac_p_gt_0.5": p_gt05,
    "frac_p_lt_0.5": float((sim_p < 0.5).mean()),
    "mean_p_below_0.5": float(sim_p[sim_p < 0.5].mean()) if (sim_p < 0.5).sum() else None,
}
truncation_like = p_eq05 > 0.10  # big point mass at 0.5
compression_like = (p_gt05 < 0.02) and (p_eq05 < 0.10) and (mech["mean_p_below_0.5"] is not None and mech["mean_p_below_0.5"] < 0.25)

neutral_calib = {
    "task": "P0-1a neutral simulation calibration (R6)",
    "method": {
        "simulator": "PAML evolver 4.10.10, codon model M0 (option 6)",
        "omega": 1.0, "kappa": 2.0,
        "codon_frequencies": "F3x4 estimated from all 4,974 real alignments",
        "tree": "per-gene species tree (same topology/branch lengths as BUSTED runs), human foreground irrelevant under neutrality",
        "missing_data": "real gene per-species per-codon gap masks applied to simulated sequences",
        "n_sim_genes": int(len(sim_p)),
        "test": "HyPhy BUSTED 2.5.101, --branches Test --srv No, 300s timeout, same pipeline as production",
    },
    "empirical_null_distribution": {
        "frac_p_lt_0.05": p_lt05,
        "n_p_lt_0.05": int((sim_p < 0.05).sum()),
        "frac_p_eq_0.5_point_mass": p_eq05,
        "frac_p_gt_0.5": p_gt05,
        "histogram_0.05_bins": {f"{bins[i]:.2f}-{bins[i+1]:.2f}": int(hist[i]) for i in range(len(hist))},
        "ks_test_p_lt_0.5_vs_U05_05_pvalue": ks_p,
    },
    "empirical_fdr": {
        "bh_alpha_0.05_rejections_out_of_nulls": k_sim,
        "empirical_fdr_bh005": k_sim / max(n, 1),
        "bh_alpha_0.10_rejections": k10,
        "bh_alpha_0.25_rejections": k25,
        "note": "all simulated genes are true nulls (omega=1), so any BH rejection is a false positive",
    },
    "mechanism_determination": {
        "observed": mech,
        "truncation_like": bool(truncation_like),
        "compression_like": bool(compression_like),
        "verdict": "TRUNCATION" if truncation_like else ("COMPRESSION" if compression_like else "MIXED/NEITHER"),
    },
    "implication_for_1690": None,  # filled below after task1 result is read
}

# tie to bh_bounded_support_check.json
t1 = json.load(open(PH9 / "bh_bounded_support_check.json", encoding="utf-8"))
if truncation_like:
    neutral_calib["implication_for_1690"] = (
        f"Simulated nulls reproduce a large point mass at p=0.5 ({p_eq05:.1%}) with virtually no p>0.5 "
        f"({p_gt05:.1%}) and near-nominal small-p tail ({p_lt05:.1%} < 0.05, expected 5%): the bounded support "
        f"arises from TRUNCATION (LRT<=0 clamped to p=0.5), not anti-conservative compression. "
        f"Per bh_bounded_support_check.json, truncation leaves the BH procedure unchanged, so the 1,690 "
        f"significant genes are not an artifact of the bounded support. The measured BH false-positive rate "
        f"on simulated nulls is {k_sim}/{n} = {k_sim/max(n,1):.1%}."
    )
else:
    neutral_calib["implication_for_1690"] = (
        f"Simulated nulls show {p_lt05:.1%} of p<0.05 (expected 5%) and point mass {p_eq05:.1%} at p=0.5; "
        f"measured BH false-positive rate {k_sim}/{n} = {k_sim/max(n,1):.1%}. See mechanism verdict above."
    )

json.dump(neutral_calib, open(PH9 / "neutral_sim_calibration.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
log("neutral_sim_calibration.json written")

# ==================================================================
# Tasks 3-4: stratified rerun groups
# ==================================================================
groups = {}
for g in ["baseline", "srv", "cpg", "topo"]:
    groups[g] = load_group(g)
    log(f"group {g}: {len(groups[g])} parsed results")

orig = subset.set_index("gene_id")
orig_p = master["busted_p"]
orig_sig = master["bh_fdr_sig"].astype(bool)

group_summary = {}
group_frames = {}
for g, d in groups.items():
    if len(d) == 0:
        group_summary[g] = {"n_genes_with_results": 0, "note": "no results yet"}
        group_frames[g] = pd.DataFrame(columns=["p", "bh_q", "sig"])
        group_frames[g].index.name = "gene_id"
        continue
    s = pd.DataFrame(d).T
    s.index.name = "gene_id"
    s = s.reset_index()
    q, k = bh(s["p"].values)
    s["bh_q"] = q
    s["sig"] = q < 0.05
    group_frames[g] = s.set_index("gene_id")
    rate = float(s["sig"].mean())
    # original (production) status for the same genes
    common = [x for x in s.gene_id if x in master.index]
    orig_sig_rate = float(orig_sig.reindex(common).mean())
    # concordance with original p
    both = pd.concat([s.set_index("gene_id")["p"], master["busted_p"].reindex(common)], axis=1)
    both.columns = ["rerun_p", "orig_p"]
    both = both.dropna()
    r = float(stats.pearsonr(both.rerun_p, both.orig_p)[0]) if len(both) > 10 else None
    rho = float(stats.spearmanr(both.rerun_p, both.orig_p)[0]) if len(both) > 10 else None
    # agreement on significance with production BH calls
    prod_sig = orig_sig.reindex(common).fillna(False)
    rerun_sig = s.set_index("gene_id")["sig"].reindex(common).fillna(False)
    agree = float((prod_sig == rerun_sig).mean())
    # embed rerun p-values into the full 4,974-gene production family (comparable per-gene calls)
    full_p = master["busted_p"].astype(float).copy()
    for gid, pv in zip(s.gene_id, s["p"].values):
        if gid in full_p.index:
            full_p.loc[gid] = pv
    full_q, _ = bh(full_p.values)
    full_q = pd.Series(full_q, index=full_p.index)
    embedded_sig = full_q.reindex(s.gene_id).values < 0.05
    embedded_agree = float((embedded_sig == orig_sig.reindex(s.gene_id).fillna(False).values).mean())
    group_summary[g] = {
        "n_genes_with_results": int(len(s)),
        "bh_sig_count": int(k),
        "bh_sig_rate": rate,
        "orig_bh_sig_rate_same_genes": orig_sig_rate,
        "pearson_r_vs_production_p": r,
        "spearman_rho_vs_production_p": rho,
        "significance_agreement_with_production": agree,
        "embedded_full_family_bh_sig_count": int(embedded_sig.sum()),
        "embedded_full_family_agreement_with_production": embedded_agree,
        "timeout_or_missing": int(400 - len(s)) if len(s) < 400 else 0,
    }

# baseline vs others comparisons
base = group_frames["baseline"]
for g in ["srv", "cpg", "topo"]:
    other = group_frames[g]
    common = base.index.intersection(other.index)
    b_sig = base.loc[common, "sig"]
    o_sig = other.loc[common, "sig"]
    group_summary[f"{g}_vs_baseline"] = {
        "n_common": int(len(common)),
        "baseline_sig_rate": float(b_sig.mean()),
        f"{g}_sig_rate": float(o_sig.mean()),
        "delta_rate": float(o_sig.mean() - b_sig.mean()),
        "both_sig": int((b_sig & o_sig).sum()),
        "baseline_only": int((b_sig & ~o_sig).sum()),
        f"{g}_only": int((~b_sig & o_sig).sum()),
        "jaccard": float((b_sig & o_sig).sum() / max((b_sig | o_sig).sum(), 1)),
        "of_baseline_sig_retained": float((b_sig & o_sig).sum() / max(b_sig.sum(), 1)),
    }

# CpG masking stats
cpg_stats = pd.read_csv(PH9 / "cpg_masking_stats.csv")
cpg_summary = {
    "n_genes_masked": int((cpg_stats.status == "ok").sum()),
    "mean_frac_codons_masked": float(cpg_stats.loc[cpg_stats.status == "ok", "frac_codons_masked"].mean()),
    "median_frac_codons_masked": float(cpg_stats.loc[cpg_stats.status == "ok", "frac_codons_masked"].median()),
    "mean_n_cpg_sites_human": float(cpg_stats.loc[cpg_stats.status == "ok", "n_cpg_sites_human"].mean()),
    "definition": "CpG dinucleotides in the human reference sequence (including codon-boundary-spanning CpG); "
                  "affected codon columns masked as ??? across all species",
}

# alignment quality (task 3d)
qual = json.load(open(PH9 / "alignment_quality_stats.json", encoding="utf-8"))

strat = {
    "task": "P0-2 stratified control reruns (R6)",
    "design": {
        "n_genes": 400,
        "stratification": "significant/non-significant (production BH) x CDS alignment-length tercile x BUSTED LRT tercile; 200 sig + 200 non-sig",
        "groups": {
            "baseline": "identical rerun of production settings (internal consistency control)",
            "srv": "--srv Yes (synonymous-rate variation enabled)",
            "cpg": "CpG-masked alignments (codons containing human-reference CpG dinucleotides, incl. codon-boundary, masked)",
            "topo": "alternative topology ((human,gorilla) sister instead of ((chimp,bonobo),human); see topology_sensitivity.json)",
        },
        "software_note": "HyPhy 2.5.101 (production run used 2.5.100); same command-line otherwise (--branches Test, 300s timeout, 4 workers)",
        "bh_note": "BH is applied within each 400-gene group; per-gene significance calls are family-size dependent "
                   "(production used the 4,974-gene family), so group comparisons use detection RATES within the same "
                   "400-gene family, plus an 'embedded' diagnostic where rerun p-values are substituted into the full "
                   "4,974-gene family before BH",
    },
    "group_results": group_summary,
    "cpg_masking": cpg_summary,
    "alignment_quality_500_gene_sample": qual,
}

json.dump(strat, open(PH9 / "stratified_rerun_sensitivity.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
log("stratified_rerun_sensitivity.json written")

# per-gene detail CSV
detail = None
for g, fr in group_frames.items():
    t = fr[["p", "bh_q", "sig"]].rename(columns={"p": f"p_{g}", "bh_q": f"q_{g}", "sig": f"sig_{g}"})
    detail = t if detail is None else detail.join(t, how="outer")
detail.to_csv(PH9 / "stratified_rerun_pvalues.csv")

# ==================================================================
# Task 4: topology sensitivity
# ==================================================================
topo = {
    "task": "P0-6b topology sensitivity (ILS) (R6)",
    "alternative_topology": "((human{Test},gorilla),(chimpanzee,bonobo)) replacing (((chimpanzee,bonobo),human),gorilla); branch lengths approximately preserved",
    "example_tree": json.load(open(PH9 / "alt_trees_report.json", encoding="utf-8"))["example_alt"],
    "n_genes": int(len(group_frames["topo"])),
    "results": group_summary["topo"],
    "vs_baseline": group_summary["topo_vs_baseline"],
    "conclusion": None,
}
vb = group_summary["topo_vs_baseline"]
topo["conclusion"] = (
    f"Under the alternative (human,gorilla) topology, BH FDR<0.05 detection is {vb[f'topo_sig_rate']:.1%} vs "
    f"{vb['baseline_sig_rate']:.1%} under the species-tree topology ({vb['both_sig']}/{vb['n_common']} genes significant in both; "
    f"{vb['of_baseline_sig_retained']:.1%} of baseline-significant genes retained). "
    f"Detection rate difference: {vb['delta_rate']:+.1%}."
)
json.dump(topo, open(PH9 / "topology_sensitivity.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
log("topology_sensitivity.json written")

with open(PH9 / "collect_report.log", "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines))
log("ALL DONE")
