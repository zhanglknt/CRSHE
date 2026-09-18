# -*- coding: utf-8 -*-
"""
Phase 9 hardening (R6 fix, P0): hCONDEL hg18->hg38 liftOver repair
===================================================================
Problem (r6_provenance_audit.md §1.3): phase8_hcondels.py Method 2 matched
McLean 2011 hg18 coordinates directly against the hg38 (GRCh38/GENCODE v47)
gene-body BED, without lift-over. Gold-standard check: only 7/44 within-gene
hCONDELs were recovered; 32/44 landed at wrong loci.

Fix (plan A):
  1. liftOver all 583 hCONDEL hg18 coordinates to hg38 (pyliftover;
     fallback panTro2->hg38 where hg18 fails, if chain available).
  2. Re-run Method 2 (coordinate overlap, +/-50 kb, first-hit convention
     identical to the original script) on lifted coordinates.
  3. Recompute every downstream number (full old->new comparison table):
     - hCONDEL gene counts (name / coord / union), coverage
     - RD & GD Fisher enrichment (Table-3 convention: RD only vs rest)
     - HAR x hCONDEL overlap (29 genes, phi, chi-square)
     - LOO variant matrix (full / LOO-caMPRA / LOO-tau / LOO-nc / LOO-brain)
     - Tier 4 subset table (subset-vs-complement, true Newcombe hybrid-score
       CI, same columns as fix-subsets' tier4_empirical_disjointness.csv)
  4. Validation: within-gene recovery rate before vs after liftOver;
     name-method vs fixed-coord-method agreement (Jaccard).

Outputs (results/phase9_hardening/):
  hcondels_liftover_fix.json              -- key numbers as a python dict
  hcondels_liftover_fix_comparison.csv    -- old -> new comparison table
  hcondel_gene_mapping_fixed.csv          -- replacement gene mapping (same
                                             schema as phase8 hcondel_gene_mapping.csv)
  tier4_empirical_disjointness.csv        -- regenerated (old file backed up
                                             as .bak_pre_hcondel_fix)
  log_hcondels_liftover_fix.txt           -- run log

Run: C:/Users/admin/miniconda3/envs/human_selection/python.exe scripts/phase9_hcondels_liftover_fix.py
"""
import json
import os
import re
import shutil
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio

BASE = "d:/人类正选择基因项目"
OUT_DIR = BASE + "/results/phase9_hardening"
LOG_PATH = OUT_DIR + "/log_hcondels_liftover_fix.txt"

_lines = []


def log(msg=""):
    print(msg)
    _lines.append(str(msg))


log("=" * 72)
log("Phase 9 hardening: hCONDEL hg18->hg38 liftOver repair")
log("=" * 72)

# ------------------------------------------------------------------ constants
EXTENSION = 50000
HCONDELS_XLS = BASE + "/data/downloads/hcondels/hCONDELs_supplementary_table2.xls"
GENE_BED = BASE + "/results/phase4_conservation/gene_body.bed"
V7_CSV = BASE + "/results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv"

# old (published) values for the comparison table ------------------------------
OLD = {
    "union_genes": 183, "name_genes": 80, "coord_genes": 121,
    "rd_or": 2.9442596876807405, "rd_p": 6.654409160391586e-06,
    "gd_or": 0.5255006156438339, "gd_p": 0.9995933112970327,
    "class_counts": {"gd": 27, "rd": 27, "neutral": 125},
    "har_overlap": 29, "phi": 0.0417, "chi2_p": 0.0049,
    "coverage_pct": 100 * 183 / 4974,
    "loo_hcondel": {
        "full":      dict(OR=2.9442596876807405, p=6.654409160391586e-06, k_over_n="27/293"),
        "LOO-caMPRA": dict(OR=2.269779911465802, p=2.473308121629363e-06, k_over_n="54/800"),
        "LOO-tau":   dict(OR=3.9113286713286715, p=6.2016958532190185e-06, k_over_n="18/148"),
        "LOO-nc":    dict(OR=3.360544217687075, p=0.000642588903502502, k_over_n="12/110"),
        "LOO-brain": dict(OR=2.2744736842105264, p=1.4227057053097574e-05, k_over_n="43/613"),
    },
    "tier4_hcondel": dict(n=183, busted_sig=56, busted_rate=0.306,
                          complement_n=4791, complement_sig=1634),
    "tier4_combined": dict(n=630, busted_sig=220),
    "within_recovery_direct": 7, "within_universe": 44, "within_wrong_locus": 32,
}

# ------------------------------------------------------------------ load universe
v7 = pd.read_csv(V7_CSV, low_memory=False)
n = len(v7)
all_genes = set(v7["gene_id"])
sym2id = {}
for _, r in v7.iterrows():
    s = str(r["gene_symbol"]).strip()
    if s and s != "nan":
        sym2id[s] = r["gene_id"]
cls = v7["classification_v7"].astype(str).values
sig = v7["bh_fdr_sig"].astype(str).str.strip().str.lower().eq("true").values
relaxed = v7["relax_relaxed"].astype(str).str.strip().str.lower().eq("true").values
rd_set = set(v7.loc[v7["classification_v7"] == "regulation-driven", "gene_id"])
gd_set = set(v7.loc[v7["classification_v7"].isin(["gene-driven", "gene-driven (relaxed)"]), "gene_id"])
dual_set = set(v7.loc[v7["classification_v7"] == "dual-driven", "gene_id"])
har_genes = set(v7.loc[pd.to_numeric(v7["n_hars"], errors="coerce").fillna(0) > 0, "gene_id"])
log(f"Universe: {n} genes | GD(+relaxed)={len(gd_set)}, RD={len(rd_set)}, dual={len(dual_set)}, "
    f"HAR-proximate={len(har_genes)}")

# gene bodies (hg38)
from collections import defaultdict
gene_by_chrom = defaultdict(list)
with open(GENE_BED, encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            gene_by_chrom[p[0]].append((int(p[1]), int(p[2]), p[3]))

# ------------------------------------------------------------------ parse xls
df = pd.read_excel(HCONDELS_XLS)
df = df.iloc[1:].reset_index(drop=True)  # drop sub-header
col_map = {
    df.columns[0]: "name", df.columns[1]: "type", df.columns[2]: "pt2_size",
    df.columns[3]: "pt2_coords", df.columns[4]: "hg18_size", df.columns[5]: "hg18_coords",
    df.columns[6]: "cons", df.columns[7]: "upstream", df.columns[8]: "upstream_dist",
    df.columns[9]: "within", df.columns[10]: "downstream", df.columns[11]: "downstream_dist",
}
df = df.rename(columns=col_map)
log(f"\nhCONDELs parsed from Supplementary Table 2: {len(df)}")


def parse_coord(s):
    m = re.match(r"(chr[\w]+):(\d+)-(\d+)", str(s))
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


# ------------------------------------------------------------------ liftOver
from pyliftover import LiftOver

lo_hg18 = LiftOver("hg18", "hg38")
try:
    lo_pt2 = LiftOver("panTro2", "hg38")
    pt2_chain_available = True
except Exception as e:
    lo_pt2 = None
    pt2_chain_available = False
    log(f"panTro2->hg38 chain unavailable ({e}); fallback disabled")

lift_ok = lift_fail = 0
lift_fail_names = []
pt2_fallback_used = 0
hc_records = []
for _, row in df.iterrows():
    name = str(row["name"])
    hg18 = parse_coord(row["hg18_coords"])
    chrom = s = e = None
    src = None
    if hg18:
        c0, s0, e0 = hg18
        e0 = max(e0, s0)
        r_s = lo_hg18.convert_coordinate(c0, s0)
        r_e = lo_hg18.convert_coordinate(c0, e0) if e0 > s0 else r_s
        if r_s:
            chrom, s = r_s[0][0], r_s[0][1]
            e = r_e[0][1] if (r_e and e0 > s0) else s + 1
            e = max(e, s + 1)
            src = "hg18"
    if chrom is None and pt2_chain_available:
        pt2 = parse_coord(row["pt2_coords"])
        if pt2:
            c0, s0, e0 = pt2
            e0 = max(e0, s0)
            r_s = lo_pt2.convert_coordinate(c0, s0)
            r_e = lo_pt2.convert_coordinate(c0, e0) if e0 > s0 else r_s
            if r_s:
                chrom, s = r_s[0][0], r_s[0][1]
                e = r_e[0][1] if (r_e and e0 > s0) else s + 1
                e = max(e, s + 1)
                src = "panTro2"
                pt2_fallback_used += 1
    if chrom is None:
        lift_fail += 1
        lift_fail_names.append(name)
        continue
    lift_ok += 1
    hc_records.append(dict(name=name, chrom=chrom, start=s, end=e, src=src,
                           within=str(row["within"]).strip(),
                           upstream=str(row["upstream"]).strip(),
                           upstream_dist=row["upstream_dist"],
                           downstream=str(row["downstream"]).strip(),
                           downstream_dist=row["downstream_dist"]))

log(f"liftOver hg18->hg38: success {lift_ok}/{len(df)} "
    f"(panTro2 fallback used: {pt2_fallback_used}); failures: {lift_fail}")
if lift_fail_names:
    log("  failed hCONDELs: " + ", ".join(lift_fail_names))

# ------------------------------------------------------------------ Method 1 (name)
name_genes = set()
for _, row in df.iterrows():
    w = str(row["within"]).strip()
    if w not in ("nan", "") and w in sym2id:
        name_genes.add(sym2id[w])
    for gcol, dcol in (("upstream", "upstream_dist"), ("downstream", "downstream_dist")):
        g = str(row[gcol]).strip()
        if g in ("nan", "") or g not in sym2id:
            continue
        try:
            d = int(float(row[dcol]))
        except Exception:
            continue
        if d <= EXTENSION:
            name_genes.add(sym2id[g])
name_genes &= all_genes

# ------------------------------------------------------------------ Method 2 (fixed coord)
coord_genes = set()
for r in hc_records:
    chrom_genes = gene_by_chrom.get(r["chrom"], [])
    for g_start, g_end, gene_id in chrom_genes:
        if r["start"] < g_end + EXTENSION and r["end"] > max(0, g_start - EXTENSION):
            coord_genes.add(gene_id)
            break  # same first-hit convention as original phase8_hcondels.py
coord_genes &= all_genes

union_genes = (name_genes | coord_genes) & all_genes
log(f"\nGene mapping (fixed): name-only={len(name_genes)}, coord-only-fixed={len(coord_genes)}, "
    f"overlap={len(name_genes & coord_genes)}, union={len(union_genes)} "
    f"(old: 80 / 121 / 18 / 183)")
log(f"  coverage: {len(union_genes)}/{n} = {100*len(union_genes)/n:.2f}% "
    f"(old 183/4974 = {OLD['coverage_pct']:.2f}%)")

# ------------------------------------------------------------------ validation: within-gene recovery
within_universe = within_direct = within_ext = 0
for r in hc_records:
    w = r["within"]
    if w in ("nan", "") or w not in sym2id:
        continue
    within_universe += 1
    eid = sym2id[w]
    for g_start, g_end, gene_id in gene_by_chrom.get(r["chrom"], []):
        if gene_id == eid:
            if r["start"] < g_end and r["end"] > g_start:
                within_direct += 1
            elif r["start"] < g_end + EXTENSION and r["end"] > max(0, g_start - EXTENSION):
                within_ext += 1
            break
log(f"\nGold-standard within-gene recovery (n={within_universe} in universe): "
    f"direct body overlap={within_direct} (old 7), "
    f"+/-50kb only={within_ext} (old 5), wrong locus={within_universe-within_direct-within_ext} (old 32)")

# old coord set for recovery comparison
old_map = pd.read_csv(BASE + "/results/phase8_tissue_analysis/hcondel_gene_mapping.csv")
old_union = set(old_map["gene_id"].dropna().astype(str))
old_coord_only = {g for g, m in zip(old_map["gene_id"], old_map["mapping_method"])
                  if str(m) == "coord"}
kept_old = len(old_coord_only & coord_genes)
log(f"Old coord-only genes ({len(old_coord_only)}) retained by fixed coord method: {kept_old}")

# ------------------------------------------------------------------ enrichment (Table-3 convention: RD vs rest)
def fisher_greater(group_set, feature_set):
    a = len(group_set & feature_set); b = len(group_set - feature_set)
    c = len(feature_set - group_set); d = n - a - b - c
    table = [[a, b], [c, d]]
    orr, p = stats.fisher_exact(table, alternative="greater")
    return orr, p, a


def fisher_two(group_set, feature_set):
    a = len(group_set & feature_set); b = len(group_set - feature_set)
    c = len(feature_set - group_set); d = n - a - b - c
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    return orr, p, a


enrich = {}
for label, fset in (("union_fixed", union_genes), ("name_only", name_genes),
                    ("coord_fixed", coord_genes)):
    orr, p, k = fisher_greater(rd_set, fset)
    orr_g, p_g, k_g = fisher_two(gd_set, fset)
    enrich[label] = dict(n_genes=len(fset), rd_k=k, rd_or=orr, rd_p=p,
                         gd_k=k_g, gd_or=orr_g, gd_p=p_g)
    log(f"  RD enrichment [{label:12s}]: k={k}/{len(rd_set)}, OR={orr:.3f}, p={p:.3g}"
        f" | GD: k={k_g}, OR={orr_g:.3f}, p={p_g:.3g}")
log(f"  (old union: RD OR={OLD['rd_or']:.3f}, p={OLD['rd_p']:.3g}; "
    f"GD OR={OLD['gd_or']:.3f}, p={OLD['gd_p']:.3g})")

# name vs coord agreement
jac = len(name_genes & coord_genes) / len(name_genes | coord_genes)
log(f"  name vs fixed-coord Jaccard: {jac:.3f}")

# class counts of union genes
cc = {c: int(sum(1 for g in union_genes if g in s))
      for c, s in (("gd", gd_set), ("rd", rd_set), ("neutral",
          all_genes - gd_set - rd_set - dual_set), ("dual", dual_set))}
log(f"  class counts (union fixed): {cc} (old gd=27, rd=27, neutral=125)")

# ------------------------------------------------------------------ HAR x hCONDEL overlap
har_arr = v7["gene_id"].isin(har_genes).values
hc_arr = v7["gene_id"].isin(union_genes).values
both = int((har_arr & hc_arr).sum())
phi = float(np.corrcoef(har_arr.astype(float), hc_arr.astype(float))[0, 1])
chi2_p = float(stats.chi2_contingency(np.histogram2d(har_arr.astype(int),
                                                     hc_arr.astype(int), bins=2)[0])[1])
log(f"\nHAR x hCONDEL overlap: {both} genes (old 29); "
    f"{both}/{len(union_genes)} = {100*both/max(len(union_genes),1):.1f}% of hCONDEL genes "
    f"(old 15.8%); {both}/{len(har_genes)} = {100*both/len(har_genes):.1f}% of HAR genes (old 6.1%)")
log(f"  phi={phi:.4f} (old 0.0417), chi2 p={chi2_p:.4g} (old 0.0049)")

# ------------------------------------------------------------------ LOO matrix reproduction
X_gds = v7[["gds_p_pct_resid", "gds_lrt_pct_resid", "gds_relax_pct", "gds_selectome"]].apply(
    pd.to_numeric, errors="coerce").fillna(0).values
rds_doan = pd.to_numeric(v7["rds_doan"], errors="coerce").fillna(0).values
rds_tau = pd.to_numeric(v7["rds_tau"], errors="coerce").fillna(0).values
rds_nc = pd.to_numeric(v7["rds_nc"], errors="coerce").fillna(0).values
rds_brain = pd.to_numeric(v7["rds_brain"], errors="coerce").fillna(0).values
gds = X_gds @ np.array([0.35, 0.30, 0.20, 0.15])
baseline_cls = cls
rds_full = 0.30 * rds_doan + 0.25 * rds_tau + 0.25 * rds_nc + 0.20 * rds_brain


def classify(gds_v, rds_v, threshold=0.15):
    out = np.empty(len(gds_v), dtype=object)
    for i in range(len(gds_v)):
        g, r, s = gds_v[i], rds_v[i], sig[i]
        if g > r + threshold and s:
            c = "gene-driven"
        elif r > g + threshold and not s:
            c = "regulation-driven"
        elif g > r + threshold and r > 0.4 and s:
            c = "dual-driven"
        elif r > g + threshold and s and g > 0.3:
            c = "dual-driven"
        else:
            c = "neutral"
        if relaxed[i] and c == "gene-driven":
            c = "gene-driven (relaxed)"
        out[i] = c
    return out


cls_base = classify(gds, rds_full)
assert np.mean(cls_base == baseline_cls) == 1.0, "baseline reproduction failed"
log("\nBaseline v7 classification reproduced (agreement 100%).")

hc_feat = v7["gene_id"].isin(union_genes).values
hc_feat_old = v7["gene_id"].isin(old_union).values
LOO = {
    "full":       dict(w=None),
    "LOO-caMPRA": dict(w=dict(doan=None, tau=0.35, nc=0.35, brain=0.30)),
    "LOO-tau":    dict(w=dict(doan=0.40, tau=None, nc=0.33, brain=0.27)),
    "LOO-nc":     dict(w=dict(doan=0.40, tau=0.33, nc=None, brain=0.27)),
    "LOO-brain":  dict(w=dict(doan=0.375, tau=0.3125, nc=0.3125, brain=None)),
}
loo_rows = []
for name_, spec in LOO.items():
    if spec["w"] is None:
        c = cls_base
    else:
        w = spec["w"]
        rds_v = (0 if w["doan"] is None else w["doan"] * rds_doan) + \
                (0 if w["tau"] is None else w["tau"] * rds_tau) + \
                (0 if w["nc"] is None else w["nc"] * rds_nc) + \
                (0 if w["brain"] is None else w["brain"] * rds_brain)
        c = classify(gds, rds_v)
    rd_idx = np.where(c == "regulation-driven")[0]
    a = int(hc_feat[rd_idx].sum()); b = int(len(rd_idx) - a)
    bg_a = int(hc_feat.sum()) - a; bg_b = n - len(rd_idx) - bg_a
    orr, p = stats.fisher_exact([[a, b], [bg_a, bg_b]], alternative="greater")
    # HAR for sanity (should reproduce old values)
    a_h = int(har_arr[rd_idx].sum()); b_h = int(len(rd_idx) - a_h)
    bg_ah = int(har_arr.sum()) - a_h; bg_bh = n - len(rd_idx) - bg_ah
    orr_h, p_h = stats.fisher_exact([[a_h, b_h], [bg_ah, bg_bh]], alternative="greater")
    loo_rows.append(dict(variant=name_, RD_n=len(rd_idx), hc_k=a,
                         hc_or=orr, hc_p=p, har_or=orr_h, har_p=p_h))
    log(f"  {name_:11s}: RD={len(rd_idx):4d}, hCONDEL k={a}/{len(rd_idx)}, "
        f"OR={orr:.3f}, p={p:.3g} (old OR={OLD['loo_hcondel'][name_]['OR']:.3f}, "
        f"p={OLD['loo_hcondel'][name_]['p']:.3g}, k={OLD['loo_hcondel'][name_]['k_over_n']}) | "
        f"HAR OR={orr_h:.3f} (unchanged)")
new_loo_min = min(r["hc_or"] for r in loo_rows)
new_loo_max = max(r["hc_or"] for r in loo_rows)
new_p_min = min(r["hc_p"] for r in loo_rows)
new_p_max = max(r["hc_p"] for r in loo_rows)
log(f"  new LOO hCONDEL OR range: {new_loo_min:.2f}-{new_loo_max:.2f} (old 2.27-3.91); "
    f"p range: {new_p_min:.3g}-{new_p_max:.3g} (old 2.47e-06-6.43e-04)")

# ------------------------------------------------------------------ Tier 4 subset table
busted = v7["bh_fdr_sig"].astype(str).str.strip().str.lower().eq("true").values
rd_or_dual = v7["classification_v7"].isin(["regulation-driven", "dual-driven"]).values
genome_sig = int(busted.sum())


def wilson(k, nn, z=1.959963984540054):
    p = k / nn
    d = 1 + z * z / nn
    c = (p + z * z / (2 * nn)) / d
    h = z * np.sqrt(p * (1 - p) / nn + z * z / (4 * nn * nn)) / d
    return c - h, c + h


def newcombe_diff(k1, n1, k2, n2):
    p1, p2 = k1 / n1, k2 / n2
    l1, u1 = wilson(k1, n1)
    l2, u2 = wilson(k2, n2)
    low = (p1 - p2) - np.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    high = (p1 - p2) + np.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    wbs_low, wbs_high = l1 - u2, u1 - l2
    se = np.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    wald_low, wald_high = (p1 - p2) - 1.959963984540054 * se, (p1 - p2) + 1.959963984540054 * se
    return low, high, wbs_low, wbs_high, wald_low, wald_high


# sanity: reproduce old hCONDEL row with old set
_t = newcombe_diff(56, 183, 1634, 4791)
assert abs(_t[0] - (-0.0987)) < 1e-3 and abs(_t[1] - 0.0364) < 1e-3, \
    f"Newcombe reproduction failed: {_t}"
log("\nNewcombe method-10 reproduction of old hCONDEL row verified.")

subsets = {
    "HAR-proximate (n_hars>0)": har_arr,
    "hCONDEL-proximate": hc_arr,
    "HAR or hCONDEL": har_arr | hc_arr,
}
tier4_rows = []
for sname, mask in subsets.items():
    k = int(mask.sum())
    s_sig = int(busted[mask].sum())
    comp_n = n - k
    comp_sig = genome_sig - s_sig
    table = [[s_sig, k - s_sig], [comp_sig, comp_n - comp_sig]]
    orr, p = stats.fisher_exact(table)
    rd = int(rd_or_dual[mask].sum())
    rd_all = int(rd_or_dual.sum())
    tab_rd = [[rd, k - rd], [rd_all - rd, (n - k) - (rd_all - rd)]]
    orr_rd, p_rd = stats.fisher_exact(tab_rd)
    gd = int(v7["classification_v7"].isin(["gene-driven", "gene-driven (relaxed)"]).values[mask].sum())
    low, high, wbs_l, wbs_h, wald_l, wald_h = newcombe_diff(s_sig, k, comp_sig, comp_n)
    tier4_rows.append({
        "subset": sname, "n": k, "busted_sig": s_sig,
        "busted_rate": round(s_sig / k, 4), "busted_OR": round(orr, 3), "busted_p": p,
        "RD_or_dual": rd, "RD_rate": round(rd / k, 4), "RD_OR": round(orr_rd, 3), "RD_p": p_rd,
        "GD_class": gd, "GD_rate": round(gd / k, 4),
        "complement_n": comp_n, "complement_sig": comp_sig,
        "complement_rate": round(comp_sig / comp_n, 4),
        "rate_diff": round(s_sig / k - comp_sig / comp_n, 4),
        "rate_diff_ci95_newcombe_low": round(low, 4),
        "rate_diff_ci95_newcombe_high": round(high, 4),
        "rate_diff_ci95_wilson_bound_subtraction_low": round(wbs_l, 4),
        "rate_diff_ci95_wilson_bound_subtraction_high": round(wbs_h, 4),
        "rate_diff_ci95_wald_low": round(wald_l, 4),
        "rate_diff_ci95_wald_high": round(wald_h, 4),
        "fisher_p_vs_complement": p,
    })
    log(f"  {sname}: n={k}, BUSTED={s_sig} ({s_sig/k:.1%}) vs complement "
        f"{comp_sig}/{comp_n} ({comp_sig/comp_n:.1%}), diff="
        f"{s_sig/k-comp_sig/comp_n:+.4f} [{low:.4f}, {high:.4f}], fisher p={p:.3g}")

tier4 = pd.DataFrame(tier4_rows)
# backup + rewrite
T4 = OUT_DIR + "/tier4_empirical_disjointness.csv"
if not os.path.exists(T4 + ".bak_pre_hcondel_fix"):
    shutil.copy2(T4, T4 + ".bak_pre_hcondel_fix")
    log("  (old tier4 CSV backed up to tier4_empirical_disjointness.csv.bak_pre_hcondel_fix)")
tier4.to_csv(T4, index=False)

# verify HAR row reproduced unchanged
old_t4 = pd.read_csv(T4 + ".bak_pre_hcondel_fix")
_h = tier4.iloc[0]
_oh = old_t4.iloc[0]
assert _h["n"] == _oh["n"] and _h["busted_sig"] == _oh["busted_sig"] and \
    abs(_h["rate_diff"] - _oh["rate_diff"]) < 1e-9, "HAR row reproduction failed"
log("  HAR row reproduced identically (n=476, sig=174, diff=0.0285).")

# ------------------------------------------------------------------ outputs
os.makedirs(OUT_DIR, exist_ok=True)

# fixed mapping CSV (same schema as phase8 hcondel_gene_mapping.csv)
rows_map = []
cls_lookup = dict(zip(v7["gene_id"], v7["classification_v7"]))
sym_lookup = dict(zip(v7["gene_id"], v7["gene_symbol"]))
for g in sorted(union_genes):
    m = []
    if g in name_genes:
        m.append("name")
    if g in coord_genes:
        m.append("coord")
    rows_map.append({"gene_id": g, "gene_symbol": sym_lookup.get(g, ""),
                     "classification_v7": cls_lookup.get(g, ""),
                     "mapping_method": "+".join(m)})
pd.DataFrame(rows_map).to_csv(OUT_DIR + "/hcondel_gene_mapping_fixed.csv", index=False)

# comparison table
comp_rows = [
    ("liftOver success (of 583)", "0 (not attempted)", f"{lift_ok}"),
    ("name-method genes (universe)", "80", f"{len(name_genes)}"),
    ("coord-method genes (universe, hg18-as-hg38 BUG)", "121", f"{len(coord_genes)}"),
    ("name-coord overlap", "18", f"{len(name_genes & coord_genes)}"),
    ("union hCONDEL genes", "183", f"{len(union_genes)}"),
    ("coverage of 4,974 universe", f"{OLD['coverage_pct']:.2f}%",
     f"{100*len(union_genes)/n:.2f}%"),
    ("RD enrichment OR (union)", f"{OLD['rd_or']:.3f}", f"{enrich['union_fixed']['rd_or']:.3f}"),
    ("RD enrichment p (union)", f"{OLD['rd_p']:.3g}", f"{enrich['union_fixed']['rd_p']:.3g}"),
    ("GD enrichment OR (union)", f"{OLD['gd_or']:.3f}", f"{enrich['union_fixed']['gd_or']:.3f}"),
    ("GD enrichment p (union)", f"{OLD['gd_p']:.3g}", f"{enrich['union_fixed']['gd_p']:.3g}"),
    ("RD class count among hCONDEL genes", "27", f"{cc['rd']}"),
    ("GD class count among hCONDEL genes", "27", f"{cc['gd']}"),
    ("neutral count among hCONDEL genes", "125", f"{cc['neutral']}"),
    ("HAR & hCONDEL double-evidence genes", "29", f"{both}"),
    ("% of hCONDEL genes also HAR-proximate", "15.8%", f"{100*both/max(len(union_genes),1):.1f}%"),
    ("% of HAR genes also hCONDEL-proximate", "6.1%", f"{100*both/len(har_genes):.1f}%"),
    ("phi (HAR x hCONDEL)", "0.0417", f"{phi:.4f}"),
    ("chi2 p (HAR x hCONDEL)", "0.0049", f"{chi2_p:.4g}"),
    ("LOO hCONDEL OR range (5 variants)", "2.27-3.91", f"{new_loo_min:.2f}-{new_loo_max:.2f}"),
    ("LOO hCONDEL p range (5 variants)", "2.47e-06-6.43e-04",
     f"{new_p_min:.3g}-{new_p_max:.3g}"),
    ("Tier4 hCONDEL n", "183", f"{int(tier4.iloc[1]['n'])}"),
    ("Tier4 hCONDEL BUSTED rate", "0.306", f"{tier4.iloc[1]['busted_rate']}"),
    ("Tier4 hCONDEL rate_diff vs complement", "-0.035",
     f"{tier4.iloc[1]['rate_diff']}"),
    ("Tier4 hCONDEL Newcombe 95% CI", "[-0.0987, 0.0364]",
     f"[{tier4.iloc[1]['rate_diff_ci95_newcombe_low']}, {tier4.iloc[1]['rate_diff_ci95_newcombe_high']}]"),
    ("Tier4 hCONDEL RD_or_dual", "29", f"{int(tier4.iloc[1]['RD_or_dual'])}"),
    ("Tier4 'HAR or hCONDEL' n", "630", f"{int(tier4.iloc[2]['n'])}"),
    ("Tier4 'HAR or hCONDEL' BUSTED rate", "0.3492", f"{tier4.iloc[2]['busted_rate']}"),
    ("within-gene direct recovery (of 44)", "7", f"{within_direct}"),
    ("within-gene wrong locus (of 44)", "32", f"{within_universe-within_direct-within_ext}"),
    ("old coord-only genes retained", "-", f"{kept_old}/{len(old_coord_only)}"),
    ("name vs coord Jaccard", "-", f"{jac:.3f}"),
]
pd.DataFrame(comp_rows, columns=["metric", "old", "new_fixed"]).to_csv(
    OUT_DIR + "/hcondels_liftover_fix_comparison.csv", index=False)

# JSON summary
result = {
    "task": "R6 P0 fix: hCONDEL hg18->hg38 liftOver repair (provenance audit 1.3, plan A)",
    "universe_n": n,
    "liftover": {"total": int(len(df)), "success": int(lift_ok), "fail": int(lift_fail),
                 "failed_names": lift_fail_names,
                 "pantro2_fallback_used": int(pt2_fallback_used),
                 "pantro2_chain_available": bool(pt2_chain_available)},
    "gene_sets": {
        "name_only": len(name_genes), "coord_fixed": len(coord_genes),
        "name_coord_overlap": len(name_genes & coord_genes),
        "union_fixed": len(union_genes),
        "union_old": OLD["union_genes"],
        "coverage_pct_new": round(100 * len(union_genes) / n, 2),
        "coverage_pct_old": round(OLD["coverage_pct"], 2),
        "jaccard_name_vs_coord": round(jac, 4),
        "old_coord_only_retained": f"{kept_old}/{len(old_coord_only)}",
    },
    "validation_within_gene": {
        "n_in_universe": within_universe, "direct_recovery_new": within_direct,
        "direct_recovery_old": OLD["within_recovery_direct"],
        "ext50kb_only_new": within_ext,
        "wrong_locus_new": within_universe - within_direct - within_ext,
        "wrong_locus_old": OLD["within_wrong_locus"],
    },
    "enrichment": {
        "old_union": {"rd_or": OLD["rd_or"], "rd_p": OLD["rd_p"],
                      "gd_or": OLD["gd_or"], "gd_p": OLD["gd_p"]},
        "new": {k: {kk: (float(vv) if isinstance(vv, (int, float, np.floating)) else vv)
                    for kk, vv in v.items()} for k, v in enrich.items()},
    },
    "class_counts_union": cc,
    "har_hcondel_overlap": {
        "both": both, "pct_of_hcondel": round(100 * both / max(len(union_genes), 1), 1),
        "pct_of_har": round(100 * both / len(har_genes), 1),
        "phi": round(phi, 4), "chi2_p": chi2_p,
        "old": {"both": 29, "phi": 0.0417, "chi2_p": 0.0049},
    },
    "loo_hcondel": {
        r["variant"]: {"OR": round(r["hc_or"], 3), "p": r["hc_p"],
                       "k_over_n": f"{r['hc_k']}/{r['RD_n']}",
                       "old_OR": round(OLD["loo_hcondel"][r["variant"]]["OR"], 3),
                       "old_p": OLD["loo_hcondel"][r["variant"]]["p"],
                       "old_k_over_n": OLD["loo_hcondel"][r["variant"]]["k_over_n"]}
        for r in loo_rows},
    "loo_hcondel_or_range": {"new": [round(new_loo_min, 2), round(new_loo_max, 2)],
                             "old": [2.27, 3.91]},
    "loo_hcondel_p_range": {"new": [new_p_min, new_p_max], "old": [2.47e-06, 6.43e-04]},
    "tier4_rows": tier4_rows,
    "files": {
        "comparison_csv": OUT_DIR + "/hcondels_liftover_fix_comparison.csv",
        "fixed_mapping_csv": OUT_DIR + "/hcondel_gene_mapping_fixed.csv",
        "tier4_csv_updated": T4,
        "tier4_backup": T4 + ".bak_pre_hcondel_fix",
        "log": LOG_PATH,
    },
}
with open(OUT_DIR + "/hcondels_liftover_fix.json", "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False, default=float)

with open(LOG_PATH, "w", encoding="utf-8") as f:
    f.write("\n".join(_lines))

print("\nSaved:", OUT_DIR + "/hcondels_liftover_fix.json")
