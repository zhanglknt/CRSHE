# -*- coding: utf-8 -*-
"""Round-2 computational fixes:
(A) k-means multi-start stability (100 seeds) for the unsupervised control
(B) GC3 / expression covariate sensitivity for CDS-length residualization
Outputs -> results/paper/revision_v7/
  kmeans_multistart.json, covariate_sensitivity.csv, covariate_sensitivity_summary.json
"""
import io, os, json
import numpy as np
import pandas as pd
from scipy.stats import rankdata, zscore
from sklearn.cluster import KMeans

BASE = os.environ.get("HSD_BASE") or "D:/人类正选择基因项目"
OUT = f"{BASE}/results/paper/revision_v7"

v7 = pd.read_csv(f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv")
n = len(v7)

# ============================================================
# A) k-means multi-start stability
# ============================================================
feats = ["gds_p_pct_resid", "gds_lrt_pct_resid", "gds_relax_pct", "gds_selectome",
         "rds_doan", "rds_tau", "rds_nc", "rds_brain"]
X = v7[feats].astype(float).values
Xz = (X - X.mean(0)) / X.std(0)
rd_mask = (v7["classification_v7"] == "regulation-driven").values
n_rd = int(rd_mask.sum())

rec, pur = [], []
for seed in range(100):
    km = KMeans(n_clusters=4, n_init=10, random_state=seed).fit(Xz)
    lab = km.labels_
    # cluster with most RD genes
    counts = np.array([(lab[rd_mask] == c).sum() for c in range(4)])
    best = int(counts.max())
    rec.append(best / n_rd)
    # purity of that cluster: RD share within it
    bc = int(np.argmax(counts))
    pur.append(best / max(int((lab == bc).sum()), 1))
rec = np.array(rec); pur = np.array(pur)
km_res = {
    "n_seeds": 100, "n_init_per_seed": 10, "k": 4, "n_rd": n_rd,
    "rd_recovery_mean": float(rec.mean()), "rd_recovery_sd": float(rec.std()),
    "rd_recovery_min": float(rec.min()), "rd_recovery_max": float(rec.max()),
    "rd_recovery_median": float(np.median(rec)),
    "frac_seeds_recovery_ge_0.80": float((rec >= 0.80).mean()),
    "frac_seeds_recovery_ge_0.85": float((rec >= 0.85).mean()),
    "cluster_purity_median": float(np.median(pur)),
    "note": "recovery = fraction of the 293 regulation-driven genes captured by the best cluster; "
            "features are the same 8 standardized component scores used by the classifier",
}
with open(f"{OUT}/kmeans_multistart.json", "w") as f:
    json.dump(km_res, f, indent=2)
print("k-means multi-start:", json.dumps(km_res, indent=2)[:600])

# ============================================================
# B) GC3 + median-expression covariate sensitivity
# ============================================================
# --- B1: GC3 from human CDS FASTA ---
gc3_map, gc_map = {}, {}
gene_key = v7["gene_id"].astype(str).str.split(".").str[0]
want = set(gene_key)
cur_id, seq = None, []
with open(f"{BASE}/results/phase1_gene_extraction/human_genes_cds.fasta", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith(">"):
            if cur_id is not None and cur_id in want and cur_id not in gc3_map:
                s = "".join(seq)
                cod = s[: (len(s) // 3) * 3]
                third = cod[2::3].upper()
                if len(third) >= 30:
                    gc3_map[cur_id] = (third.count("G") + third.count("C")) / len(third)
                    gc_map[cur_id] = (s.upper().count("G") + s.upper().count("C")) / len(s)
            cur_id = line[1:].strip().split(".")[0].split()[0]
            seq = []
        else:
            seq.append(line.strip())
if cur_id is not None and cur_id in want and cur_id not in gc3_map:
    s = "".join(seq)
    cod = s[: (len(s) // 3) * 3]
    third = cod[2::3].upper()
    if len(third) >= 30:
        gc3_map[cur_id] = (third.count("G") + third.count("C")) / len(third)
        gc_map[cur_id] = (s.upper().count("G") + s.upper().count("C")) / len(s)

v7["gc3"] = gene_key.map(gc3_map)
v7["gc"] = gene_key.map(gc_map)
print(f"GC3 matched: {v7['gc3'].notna().sum()}/{n}")

# --- B2: median GTEx expression covariate ---
tau = pd.read_csv(f"{BASE}/results/gtex_v11_tau/gtex_v11_tau.csv")
tissue_cols = [c for c in tau.columns if c not in ("gene_id", "tau")]
tau["median_tpm"] = tau[tissue_cols].median(axis=1)
tau["gene_key"] = tau["gene_id"].astype(str).str.split(".").str[0]
med = tau.set_index("gene_key")["median_tpm"].to_dict()
v7["med_expr"] = gene_key.map(med)
print(f"median expr matched: {v7['med_expr'].notna().sum()}/{n}")

# --- B3: replicate baseline residualization & classification ---
import numpy.polynomial.polynomial as P

def resid_pct(y, X_):
    """OLS residuals of y on X_ (with intercept), converted to percentile ranks."""
    Xd = np.column_stack([np.ones(len(y))] + [X_[:, j] for j in range(X_.shape[1])])
    beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    r = y - Xd @ beta
    return rankdata(r) / len(r)

def classify(gds, rds, sig, relaxed, thr=0.15):
    cls = np.where((gds > rds + thr) & sig, "gene-driven",
          np.where((rds > gds + thr) & ~sig, "regulation-driven",
          np.where((rds > gds + thr) & sig & (gds > 0.3), "dual-driven", "neutral")))
    cls = np.where(relaxed & (cls == "gene-driven"), "gene-driven (relaxed)", cls)
    return cls

cds = np.log10(v7["cds_length"].astype(float).values)
neglogp = -np.log10(np.maximum(v7["busted_p"].astype(float), 1e-300)).values
lrt = v7["busted_lrt"].astype(float).values
sig = v7["bh_fdr_sig"].astype(bool).values
relaxed = v7["relax_relaxed"].astype(bool).values
rds_full = v7["rds_v7"].astype(float).values
relax_pct = v7["gds_relax_pct"].astype(float).values
sel = v7["gds_selectome"].astype(float).values

# baseline replication (length only)
p_pct0 = resid_pct(neglogp, cds[:, None])
l_pct0 = resid_pct(lrt, cds[:, None])
gds0 = 0.35 * p_pct0 + 0.30 * l_pct0 + 0.20 * relax_pct + 0.15 * sel
cls0 = classify(gds0, rds_full, sig, relaxed)
agree0 = (cls0 == v7["classification_v7"].values).mean()
print(f"baseline replication agreement: {agree0:.4f}")

gc3 = v7["gc3"].fillna(v7["gc3"].median()).values
log_expr = np.log10(np.maximum(v7["med_expr"].fillna(v7["med_expr"].median()).values, 0.01))

variants = {
    "baseline (log10 CDS length only)": cds[:, None],
    "+ GC3": np.column_stack([cds, gc3]),
    "+ log10 median expression": np.column_stack([cds, log_expr]),
    "+ GC3 + log10 median expression": np.column_stack([cds, gc3, log_expr]),
}
rows = []
cls_base = v7["classification_v7"].values
for name, Xc in variants.items():
    p_pct = resid_pct(neglogp, Xc)
    l_pct = resid_pct(lrt, Xc)
    gds = 0.35 * p_pct + 0.30 * l_pct + 0.20 * relax_pct + 0.15 * sel
    cls = classify(gds, rds_full, sig, relaxed)
    ag = (cls == cls_base).mean()
    gd = int((cls == "gene-driven").sum()); gdr = int((cls == "gene-driven (relaxed)").sum())
    rd = int((cls == "regulation-driven").sum()); du = int((cls == "dual-driven").sum())
    ne = int((cls == "neutral").sum())
    rows.append({"variant": name, "agreement_with_baseline": round(ag, 4),
                 "gene_driven": gd, "gd_relaxed": gdr, "regulation_driven": rd,
                 "dual": du, "neutral": ne})
    print(rows[-1])

pd.DataFrame(rows).to_csv(f"{OUT}/covariate_sensitivity.csv", index=False)
with open(f"{OUT}/covariate_sensitivity_summary.json", "w") as f:
    json.dump({"baseline_replication_agreement": float(agree0),
               "gc3_matched": int(v7["gc3"].notna().sum()),
               "median_expr_matched": int(v7["med_expr"].notna().sum()),
               "variants": rows}, f, indent=2)
print("saved covariate sensitivity")
