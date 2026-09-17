# -*- coding: utf-8 -*-
"""
Tier 3 statistical hardening: bootstrap CI for the crossing point n*, the
equal-coverage point, and endpoint extrapolation checks
================================================================================
Inputs : results/phase9_hardening/tier3_campra_sensitivity.{csv,json}
         (k-sweep, 100 replicates/level, saved as mean/SD only -> we re-run the
          identical sweep here and KEEP per-replicate counts for the bootstrap)
Method : replicate-level bootstrap — for each of the 100 replicate series
         (RD counts across the 11 coverage levels), fit OLS RD(n_pos), and
         solve n* = (GD_ref - intercept)/slope; the 2.5/97.5 percentiles of the
         100 n* values give the CI (nonparametric bootstrap over subsampling
         randomness; same classifier logic as phase9_tier3_campra_sensitivity.py).
Also   : equal-coverage point n_eq solving RD(n) = GD(n) with both lines fitted
         per replicate; endpoint validation (extrapolated vs measured ceilings).

Output : results/phase9_hardening/tier3_bootstrap.json ; log _audit_t3.txt
"""
import io
import json

import numpy as np
import pandas as pd

BASE = "d:/人类正选择基因项目/"
LOG = io.open(BASE + "_audit_t3.txt", "w", encoding="utf-8")


def P(*a):
    print(*a, file=LOG)


df = pd.read_csv(BASE + "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv")
n = len(df)
gds = df["gds_v7"].values
sig = df["bh_fdr_sig"].values.astype(bool)
relaxed = df["relax_relaxed"].values.astype(bool)
doan = df["has_doan_campra"].fillna(False).values.astype(bool)
rds_tau, rds_nc, rds_brain = df["rds_tau"].values, df["rds_nc"].values, df["rds_brain"].values
n_hars = df["n_hars"].values

THRESH = 0.15
W_DOAN, W_TAU, W_NC, W_BRAIN = 0.30, 0.25, 0.25, 0.20


def classify_counts(doan_mask):
    rds = W_DOAN * doan_mask + W_TAU * rds_tau + W_NC * rds_nc + W_BRAIN * rds_brain
    gd = (gds > rds + THRESH) & sig
    rd = (rds > gds + THRESH) & (~sig)
    dual = (rds > gds + THRESH) & sig & (gds > 0.3)
    cls = np.where(gd, "gene-driven",
                   np.where(rd, "regulation-driven",
                            np.where(dual, "dual-driven", "neutral")))
    cls = np.where(relaxed & (cls == "gene-driven"), "gene-driven (relaxed)", cls)
    vals, cnts = np.unique(cls, return_counts=True)
    return dict(zip(vals, cnts))


base = classify_counts(doan)
gd_ref = base.get("gene-driven", 0)
P(f"replicated base counts: {base}; GD reference = {gd_ref}")

pos_idx = np.where(doan)[0]
n_pos_total = len(pos_idx)
ks = list(range(0, 101, 10))
rng = np.random.default_rng(20260915)  # same seed as the original script

n_star_reps, n_eq_reps = [], []
RD_reps, GD_reps = [[] for _ in range(100)], [[] for _ in range(100)]
x = np.array([int(round(k / 100 * n_pos_total)) for k in ks], dtype=float)
for k in ks:  # k-major, identical draw order to the original published sweep
    n_keep = int(round(k / 100 * n_pos_total))
    for rep in range(100):
        if n_keep >= n_pos_total:
            keep = pos_idx
        elif n_keep == 0:
            keep = np.array([], dtype=int)
        else:
            keep = rng.choice(pos_idx, size=n_keep, replace=False)
        d = np.zeros(n, dtype=bool)
        d[keep] = True
        c = classify_counts(d)
        RD_reps[rep].append(c.get("regulation-driven", 0))
        GD_reps[rep].append(c.get("gene-driven", 0) + c.get("gene-driven (relaxed)", 0))
for rep in range(100):
    rd_series = np.array(RD_reps[rep], float)
    gd_series = np.array(GD_reps[rep], float)
    s_rd, b_rd = np.polyfit(x, rd_series, 1)
    n_star_reps.append((gd_ref - b_rd) / s_rd)
    s_gd, b_gd = np.polyfit(x, gd_series, 1)
    if (s_rd - s_gd) != 0:
        n_eq_reps.append((b_gd - b_rd) / (s_rd - s_gd))

RD_reps, GD_reps = np.array(RD_reps), np.array(GD_reps)
n_star_reps, n_eq_reps = np.array(n_star_reps), np.array(n_eq_reps)
mean_curve = RD_reps.mean(axis=0)
P("\nmean RD curve matches published CSV:",
  np.allclose(mean_curve, [261.0, 264.38, 267.7, 270.62, 274.14, 276.9, 279.8, 283.46, 286.7, 289.52, 293.0], atol=0.02))

s_point, b_point = np.polyfit(x, mean_curve, 1)
n_star_point = (gd_ref - b_point) / s_point
P(f"\npoint fit (mean curve): slope={s_point:.5f}, intercept={b_point:.2f}, r={np.corrcoef(x, mean_curve)[0,1]:.5f}")
P(f"n* point = {n_star_point:.1f} ({100*n_star_point/n:.2f}% of universe)")

def ci(a):
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))


lo, hi = ci(n_star_reps)
P(f"n* bootstrap (100 replicate fits): median={np.median(n_star_reps):.0f}, "
  f"mean={n_star_reps.mean():.0f}, 95% CI=({lo:.0f}, {hi:.0f}), "
  f"pct of universe CI=({100*lo/n:.1f}%, {100*hi/n:.1f}%)")

# equal-coverage point: RD(n) = GD(n)
mean_gd_curve = GD_reps.mean(axis=0)
sg, bg = np.polyfit(x, mean_gd_curve, 1)
n_eq_point = (bg - b_point) / (s_point - sg)
lo_e, hi_e = ci(n_eq_reps)
P(f"\nGD(n) fit: slope={sg:.5f} (per +1 positive, GD count changes {sg:+.4f}), intercept={bg:.2f}")
P(f"equal-coverage point RD(n)=GD(n): n_eq = {n_eq_point:.0f} ({100*n_eq_point/n:.1f}% of universe), "
  f"bootstrap CI=({lo_e:.0f}, {hi_e:.0f})")

# endpoint validation: extrapolated RD at n_pos = 476 (HAR ceiling) and 4,974 (all)
ceil_all = classify_counts(np.ones(n, dtype=bool))
har_mask = n_hars > 0
ceil_har = classify_counts(har_mask)
rd_476 = b_point + s_point * 476
rd_4974 = b_point + s_point * n
meas_476 = ceil_har.get("regulation-driven", 0)
meas_4974 = ceil_all.get("regulation-driven", 0)
P(f"\nendpoint check: n_pos=476 -> extrapolated RD={rd_476:.0f} vs measured {meas_476} "
  f"(deviation {100*(rd_476-meas_476)/meas_476:+.1f}%)")
P(f"endpoint check: n_pos={n} -> extrapolated RD={rd_4974:.0f} vs measured {meas_4974} "
  f"(deviation {100*(rd_4974-meas_4974)/meas_4974:+.1f}%)")

# also GD at endpoints for completeness
gd_476 = bg + sg * 476
gd_4974 = bg + sg * n
P(f"GD extrapolated at 476: {gd_476:.0f} vs measured {ceil_har.get('gene-driven',0)+ceil_har.get('gene-driven (relaxed)',0)}; "
  f"at {n}: {gd_4974:.0f} vs measured {ceil_all.get('gene-driven',0)+ceil_all.get('gene-driven (relaxed)',0)}")

result = {
    "task": "Tier 3 hardening: bootstrap CI of crossing point, equal-coverage point, endpoint checks",
    "method_note": "replicate-level bootstrap: 100 independent sweep series (same classifier, same "
                   "downsampling design as phase9_tier3_campra_sensitivity.py, seed 20260915); per-series "
                   "OLS fit of RD(n_pos) over the 11 observable levels; percentile CI of the solved roots",
    "gd_reference": int(gd_ref),
    "rd_fit": {"slope": float(s_point), "intercept": float(b_point),
               "r": float(np.corrcoef(x, mean_curve)[0, 1])},
    "crossing_n_star": {
        "point": float(n_star_point),
        "pct_of_universe": float(100 * n_star_point / n),
        "bootstrap_median": float(np.median(n_star_reps)),
        "ci95": [float(lo), float(hi)],
        "ci95_pct_of_universe": [float(100 * lo / n), float(100 * hi / n)],
    },
    "equal_coverage_point": {
        "definition": "n where RD(n) = GD(n) with both counts linear in n_pos (equal-class-size coverage)",
        "gd_fit": {"slope": float(sg), "intercept": float(bg)},
        "point": float(n_eq_point),
        "pct_of_universe": float(100 * n_eq_point / n),
        "bootstrap_ci95": [float(lo_e), float(hi_e)],
    },
    "endpoint_validation": [
        {"scenario": "caMPRA-positive for all 476 HAR-proximate genes",
         "n_pos": 476, "extrapolated_RD": round(float(rd_476)), "measured_RD": int(meas_476),
         "deviation_pct": round(float(100 * (rd_476 - meas_476) / meas_476), 1)},
        {"scenario": "caMPRA-positive for all 4,974 genes",
         "n_pos": int(n), "extrapolated_RD": round(float(rd_4974)), "measured_RD": int(meas_4974),
         "deviation_pct": round(float(100 * (rd_4974 - meas_4974) / meas_4974), 1)},
    ],
    "mean_curve_matches_published": bool(np.allclose(
        mean_curve, pd.read_csv(BASE + "results/phase9_hardening/tier3_campra_sensitivity.csv")["RD_mean"], atol=0.02)),
}
with io.open(BASE + "results/phase9_hardening/tier3_bootstrap.json", "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False)
P("\nsaved results/phase9_hardening/tier3_bootstrap.json")
LOG.close()
print("done")
