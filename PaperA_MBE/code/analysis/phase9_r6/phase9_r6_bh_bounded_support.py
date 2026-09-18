# -*- coding: utf-8 -*-
"""Task 1 (P0-1b): BH bounded-support correction check.

Two mechanisms for a bounded (support-limited) null p-value distribution:
  A. Anti-conservative compression: null ~ U(0, 0.5)  -> corrected p* = p / 0.5 (capped at 1)
  B. Truncation/clipping: p > 0.5 clipped to 0.5      -> BH ordering unchanged, BH result unchanged

Baseline: BH FDR<0.05 on original 4,974 BUSTED p-values = 1,690 significant genes.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(r"d:\人类正选择基因项目")
CSV = BASE / "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv"
OUT_DIR = BASE / "results/phase9_hardening"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "bh_bounded_support_check.json"

log_lines = []


def log(s):
    print(s)
    log_lines.append(str(s))


df = pd.read_csv(CSV)
log(f"n_genes={len(df)}")
p = df["busted_p"].astype(float).values
n = len(p)
n_nan = int(np.isnan(p).sum())
log(f"nan_p={n_nan}")
p = p[~np.isnan(p)]
n_valid = len(p)

# ---- Baseline BH (verify 1,690) ----
def bh(pvals, alpha=0.05):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals, kind="mergesort")
    ranked = pvals[order]
    thresh = np.arange(1, n + 1) / n * alpha
    below = ranked <= thresh
    if not below.any():
        return 0, np.full(n, np.nan), np.zeros(n, dtype=bool)
    k = np.max(np.nonzero(below)[0]) + 1
    q = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.clip(q, 0, 1)
    rej = np.zeros(n, dtype=bool)
    rej[order[:k]] = True
    return int(k), q, rej

k_base, q_base, rej_base = bh(p)
log(f"baseline_BH_sig={k_base}")

# ---- Mechanism A: anti-conservative compression, null U(0,0.5) -> p* = p/0.5 ----
p_star = np.minimum(p / 0.5, 1.0)
k_star, _, _ = bh(p_star)
# equivalently BH at alpha/2 on original p
k_half, _, _ = bh(p, alpha=0.025)
log(f"mechanismA_compression_BH_sig={k_star} (check alpha/2 on raw p: {k_half})")
n_p_lt_05_raw = int((p < 0.05).sum())
n_pstar_lt_05 = int((p_star < 0.05).sum())
log(f"nominal_p<0.05: raw={n_p_lt_05_raw}, p*={n_pstar_lt_05} (p*<0.05 <=> p<0.025)")

# ---- Mechanism B: truncation, p>0.5 clipped to 0.5 ----
p_clip = np.minimum(p, 0.5)
k_clip, _, _ = bh(p_clip)
log(f"mechanismB_truncation_BH_sig={k_clip}")
# analytic argument: clipping to 0.5 cannot create rejections because 0.5 <= i/n*0.05 needs i/n>=10
max_thresh = 0.05  # i/n * alpha <= alpha = 0.05 << 0.5
log(f"max BH threshold (i=n) = {max_thresh} < 0.5, so clipped mass at 0.5 is never rejected; "
    f"and ordering of p<0.5 unchanged -> BH identical")

# point mass at exactly 0.5 in observed data
n_at_05 = int((p == 0.5).sum())
n_gt_05 = int((p > 0.5).sum())
log(f"observed p==0.5: {n_at_05}, p>0.5: {n_gt_05}")

# p-value distribution summary
qs = np.percentile(p, [0, 5, 25, 50, 75, 90, 95, 100])
log("p quantiles [0,5,25,50,75,90,95,100]: " + ", ".join(f"{x:.4g}" for x in qs))

result = {
    "task": "P0-1b bounded-support correction check (R6)",
    "n_genes_total": int(len(df)),
    "n_valid_p": int(n_valid),
    "n_nan_p": n_nan,
    "baseline": {
        "bh_alpha": 0.05,
        "n_significant": int(k_base),
        "note": "BH on original BUSTED p-values (busted_results_v3 pipeline, n=4,974)"
    },
    "mechanism_A_anticonservative_compression": {
        "assumption": "null p ~ U(0, 0.5); observed p deflated 2x relative to Uniform(0,1)",
        "correction": "p* = p / 0.5, capped at 1",
        "bh_alpha": 0.05,
        "n_significant_corrected": int(k_star),
        "equivalent_to": "BH at alpha=0.025 on original p (monotone transform preserves order)",
        "change_vs_baseline": int(k_star) - int(k_base),
        "nominal_p_lt_005_raw": n_p_lt_05_raw,
        "nominal_pstar_lt_005": n_pstar_lt_05,
    },
    "mechanism_B_truncation_clipping": {
        "assumption": "p > 0.5 values are truncated/clipped to the 0.5 boundary (point mass at 0.5)",
        "correction": "p* = min(p, 0.5)",
        "bh_alpha": 0.05,
        "n_significant_corrected": int(k_clip),
        "change_vs_baseline": int(k_clip) - int(k_base),
        "argument": "BH critical values i/n*alpha <= alpha = 0.05 << 0.5, so the clipped point mass at 0.5 "
                    "can never cross a BH threshold; and the relative order of all p<0.5 values is "
                    "unchanged by clipping the upper tail, so the BH rejection set is identical to baseline.",
    },
    "observed_p_distribution": {
        "n_p_eq_0.5": n_at_05,
        "n_p_gt_0.5": n_gt_05,
        "quantiles": {str(q): float(v) for q, v in zip([0, 5, 25, 50, 75, 90, 95, 100], qs)},
    },
    "conclusion": (
        f"Under the truncation mechanism (B), the 1,690 BH-significant genes are entirely unchanged "
        f"({k_clip} = 1,690): clipping cannot create or destroy rejections. Under the anti-conservative "
        f"compression mechanism (A, null U(0,0.5)), correcting p* = 2p reduces BH significance from "
        f"1,690 to {k_star} genes ({100.0*k_star/max(k_base,1):.1f}% of baseline). Which mechanism "
        f"actually operates is determined empirically by the neutral simulation calibration "
        f"(neutral_sim_calibration.json)."
    ),
}

with open(OUT, "w", encoding="utf-8") as f:
    json.dump(result, f, ensure_ascii=False, indent=2)

with open(OUT_DIR / "bh_bounded_support_check.log", "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines))

print("WROTE", OUT)
