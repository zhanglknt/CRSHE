# -*- coding: utf-8 -*-
"""
Tier 3: caMPRA coverage-sensitivity curve and extrapolated crossing point k*
=============================================================================
Question (associate researcher critique #4): regulatory evidence covers only
~2% of the universe (97 caMPRA-positive genes) vs ~100% for coding evidence.
How much regulatory-assay coverage would be needed for regulation-driven (RD)
to match gene-driven (GD = 1,214)?

Design:
  1. Replicate the v7 classifier exactly (sanity check against the saved
     classification_v7 column).
  2. Downsample the 97 caMPRA-positive genes: keep k% (k = 0..100, 100 reps),
     zero out rds_doan for the rest, recompute RDS and reclassify. This maps
     RD count as a function of the number of assay-positive genes (n_pos).
  3. Test linearity of RD(n_pos) on the observable range, fit a line, and
     extrapolate to the crossing point where RD = GD = 1,214.
  4. Theoretical ceilings:
     (a) doan=1 for ALL genes (perfect assay, every gene active);
     (b) doan=1 for the 476 HAR-proximate genes only (realistic ceiling:
         caMPRA targets HARs, so only HAR-proximate genes can be positive).

Output:
  results/phase9_hardening/tier3_campra_sensitivity.json
  results/phase9_hardening/tier3_campra_sensitivity.csv
  results/phase9_hardening/figures/Figure_Tier3_coverage_sensitivity.png/.pdf
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = 'd:/人类正选择基因项目/'
OUT = BASE + 'results/phase9_hardening/'
FIG = OUT + 'figures/'
os.makedirs(FIG, exist_ok=True)

THRESH = 0.15
rng = np.random.default_rng(20260915)

# ------------------------------------------------------------------ load
df = pd.read_csv(BASE + 'results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv')
n = len(df)
gds = df['gds_v7'].values
sig = df['bh_fdr_sig'].values.astype(bool)
relaxed = df['relax_relaxed'].values.astype(bool)
doan = df['has_doan_campra'].fillna(False).values.astype(bool)
rds_tau = df['rds_tau'].values
rds_nc = df['rds_nc'].values
rds_brain = df['rds_brain'].values
n_hars = df['n_hars'].values
cls_saved = df['classification_v7'].values

W_DOAN, W_TAU, W_NC, W_BRAIN = 0.30, 0.25, 0.25, 0.20


def classify_counts(doan_mask):
    """Vectorized v7 classification; returns class-count dict."""
    rds = W_DOAN * doan_mask + W_TAU * rds_tau + W_NC * rds_nc + W_BRAIN * rds_brain
    gd = (gds > rds + THRESH) & sig
    rd = (rds > gds + THRESH) & (~sig)
    dual = (rds > gds + THRESH) & sig & (gds > 0.3)
    # precedence identical to original if/elif chain: gd first, then rd, then dual
    cls = np.where(gd, 'gene-driven',
                   np.where(rd, 'regulation-driven',
                            np.where(dual, 'dual-driven', 'neutral')))
    cls = np.where(relaxed & (cls == 'gene-driven'), 'gene-driven (relaxed)', cls)
    vals, cnts = np.unique(cls, return_counts=True)
    return dict(zip(vals, cnts))


# ------------------------------------------------------------------ sanity
base_counts = classify_counts(doan)
print('replicated v7 counts:', base_counts)
gd_ref = base_counts.get('gene-driven', 0)
match = 0
# full per-gene replication check
rds_full = W_DOAN * doan + W_TAU * rds_tau + W_NC * rds_nc + W_BRAIN * rds_brain
gd_m = (gds > rds_full + THRESH) & sig
rd_m = (rds_full > gds + THRESH) & (~sig)
du_m = (rds_full > gds + THRESH) & sig & (gds > 0.3)
cls_rep = np.where(gd_m, 'gene-driven', np.where(rd_m, 'regulation-driven', np.where(du_m, 'dual-driven', 'neutral')))
cls_rep = np.where(relaxed & (cls_rep == 'gene-driven'), 'gene-driven (relaxed)', cls_rep)
match = float((cls_rep == cls_saved).mean())
print(f'per-gene replication agreement: {match:.4f}')
assert match > 0.999, 'classifier replication failed'

pos_idx = np.where(doan)[0]
n_pos_total = len(pos_idx)
print(f'caMPRA positives: {n_pos_total}, GD reference: {gd_ref}')

# ------------------------------------------------------------------ sweep
ks = list(range(0, 101, 10))
rows = []
for k in ks:
    n_keep = int(round(k / 100 * n_pos_total))
    rd_list, gd_list, dual_list = [], [], []
    for rep in range(100):
        if n_keep >= n_pos_total:
            keep = pos_idx
        elif n_keep == 0:
            keep = np.array([], dtype=int)
        else:
            keep = rng.choice(pos_idx, size=n_keep, replace=False)
        d = np.zeros(n, dtype=float)
        d[keep] = 1.0
        c = classify_counts(d.astype(bool))
        rd_list.append(c.get('regulation-driven', 0))
        gd_list.append(c.get('gene-driven', 0) + c.get('gene-driven (relaxed)', 0))
        dual_list.append(c.get('dual-driven', 0))
    rows.append({'k_pct': k, 'n_pos': n_keep,
                 'RD_mean': float(np.mean(rd_list)), 'RD_sd': float(np.std(rd_list)),
                 'GD_mean': float(np.mean(gd_list)), 'GD_sd': float(np.std(gd_list)),
                 'dual_mean': float(np.mean(dual_list))})
    print(f'k={k:3d}% n_pos={n_keep:3d}: RD={np.mean(rd_list):7.1f}±{np.std(rd_list):5.1f} '
          f'GD={np.mean(gd_list):7.1f} dual={np.mean(dual_list):5.1f}')

sweep = pd.DataFrame(rows)
sweep.to_csv(OUT + 'tier3_campra_sensitivity.csv', index=False)

# ------------------------------------------------------------------ linearity + extrapolation
x = sweep['n_pos'].values
y = sweep['RD_mean'].values
slope, intercept = np.polyfit(x, y, 1)
r_lin = np.corrcoef(x, y)[0, 1]
# crossing where RD = GD reference (1,214) using GD count at k=100 (actual)
gd_actual = base_counts.get('gene-driven', 0)
n_cross = (gd_actual - intercept) / slope if slope > 0 else np.inf
frac_cross = n_cross / n  # fraction of universe that must be assay-positive
print(f'\nlinear fit: RD = {slope:.3f}·n_pos + {intercept:.1f} (r={r_lin:.4f})')
print(f'crossing: RD = GD({gd_actual}) at n_pos* = {n_cross:.0f} genes '
      f'= {100*frac_cross:.1f}% of universe (current caMPRA positives: {n_pos_total} = {100*n_pos_total/n:.1f}%)')
# slope in per-10-positives terms
print(f'slope: each additional 10 assay-positive genes -> {10*slope:.1f} RD genes')

# ------------------------------------------------------------------ ceilings
ceil_all = classify_counts(np.ones(n, dtype=bool))
har_mask = n_hars > 0
ceil_har = classify_counts(har_mask)
print(f"\nceiling (doan=1 ALL {n} genes): RD={ceil_all.get('regulation-driven',0)}, "
      f"GD={ceil_all.get('gene-driven',0)}, dual={ceil_all.get('dual-driven',0)}")
print(f"ceiling (doan=1 for {har_mask.sum()} HAR-proximate genes): RD={ceil_har.get('regulation-driven',0)}, "
      f"GD={ceil_har.get('gene-driven',0)}, dual={ceil_har.get('dual-driven',0)}")

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(figsize=(7.2, 5.4))
ax.errorbar(sweep['n_pos'], sweep['RD_mean'], yerr=sweep['RD_sd'], fmt='o-', color='#2C7FB8',
            capsize=3, lw=1.8, ms=6, label='RD count (mean±SD, 100 reps)')
ax.axhline(gd_actual, color='#D95F02', ls='--', lw=1.8, label=f'GD count = {gd_actual}')
# extrapolation
xx = np.linspace(0, max(n_cross * 1.05, 120), 100)
ax.plot(xx, slope * xx + intercept, color='#2C7FB8', ls=':', lw=1.5, alpha=0.8)
ax.plot([n_cross], [gd_actual], marker='*', ms=18, color='#7B3294', zorder=5)
ax.annotate(f'extrapolated crossing\nn* = {n_cross:.0f} positives\n({100*frac_cross:.0f}% of universe)',
            xy=(n_cross, gd_actual), xytext=(n_cross * 0.55, gd_actual * 0.72),
            fontsize=10, color='#7B3294', arrowprops=dict(arrowstyle='->', color='#7B3294'))
ax.axvspan(0, n_pos_total, alpha=0.08, color='green')
ax.text(n_pos_total / 2, ax.get_ylim()[0] + 30, 'observed range\n(caMPRA, 2% coverage)',
        ha='center', fontsize=9, color='green')
ax.set_xlabel('Number of caMPRA assay-positive genes (regulatory-evidence coverage)')
ax.set_ylabel('Regulation-driven (RD) gene count')
ax.set_title('Tier 3: RD count vs regulatory-assay coverage\n(v7 classifier, caMPRA downsampling)')
ax.legend(loc='upper left', frameon=False)
fig.tight_layout()
fig.savefig(FIG + 'Figure_Tier3_coverage_sensitivity.png', dpi=300)
fig.savefig(FIG + 'Figure_Tier3_coverage_sensitivity.pdf')
print('figure saved')

# ------------------------------------------------------------------ save
result = {
    'n_universe': int(n),
    'n_campra_positives': int(n_pos_total),
    'gd_reference': int(gd_actual),
    'replication_agreement': match,
    'sweep': rows,
    'linear_fit': {'slope': float(slope), 'intercept': float(intercept), 'r': float(r_lin)},
    'crossing': {'n_pos_star': float(n_cross), 'pct_of_universe': float(100 * frac_cross),
                 'fold_over_current': float(n_cross / n_pos_total)},
    'ceilings': {
        'doan_all': {k: int(v) for k, v in ceil_all.items()},
        'doan_har_only': {k: int(v) for k, v in ceil_har.items()},
        'n_har_proximate': int(har_mask.sum()),
    },
}
with open(OUT + 'tier3_campra_sensitivity.json', 'w', encoding='utf-8') as f:
    json.dump(result, f, indent=2, ensure_ascii=False)
print('saved:', OUT + 'tier3_campra_sensitivity.json')
