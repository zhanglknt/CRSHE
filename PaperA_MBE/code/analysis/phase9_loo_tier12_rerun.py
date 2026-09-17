# -*- coding: utf-8 -*-
"""
P0-A Step 1: Recompute LOO-variant per-gene classifications from scratch
=========================================================================
Round-5 review issue: Tier 1/2 enrichments were never re-run under LOO
classification variants. This script rebuilds the LOO per-gene classifications
(full, LOO-caMPRA, LOO-tau, LOO-nc, LOO-brain) from the documented RDS weight
formula, validates them against results/paper/revision_v7/loo_full_matrix.csv,
and writes a per-gene table for downstream Tier 1/2 re-runs.

RDS (documented): rds = 0.30*rds_doan + 0.25*rds_tau + 0.25*rds_nc + 0.20*rds_brain
GDS unchanged (gds_v7). Classification rules identical to v7 (threshold 0.15):
  RD   = rds > gds + 0.15 AND NOT BUSTED-significant
  dual = rds > gds + 0.15 AND BUSTED-significant AND gds > 0.3
  GD   = gds > rds + 0.15 AND BUSTED-significant
  precedence: GD > RD > dual; RELAXED flag converts gene-driven -> gene-driven (relaxed)

LOO variant weights: canonical published weights that produced
results/paper/revision_v7/loo_full_matrix.csv (paperA_revision_v7_analysis.py):
  LOO-caMPRA: tau=0.35, nc=0.35, brain=0.30   (redistributed, rounded)
  LOO-tau:    doan=0.40, nc=0.33, brain=0.27  (approx. proportional /0.75, rounded)
  LOO-nc:     doan=0.40, tau=0.33, brain=0.27 (approx. proportional /0.80, rounded)
  LOO-brain:  doan=0.375, tau=0.3125, nc=0.3125 (exact proportional /0.80)

Output:
  results/phase9_hardening/loo_variant_per_gene.csv
  results/phase9_hardening/loo_variant_check.json
"""
import json
import os

import numpy as np
import pandas as pd

BASE = 'd:/人类正选择基因项目/'
OUT = BASE + 'results/phase9_hardening/'
os.makedirs(OUT, exist_ok=True)

THRESH = 0.15
W = {'doan': 0.30, 'tau': 0.25, 'nc': 0.25, 'brain': 0.20}

df = pd.read_csv(BASE + 'results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv')
n = len(df)

gds = df['gds_v7'].values
sig = df['bh_fdr_sig'].values.astype(bool)
relaxed = df['relax_relaxed'].values.astype(bool)
doan = df['has_doan_campra'].fillna(False).values.astype(bool)

comp = {
    'doan': doan.astype(float),
    'tau': df['rds_tau'].values,
    'nc': df['rds_nc'].values,
    'brain': df['rds_brain'].values,
}


def classify(rds):
    gd = (gds > rds + THRESH) & sig
    rd = (rds > gds + THRESH) & (~sig)
    dual = (rds > gds + THRESH) & sig & (gds > 0.3)
    cls = np.where(gd, 'gene-driven',
                   np.where(rd, 'regulation-driven',
                            np.where(dual, 'dual-driven', 'neutral')))
    cls = np.where(relaxed & (cls == 'gene-driven'), 'gene-driven (relaxed)', cls)
    return cls


def rds_for(weights):
    """weights: dict of component -> weight (None = dropped)."""
    rds = np.zeros(n)
    for c, w in weights.items():
        if w is not None:
            rds += w * comp[c]
    return rds


variants = {
    'full': {'doan': 0.30, 'tau': 0.25, 'nc': 0.25, 'brain': 0.20},
    'LOO-caMPRA': {'doan': None, 'tau': 0.35, 'nc': 0.35, 'brain': 0.30},
    'LOO-tau': {'doan': 0.40, 'tau': None, 'nc': 0.33, 'brain': 0.27},
    'LOO-nc': {'doan': 0.40, 'tau': 0.33, 'nc': None, 'brain': 0.27},
    'LOO-brain': {'doan': 0.375, 'tau': 0.3125, 'nc': 0.3125, 'brain': None},
}

per_gene = pd.DataFrame({'gene_id': df['gene_id'], 'gene_symbol': df['gene_symbol'],
                         'sym_final_v7': df['gene_symbol']})
counts = {}
saved_cols = {'full': 'classification_v7', 'LOO-caMPRA': 'classification_loo_campra',
              'LOO-tau': 'classification_loo_tau', 'LOO-brain': 'classification_loo_brain'}
check = {}

for vname, keep in variants.items():
    rds = rds_for(keep)
    cls = classify(rds)
    per_gene[f'cls_{vname}'] = cls
    per_gene[f'rds_{vname}'] = rds
    vals, cnts = np.unique(cls, return_counts=True)
    counts[vname] = {k: int(v) for k, v in zip(vals, cnts)}
    if vname in saved_cols:
        agree = float((cls == df[saved_cols[vname]].values).mean())
        check[vname] = {'agreement_with_saved_column': round(agree, 6),
                        'match': bool(agree > 0.999)}

# ---- validation against loo_full_matrix.csv (documented summary) ----
loo_mat = pd.read_csv(BASE + 'results/paper/revision_v7/loo_full_matrix.csv')
matrix_ok = True
matrix_detail = {}
for _, row in loo_mat.iterrows():
    v = row['variant']
    if v not in counts:
        continue
    exp = {k: int(row[k]) for k in ['gene-driven', 'gene-driven (relaxed)',
                                    'regulation-driven', 'dual-driven', 'neutral']}
    got = counts[v]
    ok = all(got.get(k, 0) == vv for k, vv in exp.items())
    matrix_detail[v] = {'expected': exp, 'got': got, 'match': bool(ok)}
    if not ok:
        matrix_ok = False

print('=== per-variant counts (recomputed) ===')
for v, c in counts.items():
    print(v, c)
print('\n=== agreement with saved per-gene columns ===')
for v, c in check.items():
    print(v, c)
print('\n=== validation vs loo_full_matrix.csv ===')
for v, d in matrix_detail.items():
    print(v, 'MATCH' if d['match'] else 'MISMATCH', d['expected'], '->', d['got'])
print('\nALL MATRIX CHECKS PASS:', matrix_ok)

per_gene.to_csv(OUT + 'loo_variant_per_gene.csv', index=False)

result = {
    'n_universe': int(n),
    'weights_full': W,
    'loo_weights': {v: {c: w for c, w in ws.items()} for v, ws in variants.items()},
    'loo_weights_source': 'canonical published weights from paperA_revision_v7_analysis.py '
                          '(produce loo_full_matrix.csv); LOO-tau/nc are rounded, '
                          'not exact proportional renormalization',
    'counts': counts,
    'agreement_with_saved_columns': check,
    'loo_full_matrix_validation': matrix_detail,
    'loo_full_matrix_all_match': bool(matrix_ok),
}
with open(OUT + 'loo_variant_check.json', 'w', encoding='utf-8') as f:
    json.dump(result, f, indent=2, ensure_ascii=False)
print('\nsaved:', OUT + 'loo_variant_per_gene.csv')
print('saved:', OUT + 'loo_variant_check.json')
