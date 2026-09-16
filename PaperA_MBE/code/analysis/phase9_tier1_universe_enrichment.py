# -*- coding: utf-8 -*-
"""
Tier 1: Custom-background (4,974-gene universe) enrichment re-testing
=====================================================================
Problem (associate researcher critique): original phase6 Enrichr used the FULL
human gene background (~18k genes), not our 4,974 primate-ortholog universe.
Full-background p-values are conservative but the correct statistical frame
conditions on the universe actually assayed.

Fix: for every term in GO_Biological_Process_2023 and SynGO_2024, recompute
hypergeometric enrichment with:
  N = |universe genes with symbols|
  K = |term genes within universe|
  n = |class gene set|
  k = |class ∩ term|
BH across all terms tested per (class, library). Plus 1,000 label-permutation
empirical P for top terms.

Classes (v7): gene-driven (1214), regulation-driven (293), dual-driven (11),
plus sensitivity 'gene-driven+relaxed' (1266).

Output:
  results/phase9_hardening/tier1_universe_enrichment_{class}_{lib}.csv
  results/phase9_hardening/tier1_universe_enrichment_summary.json
"""
import json
import os

import gseapy
import numpy as np
import pandas as pd
from scipy import stats

BASE = 'd:/人类正选择基因项目/'
OUT = BASE + 'results/phase9_hardening/'
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ universe
df = pd.read_csv(BASE + 'results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv')
# recover symbols for ENSG-only genes from GENCODE v47 (929/929 recovered)
gmap = pd.read_csv(BASE + 'data/gene_id_to_symbol_gencode47.csv')
df = df.merge(gmap, on='gene_id', how='left', suffixes=('', '_g'))
df['sym_final'] = df['gene_symbol'].where(
    df['gene_symbol'].notna() & (df['gene_symbol'].astype(str).str.strip() != ''),
    df['gene_symbol_g'])
sym_ok = df['sym_final'].notna() & (df['sym_final'].astype(str).str.strip() != '')
uni = df[sym_ok].copy()
uni['sym'] = uni['sym_final'].astype(str).str.upper().str.strip()
uni_syms = set(uni['sym'])
N = len(uni_syms)
print(f'universe with symbols: {N} / {len(df)} (dropped {len(df)-N} unnamed; GENCODE recovered 929)')

classes = {
    'gene_driven': uni['classification_v7'] == 'gene-driven',
    'regulation_driven': uni['classification_v7'] == 'regulation-driven',
    'dual_driven': uni['classification_v7'] == 'dual-driven',
    'gene_driven_all': uni['classification_v7'].isin(['gene-driven', 'gene-driven (relaxed)']),
}

libs = {'GO_BP': 'GO_Biological_Process_2023', 'SynGO': 'SynGO_2024'}
MIN_TERM = {'GO_BP': 5, 'SynGO': 3}

summary = {'universe_n': N, 'universe_with_symbol': N, 'dropped_unnamed': int(len(df) - N)}

# ------------------------------------------------------------------ helpers
def hypergeometric_p(k, K, n, N):
    """P(X >= k) for term enrichment."""
    if k == 0:
        return 1.0
    return float(stats.hypergeom.sf(k - 1, N, K, n))

def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    return out

rng = np.random.default_rng(20260915)

for cname, mask in classes.items():
    cls_syms = set(uni.loc[mask, 'sym'])
    n_cls = len(cls_syms)
    print(f'\n=== {cname}: n={n_cls} ===')
    summary[cname] = {'n': n_cls}
    for lname, libname in libs.items():
        try:
            lib = gseapy.get_library(name=libname, organism='human')
        except Exception as e:
            print(f'  {lname}: library FAIL {e!r}')
            continue
        rows = []
        for term, genes in lib.items():
            tset = set(g.upper().strip() for g in genes) & uni_syms
            K = len(tset)
            if K < MIN_TERM[lname]:
                continue
            k = len(cls_syms & tset)
            p = hypergeometric_p(k, K, n_cls, N)
            # fold enrichment
            fe = (k / n_cls) / (K / N) if k > 0 else 0.0
            rows.append({'term': term, 'K_universe': K, 'k_overlap': k,
                         'fold_enrichment': round(fe, 3), 'p': p})
        tab = pd.DataFrame(rows)
        if tab.empty:
            continue
        tab['p_bh'] = bh(tab['p'].values)
        tab = tab.sort_values('p').reset_index(drop=True)
        n_sig = int((tab['p_bh'] < 0.05).sum())
        n_nom = int((tab['p'] < 0.05).sum())
        print(f'  {lname}: {len(tab)} terms tested, nominal p<0.05: {n_nom}, BH FDR<0.05: {n_sig}')
        summary[cname][lname] = {'terms_tested': int(len(tab)), 'nominal_005': n_nom, 'bh_005': n_sig}

        # permutation empirical P for top 10 terms
        top = tab.head(10)
        perm_res = []
        n_perm = 1000
        for _, row in top.iterrows():
            tset = set(g.upper().strip() for g in lib[row['term']]) & uni_syms
            K = len(tset)
            k_obs = row['k_overlap']
            p_obs = row['p']
            # permutation: draw n_cls random genes from universe, count overlap >= k_obs
            uni_list = np.array(sorted(uni_syms))
            tarr = np.zeros(len(uni_list), dtype=bool)
            for i, s in enumerate(uni_list):
                tarr[i] = s in tset
            hits = 0
            for _ in range(n_perm):
                idx = rng.choice(len(uni_list), size=n_cls, replace=False)
                if tarr[idx].sum() >= k_obs:
                    hits += 1
            p_emp = (hits + 1) / (n_perm + 1)
            perm_res.append({'term': row['term'], 'k': int(k_obs), 'K': int(K),
                             'p_hyper': p_obs, 'p_bh': row['p_bh'], 'p_perm': round(p_emp, 5)})
        tab['p_perm'] = tab['term'].map({r['term']: r['p_perm'] for r in perm_res})
        tab.to_csv(OUT + f'tier1_universe_enrichment_{cname}_{lname}.csv', index=False)
        summary[cname][lname]['top10_perm'] = perm_res
        for r in perm_res[:5]:
            print(f"    {r['term'][:55]:55s} k={r['k']:3d}/{r['K']:4d} p={r['p_hyper']:.2e} "
                  f"BH={r['p_bh']:.2e} perm={r['p_perm']:.4f}")

with open(OUT + 'tier1_universe_enrichment_summary.json', 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2, ensure_ascii=False, default=float)
print('\nSaved summary:', OUT + 'tier1_universe_enrichment_summary.json')
