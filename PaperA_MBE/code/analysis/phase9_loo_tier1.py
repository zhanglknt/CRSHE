# -*- coding: utf-8 -*-
"""
P0-A Step 2: Tier 1 (custom-universe GO BP / SynGO enrichment) re-run under
LOO classification variants (LOO-brain, LOO-tau)
==============================================================================
Round-5 review issue: Tier 1 enrichments were only run under the v7 (full)
classification. This script repeats the exact phase9_tier1_universe_enrichment.py
methodology for the LOO-brain and LOO-tau per-gene classifications:
  - hypergeometric P(X>=k) vs 4,974-gene universe (GENCODE v47 sym_final)
  - terms with >=5 (GO BP) / >=3 (SynGO) universe genes
  - BH per (variant, class, library) across all terms
  - 1,000 label permutations for top-10 terms per cell

Output:
  results/phase9_hardening/loo_tier1_{variant}_{class}_{lib}.csv
  results/phase9_hardening/loo_tier1_summary.json
"""
import json
import os

import gseapy
import numpy as np
import pandas as pd
from scipy import stats

BASE = 'd:/人类正选择基因项目/'
OUT = BASE + 'results/phase9_hardening/'

# ------------------------------------------------------------------ universe
pg = pd.read_csv(OUT + 'loo_variant_per_gene.csv')
gmap = pd.read_csv(BASE + 'data/gene_id_to_symbol_gencode47.csv')
pg = pg.merge(gmap, on='gene_id', how='left', suffixes=('', '_g'))
pg['sym_final'] = pg['gene_symbol'].where(
    pg['gene_symbol'].notna() & (pg['gene_symbol'].astype(str).str.strip() != ''),
    pg['gene_symbol_g'])
sym_ok = pg['sym_final'].notna() & (pg['sym_final'].astype(str).str.strip() != '')
uni = pg[sym_ok].copy()
uni['sym'] = uni['sym_final'].astype(str).str.upper().str.strip()
uni_syms = set(uni['sym'])
N = len(uni_syms)
print(f'universe with symbols: {N} / {len(pg)}')

VARIANTS = ['LOO-brain', 'LOO-tau']
libs = {'GO_BP': 'GO_Biological_Process_2023', 'SynGO': 'SynGO_2024'}
MIN_TERM = {'GO_BP': 5, 'SynGO': 3}

summary = {'universe_n': N, 'variants': VARIANTS}

def hypergeometric_p(k, K, n, N):
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
uni_list = np.array(sorted(uni_syms))

lib_cache = {}
for lname, libname in libs.items():
    lib_cache[lname] = gseapy.get_library(name=libname, organism='human')

for variant in VARIANTS:
    col = f'cls_{variant}'
    classes = {
        'gene_driven': uni[col] == 'gene-driven',
        'gene_driven_all': uni[col].isin(['gene-driven', 'gene-driven (relaxed)']),
        'regulation_driven': uni[col] == 'regulation-driven',
        'dual_driven': uni[col] == 'dual-driven',
    }
    summary[variant] = {c: {'n': int(m.sum())} for c, m in classes.items()}
    for cname, mask in classes.items():
        cls_syms = set(uni.loc[mask, 'sym'])
        n_cls = len(cls_syms)
        if n_cls == 0:
            continue
        print(f'\n=== {variant} / {cname}: n={n_cls} ===')
        for lname in libs:
            lib = lib_cache[lname]
            rows = []
            for term, genes in lib.items():
                tset = set(g.upper().strip() for g in genes) & uni_syms
                K = len(tset)
                if K < MIN_TERM[lname]:
                    continue
                k = len(cls_syms & tset)
                p = hypergeometric_p(k, K, n_cls, N)
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
            print(f'  {lname}: {len(tab)} terms, nominal p<0.05: {n_nom}, BH FDR<0.05: {n_sig}')
            summary[variant][cname][lname] = {'terms_tested': int(len(tab)),
                                              'nominal_005': n_nom, 'bh_005': n_sig}
            # permutation empirical P for top 10 terms
            top = tab.head(10)
            perm_res = []
            n_perm = 1000
            for _, row in top.iterrows():
                tset = set(g.upper().strip() for g in lib[row['term']]) & uni_syms
                k_obs = row['k_overlap']
                p_obs = row['p']
                tarr = np.array([s in tset for s in uni_list])
                hits = 0
                for _ in range(n_perm):
                    idx = rng.choice(len(uni_list), size=n_cls, replace=False)
                    if tarr[idx].sum() >= k_obs:
                        hits += 1
                p_emp = (hits + 1) / (n_perm + 1)
                perm_res.append({'term': row['term'], 'k': int(k_obs),
                                 'K': int(len(tset)), 'p_hyper': p_obs,
                                 'p_bh': row['p_bh'], 'p_perm': round(p_emp, 5)})
            tab['p_perm'] = tab['term'].map({r['term']: r['p_perm'] for r in perm_res})
            tab.to_csv(OUT + f'loo_tier1_{variant}_{cname}_{lname}.csv', index=False)
            summary[variant][cname][lname]['top10_perm'] = perm_res
            for r in perm_res[:5]:
                print(f"    {r['term'][:55]:55s} k={r['k']:3d}/{r['K']:4d} "
                      f"p={r['p_hyper']:.2e} BH={r['p_bh']:.2e} perm={r['p_perm']:.4f}")

with open(OUT + 'loo_tier1_summary.json', 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2, ensure_ascii=False, default=float)
print('\nsaved:', OUT + 'loo_tier1_summary.json')
