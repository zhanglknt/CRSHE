# -*- coding: utf-8 -*-
"""
P0-A Step 3+4: Tier 2 (GWAS Catalog enrichment) re-run under LOO variants
+ non-neurological negative-control traits
==============================================================================
Part A (LOO re-run): repeat phase9_tier2_gwas_enrichment.py methodology under
LOO-brain and LOO-tau classifications (8 neuro traits x 4 classes, Fisher
exact vs rest of universe, BH across 8x4=32 tests per variant).

Part B (negative controls): 4 non-neuro traits (Height, LDL cholesterol,
Crohn disease, Type 2 diabetes) tested under v7 (full), LOO-brain, LOO-tau
classifications. Per classification, BH across the joint family of
(8 neuro + 4 control) traits x 4 classes = 48 tests. Expectation: RD enriched
in the 8 neuro traits but NOT in the 4 controls.

Data: GWAS Catalog full associations v1.0, P < 5e-8, MAPPED_GENE symbols.

Output:
  results/phase9_hardening/loo_tier2_loo-brain.csv / .json
  results/phase9_hardening/loo_tier2_loo-tau.csv / .json
  results/phase9_hardening/tier2_negative_controls.csv / .json
"""
import json
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

BASE = 'd:/人类正选择基因项目/'
GWAS = BASE + 'data/gwas/gwas-catalog-download-associations-v1.0-full.tsv'
OUT = BASE + 'results/phase9_hardening/'

# ------------------------------------------------------------------ universe
pg = pd.read_csv(OUT + 'loo_variant_per_gene.csv')
gmap = pd.read_csv(BASE + 'data/gene_id_to_symbol_gencode47.csv')
pg = pg.merge(gmap, on='gene_id', how='left', suffixes=('', '_g'))
pg['sym_final'] = pg['gene_symbol'].where(
    pg['gene_symbol'].notna() & (pg['gene_symbol'].astype(str).str.strip() != ''),
    pg['gene_symbol_g'])
pg['sym'] = pg['sym_final'].astype(str).str.upper().str.strip()
uni_syms = set(pg['sym'])
N = len(uni_syms)
print(f'universe: {N} genes with symbols')

neuro_traits = {
    'Schizophrenia': ('schizophren',),
    'ASD': ('autis',),
    'Educational attainment': ('educational attainment',),
    'Intelligence': ('intelligence',),
    'Major depression': ('major depress',),
    'Bipolar disorder': ('bipolar',),
    'Neuroticism': ('neuroticism',),
    'Cognitive performance': ('cognitive performance',),
}
control_traits = {
    'Height': ('height',),
    'LDL cholesterol': ('ldl cholesterol', 'low density lipoprotein cholesterol'),
    'Crohn disease': ('crohn',),
    'Type 2 diabetes': ('type 2 diabetes',),
}
all_traits = {**neuro_traits, **control_traits}
is_control = {t: t in control_traits for t in all_traits}

token_re = re.compile(r'[A-Za-z][A-Za-z0-9@.\-]*')

# ------------------------------------------------------------------ single GWAS scan
trait_genes = {k: set() for k in all_traits}
trait_counts = {k: set() for k in all_traits}  # distinct reported trait strings
n_rows = 0
n_gws = 0
with open(GWAS, encoding='utf-8', errors='replace') as f:
    header = f.readline().rstrip('\n').split('\t')
    ix = {c: i for i, c in enumerate(header)}
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < len(header):
            continue
        n_rows += 1
        trait = parts[ix['DISEASE/TRAIT']].lower()
        matched = [gname for gname, keys in all_traits.items()
                   if any(k in trait for k in keys)]
        if not matched:
            continue
        try:
            p = float(parts[ix['P-VALUE']])
        except (ValueError, IndexError):
            continue
        if not p < 5e-8:
            continue
        n_gws += 1
        mg = parts[ix['MAPPED_GENE']]
        toks = set()
        for tok in token_re.findall(mg):
            t = tok.upper().strip('.-')
            if t in uni_syms:
                toks.add(t)
        for gname in matched:
            trait_counts[gname].add(parts[ix['DISEASE/TRAIT']])
            trait_genes[gname].update(toks)

print(f'scanned {n_rows} rows, {n_gws} GWS rows for target traits')
for g in all_traits:
    print(f'  {g}: {len(trait_genes[g])} universe genes '
          f'({len(trait_counts[g])} distinct reported traits)')


def bh_q(pvals):
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    return out


def run_fisher(trait_set, cls_col):
    """Fisher exact (class vs rest of universe) for every trait x class.
    trait_set: dict trait name -> set of universe gene symbols."""
    classes = {
        'gene_driven': pg[cls_col] == 'gene-driven',
        'gene_driven_all': pg[cls_col].isin(['gene-driven', 'gene-driven (relaxed)']),
        'regulation_driven': pg[cls_col] == 'regulation-driven',
        'dual_driven': pg[cls_col] == 'dual-driven',
    }
    rows = []
    for gname, gset in trait_set.items():
        K = len(gset)
        if K == 0:
            continue
        for cname, mask in classes.items():
            cls_syms = set(pg.loc[mask, 'sym'])
            n_cls = len(cls_syms)
            k = len(cls_syms & gset)
            rest_k = K - k
            rest_n = N - n_cls
            table = [[k, n_cls - k], [rest_k, rest_n - rest_k]]
            orr, p = stats.fisher_exact(table, alternative='greater')
            rows.append({'trait': gname, 'is_control': is_control.get(gname, False),
                         'class': cname, 'K_universe': K, 'n_class': n_cls,
                         'k_overlap': k, 'expected': round(K * n_cls / N, 1),
                         'OR': round(orr, 3) if np.isfinite(orr) else None, 'p': p})
    return pd.DataFrame(rows)


def bh_apply(tab):
    tab = tab.copy().reset_index(drop=True)
    tab['q_bh'] = bh_q(tab['p'].values)
    return tab.sort_values('p').reset_index(drop=True)


# ------------------------------------------------------------------ Part A: LOO re-run
for variant in ['LOO-brain', 'LOO-tau']:
    neuro_sets = {k: trait_genes[k] for k in neuro_traits}
    tab = run_fisher(neuro_sets, f'cls_{variant}')
    tab = bh_apply(tab)
    tab.to_csv(OUT + f'loo_tier2_{variant}.csv', index=False)
    m = len(tab)
    print(f'\n=== Tier2 {variant} (8 neuro traits x 4 classes, {m} tests, BH within) ===')
    for _, r in tab.iterrows():
        print(f"  {r['trait']:24s} {r['class']:20s} k={r['k_overlap']:3d}/{r['n_class']:4d} "
              f"(exp {r['expected']:.0f}) OR={r['OR']} p={r['p']:.2e} q={r['q_bh']:.3f}")
    with open(OUT + f'loo_tier2_{variant}.json', 'w', encoding='utf-8') as f:
        json.dump({'universe_n': int(N), 'variant': variant, 'n_tests': int(m),
                   'trait_gene_counts': {k: len(trait_genes[k]) for k in neuro_traits},
                   'results': tab.to_dict(orient='records')}, f, indent=2,
                  ensure_ascii=False, default=float)
    print('saved:', OUT + f'loo_tier2_{variant}.json')

# ------------------------------------------------------------------ Part B: negative controls
neg_summary = {}
neg_rows_all = []
for variant, cls_col in [('full_v7', 'cls_full'), ('LOO-brain', 'cls_LOO-brain'),
                         ('LOO-tau', 'cls_LOO-tau')]:
    all_sets = {k: trait_genes[k] for k in all_traits}
    tab = run_fisher(all_sets, cls_col)             # 12 traits x 4 classes = 48 tests
    tab = bh_apply(tab)                             # BH within classification
    tab.insert(0, 'classification', variant)
    neg_rows_all.append(tab)
    rd_neuro = tab[(tab['class'] == 'regulation_driven') & (~tab['is_control'])]
    rd_ctrl = tab[(tab['class'] == 'regulation_driven') & (tab['is_control'])]
    print(f"\n=== negative controls under {variant} ===")
    print('  RD x 8 neuro traits: significant (q<0.05): '
          f"{int((rd_neuro['q_bh'] < 0.05).sum())}/{len(rd_neuro)}")
    print('  RD x 4 control traits: significant (q<0.05): '
          f"{int((rd_ctrl['q_bh'] < 0.05).sum())}/{len(rd_ctrl)}")
    for _, r in tab[tab['is_control']].iterrows():
        print(f"  CTRL {r['trait']:18s} {r['class']:20s} k={r['k_overlap']:3d}/"
              f"{r['n_class']:4d} (exp {r['expected']:.0f}) OR={r['OR']} "
              f"p={r['p']:.3f} q={r['q_bh']:.3f}")
    neg_summary[variant] = {
        'trait_gene_counts': {k: len(trait_genes[k]) for k in all_traits},
        'n_tests': int(len(tab)),
        'RD_neuro_sig_q05': int((rd_neuro['q_bh'] < 0.05).sum()),
        'RD_neuro_total': int(len(rd_neuro)),
        'RD_control_sig_q05': int((rd_ctrl['q_bh'] < 0.05).sum()),
        'RD_control_total': int(len(rd_ctrl)),
        'results': tab.to_dict(orient='records'),
    }

neg_tab = pd.concat(neg_rows_all, ignore_index=True)
neg_tab.to_csv(OUT + 'tier2_negative_controls.csv', index=False)
with open(OUT + 'tier2_negative_controls.json', 'w', encoding='utf-8') as f:
    json.dump({'universe_n': int(N),
               'control_traits': list(control_traits.keys()),
               'gwas_scan': {'n_rows_scanned': int(n_rows), 'n_gws_rows_target_traits': int(n_gws)},
               'trait_gene_counts': {k: len(v) for k, v in trait_genes.items()},
               'trait_distinct_reported': {k: len(v) for k, v in trait_counts.items()},
               'by_classification': neg_summary}, f, indent=2,
              ensure_ascii=False, default=float)
print('\nsaved:', OUT + 'tier2_negative_controls.csv')
print('saved:', OUT + 'tier2_negative_controls.json')
