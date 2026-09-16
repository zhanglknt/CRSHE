# -*- coding: utf-8 -*-
"""
Tier 2: GWAS Catalog enrichment — external non-circular anchor
==============================================================
Question: do gene-driven / regulation-driven / dual classes show differential
enrichment for neuropsychiatric and cognitive trait loci (SCZ, ASD, EA,
intelligence, MDD, bipolar, neuroticism), against the 4,974-gene universe?

Data: GWAS Catalog full associations (v1.0), genome-wide significant
associations only (P < 5e-8), MAPPED_GENE symbols.
Background: the 4,974-gene universe itself (external to all class definitions
=> non-circular).
Tests: Fisher exact (class vs rest of universe) per trait x class, BH across
all tests.

Output:
  results/phase9_hardening/tier2_gwas_enrichment.csv
  results/phase9_hardening/tier2_gwas_enrichment.json
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
df = pd.read_csv(BASE + 'results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv')
gmap = pd.read_csv(BASE + 'data/gene_id_to_symbol_gencode47.csv')
df = df.merge(gmap, on='gene_id', how='left', suffixes=('', '_g'))
df['sym_final'] = df['gene_symbol'].where(
    df['gene_symbol'].notna() & (df['gene_symbol'].astype(str).str.strip() != ''),
    df['gene_symbol_g'])
df['sym'] = df['sym_final'].astype(str).str.upper().str.strip()
uni_syms = set(df['sym'])
N = len(uni_syms)

classes = {
    'gene_driven': df['classification_v7'] == 'gene-driven',
    'gene_driven_all': df['classification_v7'].isin(['gene-driven', 'gene-driven (relaxed)']),
    'regulation_driven': df['classification_v7'] == 'regulation-driven',
    'dual_driven': df['classification_v7'] == 'dual-driven',
}

trait_groups = {
    'Schizophrenia': ('schizophren',),
    'ASD': ('autis',),
    'Educational attainment': ('educational attainment',),
    'Intelligence': ('intelligence',),
    'Major depression': ('major depress',),
    'Bipolar disorder': ('bipolar',),
    'Neuroticism': ('neuroticism',),
    'Cognitive performance': ('cognitive performance',),
}

token_re = re.compile(r'[A-Za-z][A-Za-z0-9@.\-]*')

# ------------------------------------------------------------------ scan
trait_genes = {k: set() for k in trait_groups}
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
        matched = False
        for gname, keys in trait_groups.items():
            if any(k in trait for k in keys):
                matched = True
        if not matched:
            continue
        # parse p-value
        try:
            p = float(parts[ix['P-VALUE']])
        except (ValueError, IndexError):
            continue
        if not p < 5e-8:
            continue
        n_gws += 1
        mg = parts[ix['MAPPED_GENE']]
        for tok in token_re.findall(mg):
            t = tok.upper().strip('.-')
            if t in uni_syms:
                for gname, keys in trait_groups.items():
                    if any(k in trait for k in keys):
                        trait_genes[gname].add(t)

print(f'scanned {n_rows} rows, {n_gws} GWS rows for target traits')
for g, s in trait_genes.items():
    print(f'  {g}: {len(s)} universe genes')

# ------------------------------------------------------------------ tests
rows = []
for gname, gset in trait_genes.items():
    K = len(gset)
    if K == 0:
        continue
    for cname, mask in classes.items():
        cls_syms = set(df.loc[mask, 'sym'])
        n_cls = len(cls_syms)
        k = len(cls_syms & gset)
        rest_k = K - k
        rest_n = N - n_cls
        table = [[k, n_cls - k], [rest_k, rest_n - rest_k]]
        orr, p = stats.fisher_exact(table, alternative='greater')
        rows.append({'trait': gname, 'class': cname, 'K_universe': K, 'n_class': n_cls,
                     'k_overlap': k, 'expected': round(K * n_cls / N, 1),
                     'OR': round(orr, 3) if np.isfinite(orr) else None, 'p': p})

tab = pd.DataFrame(rows)
# BH across all tests
m = len(tab)
order = np.argsort(tab['p'].values)
p_sorted = tab['p'].values[order] * m / (np.arange(m) + 1)
q = np.minimum.accumulate(p_sorted[::-1])[::-1]
tab['q_bh'] = np.nan
tab.loc[tab.index[order], 'q_bh'] = np.minimum(q, 1.0)
tab = tab.sort_values('p').reset_index(drop=True)
tab.to_csv(OUT + 'tier2_gwas_enrichment.csv', index=False)

print('\n=== enrichment (sorted by p) ===')
for _, r in tab.head(20).iterrows():
    print(f"  {r['trait']:24s} {r['class']:20s} k={r['k_overlap']:3d}/{r['n_class']:4d} "
          f"(exp {r['expected']:.0f}, K={r['K_universe']}) OR={r['OR']} p={r['p']:.2e} q={r['q_bh']:.3f}")

with open(OUT + 'tier2_gwas_enrichment.json', 'w', encoding='utf-8') as f:
    json.dump({'universe_n': int(N), 'n_gws_rows': int(n_gws),
               'trait_gene_counts': {k: len(v) for k, v in trait_genes.items()},
               'results': tab.to_dict(orient='records')}, f, indent=2, ensure_ascii=False, default=float)
print('\nsaved:', OUT + 'tier2_gwas_enrichment.json')
