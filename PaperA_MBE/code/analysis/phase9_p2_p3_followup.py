# -*- coding: utf-8 -*-
"""
P2: GO-slim-style functional category composition of GD / RD / dual classes
============================================================================
Associate-research follow-up #2: classify the 4,974-gene universe into broad
functional categories and test category-specificity of each class.

Approach: keyword-bucket GO-slim over GO_Biological_Process_2023 term names
(Enrichr/GO BP 2023 has no slim mapping in Enrichr; buckets are documented
regexes, overlapping membership allowed). Gene belongs to a bucket if ANY of
its annotated GO BP terms matches. Per class x bucket: Fisher exact (class vs
rest of universe), BH across all tests.

P3: RDS/GDS component decomposition of the classified genes
============================================================
Associate-research follow-up #3: make the "regulatory contribution"
definition concrete — which evidence component drives each regulation-driven
gene, and do disease-anchored RD genes differ?

Outputs:
  results/phase9_hardening/p2_goslim_composition.csv / .json
  results/phase9_hardening/p3_component_decomposition.csv / .json
"""
import json
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

import gseapy

BASE = 'd:/人类正选择基因项目/'
OUT = BASE + 'results/phase9_hardening/'

# ------------------------------------------------------------------ universe
df = pd.read_csv(BASE + 'results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv')
gmap = pd.read_csv(BASE + 'data/gene_id_to_symbol_gencode47.csv')
df = df.merge(gmap, on='gene_id', how='left', suffixes=('', '_g'))
df['sym'] = df['gene_symbol'].where(
    df['gene_symbol'].notna() & (df['gene_symbol'].astype(str).str.strip() != ''),
    df['gene_symbol_g']).astype(str).str.upper().str.strip()
uni_syms = set(df['sym'])
N = len(uni_syms)

# ------------------------------------------------------------------ P2 buckets
BUCKETS = {
    'Immune & defense': r'immune|defense|defen[cs]e|inflammat|interferon|complement|antigen|lymphocyte|leukocyte|cytokine|viral|virulence|wound|phagocyt',
    'Neuronal & synaptic': r'neuro|synap|axon|dendrit|myelin|gliogenesis|dopamin|cholin|gaba|glutamat|seroton|synaptic|phototransduct|circadian|behavior|learning|memory',
    'Development & morphogenesis': r'development|morphogenesis|differentiation|patterning|organogenesis|histogenesis|embryo|axon guidance',
    'Cell cycle, DNA repair & apoptosis': r'cell cycle|mitotic|meiotic|cytokinesis|proliferat|apopt|programmed cell death|dna repair|double.strand|replication|chromosome segregation|genome stability',
    'Metabolism & proteostasis': r'metaboli|biosynth|catabol|lipid|fatty acid|carbohydrate|amino acid|nucleotide|purine|pyrimidine|mitochondr|oxidative phosphorylat|proteolysis|ubiquitin|autophag|glucose|sterol|cholesterol',
    'Signal transduction': r'signal transduct|signaling|signalling|phosphorylat|kinase cascade|g.protein|second messenger|receptor.*pathway|mapk|camp',
    'Transport & secretion': r'transport|ion homeostasis|secretion|endocyt|exocyt|vesicle|efflux|influx|transmembrane import',
    'Transcription, RNA & chromatin': r'transcription|chromatin|histone|epigenetic|rna polymerase|mrna|splicing|ribonucleoprotein|gene silencing|translation|ribose',
    'Reproduction': r'fertiliz|gamet|spermat|oocyte|reproduc|placent|parturition',
    'Multicellular & organ homeostasis': r'homeostasis|circulatory|blood vessel|angio|heart|cardiac|kidney|renal|bone|osteo|epithelium|adhesion|extracellular matrix',
}

lib = gseapy.get_library(name='GO_Biological_Process_2023', organism='human')
term_bucket = {}
for term, genes in lib.items():
    tl = term.lower()
    for b, pat in BUCKETS.items():
        if re.search(pat, tl):
            term_bucket.setdefault(b, set()).update(
                g.upper().strip() for g in genes)

bucket_genes = {b: (s & uni_syms) for b, s in term_bucket.items()}
print('bucket sizes within universe:')
for b, s in sorted(bucket_genes.items(), key=lambda x: -len(x[1])):
    print(f'  {b}: {len(s)}')

classes = {
    'gene_driven': df['classification_v7'] == 'gene-driven',
    'gene_driven_all': df['classification_v7'].isin(['gene-driven', 'gene-driven (relaxed)']),
    'regulation_driven': df['classification_v7'] == 'regulation-driven',
    'dual_driven': df['classification_v7'] == 'dual-driven',
}

rows = []
for cname, mask in classes.items():
    cls_syms = set(df.loc[mask, 'sym'])
    n_cls = len(cls_syms)
    for b, gset in bucket_genes.items():
        K = len(gset)
        if K < 10:
            continue
        k = len(cls_syms & gset)
        rest_k = K - k
        rest_n = N - n_cls
        table = [[k, n_cls - k], [rest_k, rest_n - rest_k]]
        orr, p = stats.fisher_exact(table, alternative='greater')
        # two-sided for depletion direction too
        orr2, p2 = stats.fisher_exact(table)
        rows.append({'class': cname, 'bucket': b, 'K_universe': K, 'n_class': n_cls,
                     'k_overlap': k, 'pct_in_class': round(100 * k / n_cls, 1),
                     'pct_universe': round(100 * K / N, 1),
                     'OR_greater': round(orr, 3) if np.isfinite(orr) else None,
                     'p_greater': p, 'p_two_sided': p2})
tab = pd.DataFrame(rows)
m = len(tab)
order = np.argsort(tab['p_greater'].values)
ps = tab['p_greater'].values[order] * m / (np.arange(m) + 1)
q = np.minimum.accumulate(ps[::-1])[::-1]
tab['q_bh'] = np.nan
tab.loc[tab.index[order], 'q_bh'] = np.minimum(q, 1.0)
tab = tab.sort_values(['class', 'p_greater']).reset_index(drop=True)
tab.to_csv(OUT + 'p2_goslim_composition.csv', index=False)

print('\nP2 composition (pct_in_class vs pct_universe; q_bh):')
for cname in classes:
    sub = tab[tab['class'] == cname]
    line = [f'  {cname}:']
    for _, r in sub.iterrows():
        star = '*' if r['q_bh'] < 0.05 else ''
        line.append(f"    {r['bucket']}: {r['pct_in_class']}% vs {r['pct_universe']}% "
                    f"(k={r['k_overlap']}) q={r['q_bh']:.3f}{star}")
    print('\n'.join(line))

# ------------------------------------------------------------------ P3 decomposition
def component_decomposition(cls_name, mask, comp_cols, weights):
    sub = df[mask]
    contrib = pd.DataFrame({'gene_id': sub['gene_id'], 'gene_symbol': sub['sym'],
                            'classification': sub['classification_v7']})
    for c, w in zip(comp_cols, weights):
        contrib[c] = sub[c].values * w
    total = contrib[comp_cols].sum(axis=1)
    contrib['rds_or_gds_total'] = sub['rds_v7'].values if 'rds' in cls_name.lower() or 'regulation' in cls_name or 'dual' in cls_name else sub['gds_v7'].values
    dominant = contrib[comp_cols].idxmax(axis=1)
    contrib['dominant_component'] = dominant
    contrib['dominant_share'] = (contrib[comp_cols].max(axis=1) / total).round(3)
    return contrib

RDS_COMPS = ['rds_doan', 'rds_tau', 'rds_nc', 'rds_brain']
RDS_W = [0.30, 0.25, 0.25, 0.20]
GDS_COMPS = ['gds_p_pct_resid', 'gds_lrt_pct_resid', 'gds_relax_pct', 'gds_selectome']
GDS_W = [0.35, 0.30, 0.20, 0.15]

p3 = {}
rd_mask = df['classification_v7'] == 'regulation-driven'
gd_mask = df['classification_v7'] == 'gene-driven'
dual_mask = df['classification_v7'] == 'dual-driven'

rd_dec = component_decomposition('regulation_driven', rd_mask, RDS_COMPS, RDS_W)
gd_dec = component_decomposition('gene_driven', gd_mask, GDS_COMPS, GDS_W)
dual_dec = component_decomposition('dual_driven', dual_mask, RDS_COMPS, RDS_W)
rd_dec.to_csv(OUT + 'p3_component_decomposition.csv', index=False)

# dominant component counts
summary = {}
for name, dec in [('regulation_driven', rd_dec), ('gene_driven', gd_dec), ('dual_driven', dual_dec)]:
    vc = dec['dominant_component'].value_counts().to_dict()
    summary[name] = {'n': int(len(dec)), 'dominant_component_counts': vc,
                     'mean_dominant_share': round(float(dec['dominant_share'].mean()), 3)}
    print(f"\nP3 {name} (n={len(dec)}) dominant components: {vc}")
    print(f"  mean dominant share: {dec['dominant_share'].mean():.3f}")

# GWAS-anchored RD genes: which components dominate? (EA trait as example)
t2 = pd.read_csv(OUT + 'tier2_gwas_enrichment.csv')
ea = t2[(t2['trait'] == 'Educational attainment') & (t2['class'] == 'regulation_driven')]
# reconstruct EA gene set
EA_KEYS = ('educational attainment',)
token_re = re.compile(r'[A-Za-z][A-Za-z0-9@.\-]*')
ea_genes = set()
with open(BASE + 'data/gwas/gwas-catalog-download-associations-v1.0-full.tsv', encoding='utf-8', errors='replace') as f:
    header = f.readline().rstrip('\n').split('\t')
    ix = {c: i for i, c in enumerate(header)}
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < len(header):
            continue
        trait = parts[ix['DISEASE/TRAIT']].lower()
        if 'educational attainment' not in trait:
            continue
        try:
            p = float(parts[ix['P-VALUE']])
        except ValueError:
            continue
        if p < 5e-8:
            for tok in token_re.findall(parts[ix['MAPPED_GENE']]):
                t = tok.upper().strip('.-')
                if t in uni_syms:
                    ea_genes.add(t)
rd_dec['is_ea_gwas'] = rd_dec['gene_symbol'].astype(str).str.upper().isin(ea_genes) | rd_dec['gene_symbol'].isin(ea_genes)
ea_dec = rd_dec[rd_dec['is_ea_gwas']]
ne_dec = rd_dec[~rd_dec['is_ea_gwas']]
summary['rd_ea_gwas'] = {'n': int(len(ea_dec)),
                         'dominant_component_counts': ea_dec['dominant_component'].value_counts().to_dict(),
                         'mean_component_contrib': {c: round(float(ea_dec[c].mean()), 4) for c in RDS_COMPS}}
summary['rd_non_ea'] = {'n': int(len(ne_dec)),
                        'dominant_component_counts': ne_dec['dominant_component'].value_counts().to_dict(),
                        'mean_component_contrib': {c: round(float(ne_dec[c].mean()), 4) for c in RDS_COMPS}}
print(f"\nP3 RD ∩ EA-GWAS (n={len(ea_dec)}): dominant {summary['rd_ea_gwas']['dominant_component_counts']}")
print(f"  mean contrib: {summary['rd_ea_gwas']['mean_component_contrib']}")
print(f"P3 RD non-EA (n={len(ne_dec)}): dominant {summary['rd_non_ea']['dominant_component_counts']}")
print(f"  mean contrib: {summary['rd_non_ea']['mean_component_contrib']}")

# chi-square: dominant component differs between EA+ and EA- RD genes?
from scipy.stats import chi2_contingency
cats = sorted(set(rd_dec['dominant_component']))
tab_ea = [[int((ea_dec['dominant_component'] == c).sum()) for c in cats],
          [int((ne_dec['dominant_component'] == c).sum()) for c in cats]]
chi2, p_chi, dof, _ = chi2_contingency(tab_ea)
summary['dominant_component_ea_vs_non_ea_chi2'] = {'chi2': float(chi2), 'p': float(p_chi), 'dof': int(dof)}
print(f"  dominant-component EA+ vs EA- chi2={chi2:.2f} p={p_chi:.3g}")

summary['p2_bucket_sizes_universe'] = {b: len(s) for b, s in bucket_genes.items()}
with open(OUT + 'p3_component_decomposition.json', 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print('\nsaved p2_goslim_composition.csv and p3_component_decomposition.{csv,json}')
