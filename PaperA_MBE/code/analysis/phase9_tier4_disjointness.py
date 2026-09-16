# -*- coding: utf-8 -*-
"""
Tier 4: Empirical disjointness between coding and regulatory selection signals
==============================================================================
Question (associate researcher critique #4/#1 hardening): are the coding-driven
and regulatory-driven evidence axes empirically independent, or does one
dominate everywhere?

Analyses:
  (a) BUSTED coding-selection rate within regulatory-element-proximate subsets
      (HAR-proximate n=476, hCONDEL n=183) vs the genome-wide rate (34.0%).
      If regulatory-proximate genes were ALSO preferentially coding-selected,
      the two axes would be non-independent (one process explaining both).
  (b) GDS-RDS Spearman rank correlation across all 4,974 genes.
  (c) Complementary: RD-classification rate within HAR/hCONDEL subsets
      (regulatory selection concentrating where expected) and GD rate within
      the same subsets.

Output: results/phase9_hardening/tier4_empirical_disjointness.json
        results/phase9_hardening/tier4_empirical_disjointness.csv (per-subset table)
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

BASE = 'd:/人类正选择基因项目/'
OUT_DIR = BASE + 'results/phase9_hardening/'
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------- load
df = pd.read_csv(BASE + 'results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv')
hc = pd.read_csv(BASE + 'results/phase8_tissue_analysis/hcondel_gene_mapping.csv')

n = len(df)
genome_busted = df['bh_fdr_sig'].sum()
genome_rate = genome_busted / n

# hCONDEL gene set (by gene_id, fallback gene_symbol); drop NaN/empty to avoid mass-matching
hcondel_ids = set(hc['gene_id'].dropna().astype(str)) - {''}
hcondel_syms = set(hc['gene_symbol'].dropna().astype(str)) - {''}
sym_valid = df['gene_symbol'].notna() & (df['gene_symbol'].astype(str).str.strip() != '')
df['is_hcondel'] = df['gene_id'].isin(hcondel_ids) | (sym_valid & df['gene_symbol'].isin(hcondel_syms))
print(f'hCONDEL mapping rows: {len(hc)}, unique ids: {len(hcondel_ids)}, matched in universe: {df["is_hcondel"].sum()}')

# subsets
subsets = {
    'HAR-proximate (n_hars>0)': df['n_hars'] > 0,
    'hCONDEL-proximate': df['is_hcondel'],
    'HAR or hCONDEL': (df['n_hars'] > 0) | df['is_hcondel'],
}

results = {'universe_n': n, 'genome_busted_sig': int(genome_busted),
           'genome_busted_rate': round(genome_rate, 4)}

# ---------------------------------------------------- (a) BUSTED rate in subsets
rows = []
for name, mask in subsets.items():
    sub = df[mask]
    k = len(sub)
    sig = sub['bh_fdr_sig'].sum()
    rate = sig / k
    # Fisher exact: subset sig/unsig vs genome rest sig/unsig
    rest = df[~mask]
    table = [[int(sig), int(k - sig)],
             [int(rest['bh_fdr_sig'].sum()), int((~rest['bh_fdr_sig']).sum())]]
    orr, p = stats.fisher_exact(table)
    # RD class rate in subset
    rd = sub['classification_v7'].isin(['regulation-driven', 'dual-driven']).sum()
    gd = sub['classification_v7'].isin(['gene-driven', 'gene-driven (relaxed)']).sum()
    # RD enrichment Fisher (subset vs rest)
    rd_all = df['classification_v7'].isin(['regulation-driven', 'dual-driven'])
    rd_rest = int(rd_all.sum()) - int(rd)
    tab_rd = [[int(rd), int(k - rd)],
              [rd_rest, int(n - k - rd_rest)]]
    orr_rd, p_rd = stats.fisher_exact(tab_rd)
    rows.append({
        'subset': name, 'n': k,
        'busted_sig': int(sig), 'busted_rate': round(rate, 4),
        'busted_OR': round(orr, 3), 'busted_p': p,
        'RD_or_dual': int(rd), 'RD_rate': round(rd / k, 4),
        'RD_OR': round(orr_rd, 3), 'RD_p': p_rd,
        'GD_class': int(gd), 'GD_rate': round(gd / k, 4),
    })
    print(f"{name}: n={k}, BUSTED sig={sig} ({rate:.1%}) OR={orr:.2f} p={p:.2e} | "
          f"RD/dual={rd} ({rd/k:.1%}) OR={orr_rd:.2f} p={p_rd:.2e} | GD={gd} ({gd/k:.1%})")

tab = pd.DataFrame(rows)
tab.to_csv(OUT_DIR + 'tier4_empirical_disjointness.csv', index=False)
results['subset_table'] = rows

# ---------------------------------------------------- (b) GDS-RDS Spearman
g, r = df['gds_v7'], df['rds_v7']
rho, p_rho = stats.spearmanr(g, r)
print(f"\nGDS-RDS Spearman (all {n}): rho={rho:.4f}, p={p_rho:.2e}")
results['gds_rds_spearman_all'] = {'rho': round(rho, 4), 'p': float(p_rho), 'n': n}

# among classified genes only (exclude neutral)
cls = df[df['classification_v7'] != 'neutral']
rho_c, p_c = stats.spearmanr(cls['gds_v7'], cls['rds_v7'])
print(f"GDS-RDS Spearman (classified only n={len(cls)}): rho={rho_c:.4f}, p={p_c:.2e}")
results['gds_rds_spearman_classified'] = {'rho': round(rho_c, 4), 'p': float(p_c), 'n': int(len(cls))}

# component sub-correlations (diagnostic): rds components vs gds
for comp in ['rds_doan', 'rds_tau', 'rds_nc', 'rds_brain']:
    sub = df[[comp, 'gds_v7']].dropna()
    if len(sub) > 10:
        rho_i, p_i = stats.spearmanr(sub['gds_v7'], sub[comp])
        results[f'spearman_gds_vs_{comp}'] = {'rho': round(rho_i, 4), 'p': float(p_i), 'n': int(len(sub))}
        print(f"  GDS vs {comp}: rho={rho_i:+.4f} p={p_i:.2e} (n={len(sub)})")

# mutual information style check: GD class rate among caMPRA+ vs caMPRA-
camp = df['has_doan_campra'].fillna(False)
if camp.sum() > 0:
    gd_mask = df['classification_v7'].isin(['gene-driven', 'gene-driven (relaxed)'])
    t = pd.crosstab(camp, gd_mask)
    orr_c, p_c2 = stats.fisher_exact(t.values)
    print(f"\nGD-class rate among caMPRA+ ({camp.sum()} genes): "
          f"{gd_mask[camp].mean():.1%} vs caMPRA- {gd_mask[~camp].mean():.1%}, OR={orr_c:.3f} p={p_c2:.2e}")
    results['gd_rate_campra_pos'] = {'n_pos': int(camp.sum()),
                                     'gd_rate_pos': round(float(gd_mask[camp].mean()), 4),
                                     'gd_rate_neg': round(float(gd_mask[~camp].mean()), 4),
                                     'OR': round(orr_c, 3), 'p': float(p_c2)}

with open(OUT_DIR + 'tier4_empirical_disjointness.json', 'w', encoding='utf-8') as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
print('\nSaved:', OUT_DIR + 'tier4_empirical_disjointness.json')
