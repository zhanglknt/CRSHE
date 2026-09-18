# -*- coding: utf-8 -*-
"""
R6 review fix (molevol-methods P1-1): Tier 4 reference-frame correction
=======================================================================
Problem A (non-independence): the published rate differences compared the
HAR/hCONDEL subset BUSTED rate with the GENOME-WIDE rate (1,690/4,974 = 34.0%),
but the genome-wide denominator contains the subset itself -> the two
proportions are not independent. Corrected: compare each subset with its
COMPLEMENT (universe minus subset).

Problem B (CI mislabel): tier4_empirical_disjointness.csv columns
rate_diff_ci95_newcombe_low/high actually contained the conservative
subtraction of Wilson bounds (l1 - u2, u1 - l2), not the Newcombe hybrid
score interval. Corrected: true Newcombe (Newcombe 1998, method 10,
square-and-add of Wilson score intervals), computed against the complement;
the Wilson-bound-subtraction interval is retained as a comparison column.

Outputs:
  results/phase9_hardening/tier4_complement_fix.csv / .json
  fixed results/phase9_hardening/tier4_empirical_disjointness.csv / .json
  log: _r6_task1.txt
"""
import io
import json
import math

import pandas as pd
from scipy import stats

BASE = 'd:/人类正选择基因项目/'
OUT = BASE + 'results/phase9_hardening/'
LOG = io.open(BASE + '_r6_task1.txt', 'w', encoding='utf-8')


def P(*a):
    print(*a, file=LOG)


Z = 1.959963984540054


def wilson(k, n):
    p = k / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return c - h, c + h


def newcombe_hybrid(k1, n1, k0, n0):
    """True Newcombe hybrid score interval for p1 - p0 (Newcombe 1998 method 10)."""
    l1, u1 = wilson(k1, n1)
    l0, u0 = wilson(k0, n0)
    p1, p0 = k1 / n1, k0 / n0
    d = p1 - p0
    low = d - math.sqrt((p1 - l1) ** 2 + (u0 - p0) ** 2)
    high = d + math.sqrt((u1 - p1) ** 2 + (p0 - l0) ** 2)
    return low, high


def wilson_bound_subtraction(k1, n1, k0, n0):
    l1, u1 = wilson(k1, n1)
    l0, u0 = wilson(k0, n0)
    return l1 - u0, u1 - l0


def wald(k1, n1, k0, n0):
    p1, p0 = k1 / n1, k0 / n0
    se = math.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    return p1 - p0 - Z * se, p1 - p0 + Z * se


# ------------------------------------------------------------- load universe
df = pd.read_csv(BASE + 'results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv')
hc = pd.read_csv(BASE + 'results/phase8_tissue_analysis/hcondel_gene_mapping.csv')
n = len(df)
genome_sig = int(df['bh_fdr_sig'].sum())

hcondel_ids = set(hc['gene_id'].dropna().astype(str)) - {''}
hcondel_syms = set(hc['gene_symbol'].dropna().astype(str)) - {''}
sym_valid = df['gene_symbol'].notna() & (df['gene_symbol'].astype(str).str.strip() != '')
df['is_hcondel'] = df['gene_id'].isin(hcondel_ids) | (sym_valid & df['gene_symbol'].isin(hcondel_syms))

subsets = {
    'HAR-proximate (n_hars>0)': df['n_hars'] > 0,
    'hCONDEL-proximate': df['is_hcondel'],
    'HAR or hCONDEL': (df['n_hars'] > 0) | df['is_hcondel'],
}

# old (wrong) values currently in the CSV, for the before/after table
old = pd.read_csv(OUT + 'tier4_empirical_disjointness.csv').set_index('subset')

rows = []
before_after = []
for name, mask in subsets.items():
    sub = df[mask]
    comp = df[~mask]
    k1, n1 = int(sub['bh_fdr_sig'].sum()), len(sub)
    k0, n0 = int(comp['bh_fdr_sig'].sum()), len(comp)
    r1, r0 = k1 / n1, k0 / n0
    d = r1 - r0
    ln, hn = newcombe_hybrid(k1, n1, k0, n0)
    lw, hw = wilson_bound_subtraction(k1, n1, k0, n0)
    wl, wh = wald(k1, n1, k0, n0)
    orr, p = stats.fisher_exact([[k1, n1 - k1], [k0, n0 - k0]])
    # RD / GD context on the same subset (unchanged metrics, complement-based Fisher already)
    rd = int(sub['classification_v7'].isin(['regulation-driven', 'dual-driven']).sum())
    gd = int(sub['classification_v7'].isin(['gene-driven', 'gene-driven (relaxed)']).sum())
    rd_all = int(df['classification_v7'].isin(['regulation-driven', 'dual-driven']).sum())
    orr_rd, p_rd = stats.fisher_exact([[rd, n1 - rd], [rd_all - rd, n - n1 - rd_all + rd]])
    o = old.loc[name]
    rows.append({
        'subset': name, 'n': n1, 'busted_sig': k1,
        'subset_rate': round(r1, 4),
        'complement_n': n0, 'complement_sig': k0, 'complement_rate': round(r0, 4),
        'rate_diff': round(d, 4),
        'newcombe_low': round(ln, 4), 'newcombe_high': round(hn, 4),
        'wilson_bound_subtraction_low': round(lw, 4), 'wilson_bound_subtraction_high': round(hw, 4),
        'wald_low': round(wl, 4), 'wald_high': round(wh, 4),
        'fisher_p_vs_complement': p,
        'RD_or_dual': rd, 'RD_rate': round(rd / n1, 4),
        'RD_OR': round(orr_rd, 3), 'RD_p': p_rd,
        'GD_class': gd, 'GD_rate': round(gd / n1, 4),
    })
    before_after.append({
        'subset': name,
        'old_rate_diff_vs_genome': float(o['rate_diff']),
        'old_ci_labeled_newcombe_low': float(o['rate_diff_ci95_newcombe_low']),
        'old_ci_labeled_newcombe_high': float(o['rate_diff_ci95_newcombe_high']),
        'note_old_ci_is_actually_wilson_bound_subtraction': True,
        'old_reference': f'genome-wide {genome_sig}/{n} = {genome_sig/n:.4f} (contains subset)',
        'new_rate_diff_vs_complement': round(d, 4),
        'new_comcombe_low': round(ln, 4), 'new_newcombe_high': round(hn, 4),
        'new_wilson_bound_subtraction_low': round(lw, 4),
        'new_wilson_bound_subtraction_high': round(hw, 4),
        'new_reference': f'complement {k0}/{n0} = {r0:.4f} (disjoint)',
        'fisher_p_vs_complement': p, 'fisher_OR_vs_complement': round(orr, 3),
    })
    P(f"{name}: subset {k1}/{n1} = {100*r1:.2f}% vs complement {k0}/{n0} = {100*r0:.2f}% "
      f"-> diff {100*d:+.2f} pp, true Newcombe 95% CI ({100*ln:+.2f}, {100*hn:+.2f}), "
      f"Wald ({100*wl:+.2f}, {100*wh:+.2f}), Fisher OR={orr:.3f} p={p:.4f}")

fix = pd.DataFrame(rows)
fix.to_csv(OUT + 'tier4_complement_fix.csv', index=False)

result = {
    'task': 'R6 molevol P1-1: Tier 4 subset-vs-genome reference frame correction',
    'problem': ('Previous rate_diff compared subset rate with the genome-wide rate '
                '(1,690/4,974), which includes the subset -> non-independent proportions. '
                'Additionally the "newcombe" CI columns were the Wilson-bound subtraction '
                '(l1-u2, u1-l2), not the Newcombe hybrid score interval.'),
    'universe_n': n, 'genome_busted_sig': genome_sig,
    'method': ('subset vs COMPLEMENT; rate difference with true Newcombe hybrid score '
               '95% CI (Newcombe 1998 method 10, square-and-add of Wilson intervals); '
               'Wald and Wilson-bound-subtraction intervals reported for comparison; '
               'Fisher exact two-sided on the 2x2 subset-vs-complement table.'),
    'subsets': rows,
    'before_after_comparison': before_after,
}
with io.open(OUT + 'tier4_complement_fix.json', 'w', encoding='utf-8') as f:
    json.dump(result, f, indent=2, ensure_ascii=False)

# ------------------------------------------------- fix tier4_empirical_disjointness
csv_path = OUT + 'tier4_empirical_disjointness.csv'
tab = pd.read_csv(csv_path)
keep = ['subset', 'n', 'busted_sig', 'busted_rate', 'busted_OR', 'busted_p',
        'RD_or_dual', 'RD_rate', 'RD_OR', 'RD_p', 'GD_class', 'GD_rate']
new_tab = tab[keep].copy()
# explicit mapping
new_tab['complement_n'] = [r['complement_n'] for r in rows]
new_tab['complement_sig'] = [r['complement_sig'] for r in rows]
new_tab['complement_rate'] = [r['complement_rate'] for r in rows]
new_tab['rate_diff'] = [r['rate_diff'] for r in rows]
new_tab['rate_diff_ci95_newcombe_low'] = [r['newcombe_low'] for r in rows]
new_tab['rate_diff_ci95_newcombe_high'] = [r['newcombe_high'] for r in rows]
new_tab['rate_diff_ci95_wilson_bound_subtraction_low'] = [r['wilson_bound_subtraction_low'] for r in rows]
new_tab['rate_diff_ci95_wilson_bound_subtraction_high'] = [r['wilson_bound_subtraction_high'] for r in rows]
new_tab['rate_diff_ci95_wald_low'] = [r['wald_low'] for r in rows]
new_tab['rate_diff_ci95_wald_high'] = [r['wald_high'] for r in rows]
new_tab['fisher_p_vs_complement'] = [r['fisher_p_vs_complement'] for r in rows]
new_tab.to_csv(csv_path, index=False)
P('\nrewrote tier4_empirical_disjointness.csv (complement reference, true Newcombe)')

# sync JSON: replace the genome-reference block with the complement-reference block
jpath = OUT + 'tier4_empirical_disjointness.json'
with io.open(jpath, encoding='utf-8') as f:
    j = json.load(f)
j.pop('rate_difference_vs_genome', None)
j['rate_difference_vs_complement'] = {
    'method': ('Rate difference of BUSTED FDR-significant proportion, subset vs COMPLEMENT '
               '(rest of universe, disjoint). True Newcombe hybrid score 95% CI (Newcombe 1998 '
               'method 10: square-and-add of Wilson score intervals), reported as '
               'rate_diff_ci95_newcombe_low/high. The previously reported "newcombe" columns '
               'were in fact the conservative subtraction of Wilson bounds (l1-u2, u1-l2); those '
               'values are retained as rate_diff_ci95_wilson_bound_subtraction_low/high. '
               'Wald shown as well. R6 fix: reference changed from genome-wide (which contains '
               'the subset, violating proportion independence) to the complement.'),
    'subsets': {r['subset']: {
        'n_subset': r['n'], 'busted_sig': r['busted_sig'], 'subset_rate': r['subset_rate'],
        'complement_n': r['complement_n'], 'complement_sig': r['complement_sig'],
        'complement_rate': r['complement_rate'], 'rate_diff': r['rate_diff'],
        'rate_diff_ci95_newcombe_low': r['newcombe_low'],
        'rate_diff_ci95_newcombe_high': r['newcombe_high'],
        'rate_diff_ci95_wilson_bound_subtraction_low': r['wilson_bound_subtraction_low'],
        'rate_diff_ci95_wilson_bound_subtraction_high': r['wilson_bound_subtraction_high'],
        'rate_diff_ci95_wald_low': r['wald_low'], 'rate_diff_ci95_wald_high': r['wald_high'],
        'fisher_OR_vs_complement': round(
            stats.fisher_exact([[r['busted_sig'], r['n'] - r['busted_sig']],
                                [r['complement_sig'], r['complement_n'] - r['complement_sig']]])[0], 3),
        'fisher_p_vs_complement': r['fisher_p_vs_complement'],
    } for r in rows},
}
with io.open(jpath, 'w', encoding='utf-8') as f:
    json.dump(j, f, indent=2, ensure_ascii=False)
P('synced tier4_empirical_disjointness.json')
LOG.close()
print('done')
