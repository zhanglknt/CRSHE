# -*- coding: utf-8 -*-
"""
R6 review (GWAS side) confound calibration for Paper A Tier 2 GWAS anchoring
==============================================================================
P0-3  Gene length / position-mapping confounding
      (a) Spearman: gene span (TSS-TES, GENCODE v47) vs RDS total + 4 components
      (b) Logistic per trait: RD status ~ trait-gene status + log10(gene length), 12 traits
      (c) Length-tertile-stratified Fisher (+ CMH) for EA/SCZ/intelligence/height/T2D
      (d) Verdict on neuro vs height/T2D fate after length adjustment
P0-4  Neuro family vs negative-control family formal test
      (a) Per-trait log OR + SE, inverse-variance pooled per family, between-family z
      (b) Trait-label permutation (1000x)
      (c) Consistency (LOO sig proportions) Fisher tests
P1-1  MAPPED_GENE sensitivity: (a) locus-level (1Mb cluster, lead SNP, 1 gene) and
          (b) author-REPORTED GENE(S) only, for EA/SCZ/intelligence/height
P1-2  8 neuro trait set Jaccard matrix + Li & Ji (2005) effective # tests
P1-3  12-trait power table (K, overlap, OR, 95% CI, BH q)
P1-4  Master table: 12 traits x {full, LOO-brain, LOO-tau} x {RD, GD, GD-all}

Outputs (results/phase9_hardening/):
  gwas_confound_calibration.json          (master summary)
  gwas_p03_logistic_adjusted.csv
  gwas_p03_tertile_fisher.csv
  gwas_p04_family_permutation.csv / .json
  gwas_p11_mapping_sensitivity.csv
  gwas_p12_jaccard_matrix.csv
  gwas_p13_power_table.csv
  tier2_anchoring_master_table.csv
"""
import gzip
import json
import re
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

BASE = 'd:/人类正选择基因项目/'
GWAS = BASE + 'data/gwas/gwas-catalog-download-associations-v1.0-full.tsv'
GTF = BASE + 'data/downloads/gtex_v11/gencode.v47.basic.annotation.gtf.gz'
OUT = BASE + 'results/phase9_hardening/'
RNG = np.random.default_rng(20260917)

NEURO = ['Schizophrenia', 'ASD', 'Educational attainment', 'Intelligence',
         'Major depression', 'Bipolar disorder', 'Neuroticism', 'Cognitive performance']
CONTROLS = ['Height', 'LDL cholesterol', 'Crohn disease', 'Type 2 diabetes']
ALL_TRAITS = NEURO + CONTROLS
trait_keys = {
    'Schizophrenia': ('schizophren',), 'ASD': ('autis',),
    'Educational attainment': ('educational attainment',), 'Intelligence': ('intelligence',),
    'Major depression': ('major depress',), 'Bipolar disorder': ('bipolar',),
    'Neuroticism': ('neuroticism',), 'Cognitive performance': ('cognitive performance',),
    'Height': ('height',),
    'LDL cholesterol': ('ldl cholesterol', 'low density lipoprotein cholesterol'),
    'Crohn disease': ('crohn',), 'Type 2 diabetes': ('type 2 diabetes',),
}
token_re = re.compile(r'[A-Za-z][A-Za-z0-9@.\-]*')

# ------------------------------------------------------------------ universe
pg = pd.read_csv(OUT + 'loo_variant_per_gene.csv')
# replicate the exact symbol construction of phase9_loo_tier2_negctl.py:
# loo-CSV gene_symbol, falling back to the GENCODE v47 symbol map
gmap = pd.read_csv(BASE + 'data/gene_id_to_symbol_gencode47.csv')
pg = pg.merge(gmap, on='gene_id', how='left', suffixes=('', '_g'))
pg['sym_final'] = pg['gene_symbol'].where(
    pg['gene_symbol'].notna() & (pg['gene_symbol'].astype(str).str.strip() != ''),
    pg['gene_symbol_g'])
pg['sym'] = pg['sym_final'].astype(str).str.upper().str.strip()
uni_syms = set(pg['sym'])
N_sym = len(uni_syms)

cls_v = pd.read_csv(BASE + 'results/phase7_gene_vs_regulation/phase7g_classification_v7/'
                    'gene_classification_v7.csv')
keep = ['gene_id', 'classification_v7', 'rds_v7', 'rds_doan', 'rds_tau', 'rds_nc', 'rds_brain']
pg = pg.merge(cls_v[keep], on='gene_id', how='left', suffixes=('', '_v7'))
assert pg['rds_v7_v7'].notna().all() if 'rds_v7_v7' in pg else True

# ------------------------------------------------------------------ gene spans (GTF)
gid_re = re.compile(r'gene_id "([^"]+)"')
span = {}
with gzip.open(GTF, 'rt', encoding='utf-8') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        m = gid_re.search(parts[8])
        if not m:
            continue
        gid = m.group(1).split('.')[0]
        span[gid] = int(parts[4]) - int(parts[3]) + 1
pg['gene_len'] = pg['gene_id'].map(span)
n_missing_len = int(pg['gene_len'].isna().sum())
pg['log10_len'] = np.log10(pg['gene_len'])

# ------------------------------------------------------------------ GWAS scan
# collect per GWS row: traits matched, chr, pos, p, mapped universe tokens, reported tokens
trait_genes = {t: set() for t in ALL_TRAITS}
trait_genes_locus = {t: set() for t in ALL_TRAITS}   # P1-1a convention
trait_genes_reported = {t: set() for t in ALL_TRAITS}  # P1-1b convention
rows_by_trait = {t: [] for t in ALL_TRAITS}  # (pos, p, mapped_tokens) for clustering
n_rows = n_gws = 0


def mapped_universe_tokens(field):
    toks = []
    for tok in token_re.findall(field):
        t = tok.upper().strip('.-')
        if t in uni_syms:
            toks.append(t)
    return toks


with open(GWAS, encoding='utf-8', errors='replace') as f:
    header = f.readline().rstrip('\n').split('\t')
    ix = {c: i for i, c in enumerate(header)}
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < len(header):
            continue
        n_rows += 1
        trait = parts[ix['DISEASE/TRAIT']].lower()
        matched = [t for t in ALL_TRAITS if any(k in trait for k in trait_keys[t])]
        if not matched:
            continue
        try:
            p = float(parts[ix['P-VALUE']])
        except (ValueError, IndexError):
            continue
        if not p < 5e-8:
            continue
        n_gws += 1
        mtoks = mapped_universe_tokens(parts[ix['MAPPED_GENE']])
        rtoks = mapped_universe_tokens(parts[ix['REPORTED GENE(S)']])
        try:
            pos = int(parts[ix['CHR_POS']])
        except (ValueError, TypeError):
            pos = None
        try:
            chrom = parts[ix['CHR_ID']].strip()
        except (ValueError, TypeError):
            chrom = None
        for t in matched:
            trait_genes[t].update(mtoks)
            trait_genes_reported[t].update(rtoks)
            if pos is not None:
                rows_by_trait[t].append((chrom, pos, p, mtoks))


def cluster_loci(recs, gap=1_000_000):
    """Greedy interval merge of lead-SNP positions within +/-1Mb per chromosome.
    Returns list of clusters; each cluster = list of records."""
    by_chrom = {}
    for r in recs:
        by_chrom.setdefault(r[0], []).append(r)
    clusters = []
    for chrom, rs in by_chrom.items():
        rs.sort(key=lambda x: x[1])
        cur = [rs[0]]
        for r in rs[1:]:
            if r[1] - cur[-1][1] <= gap:
                cur.append(r)
            else:
                clusters.append(cur)
                cur = [r]
        clusters.append(cur)
    return clusters


locus_stats = {}
for t in ALL_TRAITS:
    clusters = cluster_loci(rows_by_trait[t])
    locus_stats[t] = {'n_loci': len(clusters)}
    for cl in clusters:
        lead = min(cl, key=lambda x: x[2])  # lowest p-value row
        # convention (a): count ONE gene per independent locus = first universe
        # token listed in the lead SNP's MAPPED_GENE
        if lead[3]:
            trait_genes_locus[t].add(lead[3][0])

# ------------------------------------------------------------------ helpers
def bh_q(pvals):
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    return out


def log_or_se(a, b, c, d):
    """2x2 table [[k, n_cls-k], [K-k, rest]] -> (log OR, SE, OR, lo, hi).
    Haldane-Anscombe 0.5 correction when any cell is 0."""
    if min(a, b, c, d) == 0:
        a2, b2, c2, d2 = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    else:
        a2, b2, c2, d2 = a, b, c, d
    lor = np.log(a2 * d2 / (b2 * c2))
    se = np.sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
    return lor, se, np.exp(lor), np.exp(lor - 1.96 * se), np.exp(lor + 1.96 * se)


CLS_MASKS = {
    'regulation_driven': lambda s: s == 'regulation-driven',
    'gene_driven': lambda s: s == 'gene-driven',
    'gene_driven_all': lambda s: s.isin(['gene-driven', 'gene-driven (relaxed)']),
    'dual_driven': lambda s: s == 'dual-driven',
}


def fisher_rd(trait_set, cls_col, cname='regulation_driven', one_sided=False):
    """Set-based Fisher exactly like phase9_tier2 script."""
    mask = CLS_MASKS[cname](pg[cls_col])
    cls_syms = set(pg.loc[mask, 'sym'])
    n_cls = len(cls_syms)
    K = len(trait_set)
    k = len(cls_syms & trait_set)
    table = [[k, n_cls - k], [K - k, (N_sym - n_cls) - (K - k)]]
    alt = 'greater' if one_sided else 'two-sided'
    orr, p = stats.fisher_exact(table, alternative=alt)
    return {'K_universe': K, 'n_class': n_cls, 'k_overlap': k, 'OR': orr, 'p': p,
            'a': k, 'b': n_cls - k, 'c': K - k,
            'd': (N_sym - n_cls) - (K - k)}


# ================================================================== P0-3(a)
p03a = {}
sub = pg.dropna(subset=['gene_len'])
for col in ['rds_v7', 'rds_doan', 'rds_tau', 'rds_nc', 'rds_brain']:
    rho, pv = stats.spearmanr(sub['gene_len'], sub[col], nan_policy='omit')
    p03a[col] = {'spearman_rho': float(rho), 'p': float(pv), 'n': int(len(sub))}
p03a['gene_len_summary'] = {
    'n_with_span': int(len(sub)), 'n_missing_span': n_missing_len,
    'median': float(sub['gene_len'].median()), 'mean': float(sub['gene_len'].mean())}

# also: RD / GD genes vs others gene length (Wilcoxon + median ratio)
rd_mask = pg['cls_full'] == 'regulation-driven'
gd_mask = pg['cls_full'].isin(['gene-driven', 'gene-driven (relaxed)'])
len_rd = pg.loc[rd_mask, 'gene_len'].dropna()
len_rest = pg.loc[~(rd_mask | gd_mask), 'gene_len'].dropna()
u, pv = stats.mannwhitneyu(len_rd, len_rest, alternative='two-sided')
p03a['RD_vs_rest_len'] = {'median_RD': float(len_rd.median()), 'median_rest': float(len_rest.median()),
                          'ratio': float(len_rd.median() / len_rest.median()),
                          'wilcoxon_p': float(pv)}
len_gd = pg.loc[gd_mask, 'gene_len'].dropna()
u2, pv2 = stats.mannwhitneyu(len_gd, len_rest, alternative='two-sided')
p03a['GD_vs_rest_len'] = {'median_GD': float(len_gd.median()), 'median_rest': float(len_rest.median()),
                          'ratio': float(len_gd.median() / len_rest.median()),
                          'wilcoxon_p': float(pv2)}
# do trait gene sets contain longer genes?
trait_len = {}
for t in ALL_TRAITS:
    in_t = pg['sym'].isin(trait_genes[t])
    lt_in = pg.loc[in_t, 'gene_len'].dropna()
    lt_out = pg.loc[~in_t, 'gene_len'].dropna()
    if len(lt_in) > 0:
        uu, pp = stats.mannwhitneyu(lt_in, lt_out, alternative='two-sided')
        trait_len[t] = {'median_in': float(lt_in.median()), 'median_out': float(lt_out.median()),
                        'ratio': float(lt_in.median() / lt_out.median()), 'p': float(pp)}
p03a['trait_gene_set_length'] = trait_len

# ================================================================== P0-3(b) logistic
p03b_rows = []
pg_reg = pg.dropna(subset=['gene_len']).copy()
pg_reg['RD'] = (pg_reg['cls_full'] == 'regulation-driven').astype(int)
for t in ALL_TRAITS:
    pg_reg['trait'] = pg_reg['sym'].isin(trait_genes[t]).astype(int)
    # unadjusted (trait only) and adjusted (trait + log10 len)
    X0 = sm.add_constant(pg_reg[['trait']].astype(float))
    m0 = sm.Logit(pg_reg['RD'], X0).fit(disp=0)
    X1 = sm.add_constant(pg_reg[['trait', 'log10_len']].astype(float))
    m1 = sm.Logit(pg_reg['RD'], X1).fit(disp=0)
    b_t0 = m0.params['trait']; p_t0 = m0.pvalues['trait']
    b_t1 = m1.params['trait']; p_t1 = m1.pvalues['trait']
    b_l = m1.params['log10_len']; p_l = m1.pvalues['log10_len']
    or_t = np.exp(b_t1); lo = np.exp(b_t1 - 1.96 * m1.bse['trait']); hi = np.exp(b_t1 + 1.96 * m1.bse['trait'])
    p03b_rows.append({'trait': t, 'is_control': t in CONTROLS,
                      'OR_unadj': float(np.exp(b_t0)), 'p_unadj': float(p_t0),
                      'OR_adj': float(or_t), 'CI_lo': float(lo), 'CI_hi': float(hi),
                      'p_adj': float(p_t1),
                      'len_logOR': float(b_l), 'len_p': float(p_l),
                      'n_trait_genes_in_reg': int(pg_reg['trait'].sum())})
p03b = pd.DataFrame(p03b_rows)
p03b['q_adj_bh'] = bh_q(p03b['p_adj'].values)
p03b.to_csv(OUT + 'gwas_p03_logistic_adjusted.csv', index=False)

# ================================================================== P0-3(c) tertile Fisher + CMH
p03c_rows = []
pg3 = pg.dropna(subset=['gene_len']).copy()
q1, q2 = pg3['gene_len'].quantile([1 / 3, 2 / 3])
pg3['tertile'] = np.digitize(pg3['gene_len'], [q1, q2])  # 0,1,2
tertile_results = {}
for t in ['Educational attainment', 'Schizophrenia', 'Intelligence', 'Height', 'Type 2 diabetes']:
    pg3['trait'] = pg3['sym'].isin(trait_genes[t])
    tables = []
    for ter in [0, 1, 2]:
        s = pg3[pg3['tertile'] == ter]
        tab = pd.crosstab(s['trait'], s['cls_full'] == 'regulation-driven').values
        # arrange [[trait&RD, trait&~RD],[~trait&RD, ~trait&~RD]]
        a = int(((s['trait']) & (s['cls_full'] == 'regulation-driven')).sum())
        b = int((s['trait'] & (s['cls_full'] != 'regulation-driven')).sum())
        c_ = int((~s['trait'] & (s['cls_full'] == 'regulation-driven')).sum())
        d_ = int((~s['trait'] & (s['cls_full'] != 'regulation-driven')).sum())
        orr, p = stats.fisher_exact([[a, b], [c_, d_]], alternative='two-sided')
        p03c_rows.append({'trait': t, 'tertile': ter + 1,
                          'n_genes': int(len(s)), 'a_trait_RD': a, 'b_trait_nonRD': b,
                          'c_nontrait_RD': c_, 'OR': orr, 'p': p})
        tables.append(np.array([[a, b], [c_, d_]]))
    # CMH pooled
    try:
        from statsmodels.stats.contingency_tables import StratifiedTable
        st = StratifiedTable(tables)
        cmh = st.test_null_odds()
        cmh_or = float(st.oddsratio_pooled)
        try:
            cmh_lo, cmh_hi = st.oddsratio_pooled_confint()
        except Exception:
            cmh_lo, cmh_hi = np.nan, np.nan
        tertile_results[t] = {'cmh_common_OR': cmh_or,
                              'cmh_CI': [float(cmh_lo), float(cmh_hi)],
                              'cmh_p': float(cmh.pvalue), 'cmh_stat': float(cmh.statistic)}
    except Exception as e:
        tertile_results[t] = {'error': str(e)}
    p03c_rows[-3:][0]  # noop
p03c = pd.DataFrame(p03c_rows)
p03c.to_csv(OUT + 'gwas_p03_tertile_fisher.csv', index=False)

# ================================================================== P0-4(a) family pooled logOR
p04a_rows = []
for t in ALL_TRAITS:
    fr = fisher_rd(trait_genes[t], 'cls_full', one_sided=False)
    lor, se, orr, lo, hi = log_or_se(fr['a'], fr['b'], fr['c'], fr['d'])
    p04a_rows.append({'trait': t, 'is_control': t in CONTROLS, 'K': fr['K_universe'],
                      'k': fr['k_overlap'], 'logOR': lor, 'SE': se, 'OR': orr,
                      'CI_lo': lo, 'CI_hi': hi, 'w': 1 / se ** 2,
                      'fisher_p_two_sided': fr['p']})
p04a = pd.DataFrame(p04a_rows)


def pooled(df_):
    w = df_['w'].values
    beta = np.sum(w * df_['logOR'].values) / np.sum(w)
    se_b = np.sqrt(1 / np.sum(w))
    return beta, se_b


b_neuro, se_neuro = pooled(p04a[~p04a['is_control']])
b_ctrl, se_ctrl = pooled(p04a[p04a['is_control']])
z_between = (b_neuro - b_ctrl) / np.sqrt(se_neuro ** 2 + se_ctrl ** 2)
p_between = 2 * stats.norm.sf(abs(z_between))
p04a_summary = {
    'neuro_pooled_logOR': b_neuro, 'neuro_pooled_SE': se_neuro,
    'neuro_pooled_OR': np.exp(b_neuro),
    'neuro_CI': [np.exp(b_neuro - 1.96 * se_neuro), np.exp(b_neuro + 1.96 * se_neuro)],
    'control_pooled_logOR': b_ctrl, 'control_pooled_SE': se_ctrl,
    'control_pooled_OR': np.exp(b_ctrl),
    'control_CI': [np.exp(b_ctrl - 1.96 * se_ctrl), np.exp(b_ctrl + 1.96 * se_ctrl)],
    'z_between': z_between, 'p_between': p_between,
}
p04a_summary['len_adjusted_between'] = {}
# same but with length-adjusted per-trait logORs from P0-3(b)
adj = p03b.set_index('trait')
lor_adj = adj['OR_adj'].apply(lambda x: np.log(x)).values
# SE from CI
se_adj = (np.log(adj['CI_hi'].values) - np.log(adj['CI_lo'].values)) / (2 * 1.96)
is_ctrl = adj['is_control'].values.astype(bool)
w_adj = 1 / se_adj ** 2
b_n = np.sum(w_adj[~is_ctrl] * lor_adj[~is_ctrl]) / np.sum(w_adj[~is_ctrl])
s_n = np.sqrt(1 / np.sum(w_adj[~is_ctrl]))
b_c = np.sum(w_adj[is_ctrl] * lor_adj[is_ctrl]) / np.sum(w_adj[is_ctrl])
s_c = np.sqrt(1 / np.sum(w_adj[is_ctrl]))
z_ab = (b_n - b_c) / np.sqrt(s_n ** 2 + s_c ** 2)
p04a_summary['len_adjusted_between'] = {
    'neuro_pooled_OR': float(np.exp(b_n)), 'control_pooled_OR': float(np.exp(b_c)),
    'z_between': float(z_ab), 'p_between': float(2 * stats.norm.sf(abs(z_ab)))}

# ================================================================== P0-4(b) permutation
obs_diff = b_neuro - b_ctrl
lor_all = p04a['logOR'].values
w_all = p04a['w'].values
n_perm = 1000
diffs = np.empty(n_perm)
for i in range(n_perm):
    idx = RNG.permutation(12)
    neuro_idx = idx[:8]
    ctrl_idx = idx[8:]
    wn = w_all[neuro_idx]; wc = w_all[ctrl_idx]
    bn = np.sum(wn * lor_all[neuro_idx]) / np.sum(wn)
    bc = np.sum(wc * lor_all[ctrl_idx]) / np.sum(wc)
    diffs[i] = bn - bc
p_perm = float((np.sum(np.abs(diffs) >= abs(obs_diff)) + 1) / (n_perm + 1))
p04b = {'n_perm': n_perm, 'obs_diff_logOR': obs_diff,
        'perm_mean_diff': float(diffs.mean()), 'perm_sd': float(diffs.std(ddof=1)),
        'p_two_sided_empirical': p_perm,
        'perm_q025': float(np.quantile(diffs, 0.025)), 'perm_q975': float(np.quantile(diffs, 0.975))}

# ================================================================== P0-4(c) consistency across LOO
# Reproduce the original 48-test convention (12 traits x 4 classes, one-sided
# Fisher, BH within classification) then extract the RD rows.
p04c = {}
cls_cols = {'full': 'cls_full', 'LOO-brain': 'cls_LOO-brain', 'LOO-tau': 'cls_LOO-tau'}
for vname, col in cls_cols.items():
    rows = []
    for t in ALL_TRAITS:
        for cname in ['gene_driven', 'gene_driven_all', 'regulation_driven', 'dual_driven']:
            fr = fisher_rd(trait_genes[t], col, cname, one_sided=True)
            fr['trait'] = t
            fr['is_control'] = t in CONTROLS
            fr['class'] = cname
            rows.append(fr)
    tab = pd.DataFrame(rows)
    tab['q_bh'] = bh_q(tab['p'].values)  # BH across 48 tests, as in prior analysis
    rd = tab[tab['class'] == 'regulation_driven']
    ns = int(((~rd['is_control']) & (rd['q_bh'] < 0.05)).sum())
    cs = int(((rd['is_control']) & (rd['q_bh'] < 0.05)).sum())
    fet_p = stats.fisher_exact([[ns, 8 - ns], [cs, 4 - cs]], alternative='two-sided')[1]
    p04c[vname] = {'RD_neuro_sig': ns, 'RD_neuro_total': 8,
                   'RD_control_sig': cs, 'RD_control_total': 4,
                   'prop_test_fisher_p': fet_p,
                   'which_neuro_sig': rd.loc[(~rd['is_control']) & (rd['q_bh'] < 0.05), 'trait'].tolist(),
                   'which_control_sig': rd.loc[(rd['is_control']) & (rd['q_bh'] < 0.05), 'trait'].tolist()}

p04_out = pd.DataFrame(p04a_rows)
p04_out.to_csv(OUT + 'gwas_p04_family_permutation.csv', index=False)

# ================================================================== P1-1 mapping sensitivity
p11_rows = []
for t in ['Educational attainment', 'Schizophrenia', 'Intelligence', 'Height']:
    for conv, gset in [('original_MAPPED', trait_genes[t]),
                       ('locus_dedup', trait_genes_locus[t]),
                       ('author_reported', trait_genes_reported[t])]:
        fr = fisher_rd(gset, 'cls_full', one_sided=True)
        lor, se, orr, lo, hi = log_or_se(fr['a'], fr['b'], fr['c'], fr['d'])
        p11_rows.append({'trait': t, 'convention': conv, 'n_genes': fr['K_universe'],
                         'k_RD': fr['k_overlap'], 'n_class': fr['n_class'],
                         'OR': orr, 'CI_lo': lo, 'CI_hi': hi, 'p_one_sided': fr['p'],
                         'loci': locus_stats[t]['n_loci'] if conv == 'locus_dedup' else None})
p11 = pd.DataFrame(p11_rows)
p11.to_csv(OUT + 'gwas_p11_mapping_sensitivity.csv', index=False)

# ================================================================== P1-2 Jaccard + Li & Ji
M = np.zeros((8, 8))
for i, ti in enumerate(NEURO):
    for j, tj in enumerate(NEURO):
        M[i, j] = len(trait_genes[ti] & trait_genes[tj]) / len(trait_genes[ti] | trait_genes[tj])
jac = pd.DataFrame(M, index=NEURO, columns=NEURO)
jac.to_csv(OUT + 'gwas_p12_jaccard_matrix.csv')

eig = np.linalg.eigvalsh(M)  # ascending; Jaccard kernel is PSD
eig = np.sort(eig)[::-1]
eig_pos = eig[eig > 0]


def li_ji_meff(eigenvalues):
    """Li & Ji (2005): integer parts contribute directly; fractional parts are
    packed (first-fit decreasing) into whole eigenvalues, each full bin adds 1."""
    ints = np.floor(eigenvalues).astype(int)
    fracs = np.sort(eigenvalues - ints)[::-1]
    bins = []
    for f in fracs:
        placed = False
        for bi in range(len(bins)):
            if bins[bi] + f <= 1.0 + 1e-12:
                bins[bi] += f
                placed = True
                break
        if not placed:
            bins.append(f)
    extra = sum(1 for b in bins if b >= 1.0 - 1e-9)
    return int(ints.sum()) + extra


meff_liji = li_ji_meff(eig_pos)
meff_kaiser = int((eig >= 1).sum())
p12 = {'jaccard_range_offdiag': [float(M[np.triu_indices(8, 1)].min()),
                                  float(M[np.triu_indices(8, 1)].max()),
                                  float(M[np.triu_indices(8, 1)].mean())],
       'eigenvalues': [float(x) for x in eig],
       'Meff_LiJi2005': meff_liji, 'Meff_Kaiser_ge1': meff_kaiser,
       'pairwise_max_jaccard_traitpair': None}
iu = np.triu_indices(8, 1)
mx = np.argmax(M[iu])
p12['pairwise_max_jaccard_traitpair'] = [NEURO[iu[0][mx]], NEURO[iu[1][mx]], float(M[iu][mx])]

# ================================================================== P1-3 power table
p13_rows = []
for t in ALL_TRAITS:
    fr = fisher_rd(trait_genes[t], 'cls_full', one_sided=False)
    lor, se, orr, lo, hi = log_or_se(fr['a'], fr['b'], fr['c'], fr['d'])
    p13_rows.append({'trait': t, 'is_control': t in CONTROLS,
                     'n_universe_mapped_genes': fr['K_universe'],
                     'k_RD_overlap': fr['k_overlap'], 'n_RD': fr['n_class'],
                     'expected': fr['K_universe'] * fr['n_class'] / N_sym,
                     'OR': orr, 'CI_lo': lo, 'CI_hi': hi,
                     'CI_width_log': float(np.log(hi) - np.log(lo)),
                     'fisher_p': fr['p']})
p13 = pd.DataFrame(p13_rows)
p13['q_bh'] = bh_q(p13['fisher_p'].values)
p13.to_csv(OUT + 'gwas_p13_power_table.csv', index=False)

# ================================================================== P1-4 master table
master_rows = []
class_defs = {'RD': ('regulation_driven', None),
              'GD_strict': ('gene_driven_strict', None),
              'GD_all': ('gene_driven_all', None)}
for vname, col in cls_cols.items():
    rows = []
    for t in ALL_TRAITS:
        for cname in ['regulation_driven', 'gene_driven_strict', 'gene_driven_all']:
            if cname == 'regulation_driven':
                mask = pg[col] == 'regulation-driven'
            elif cname == 'gene_driven_strict':
                mask = pg[col] == 'gene-driven'
            else:
                mask = pg[col].isin(['gene-driven', 'gene-driven (relaxed)'])
            cls_syms = set(pg.loc[mask, 'sym'])
            n_cls = len(cls_syms)
            K = len(trait_genes[t])
            k = len(cls_syms & trait_genes[t])
            table = [[k, n_cls - k], [K - k, (N_sym - n_cls) - (K - k)]]
            orr, p = stats.fisher_exact(table, alternative='two-sided')
            lor, se, orr2, lo, hi = log_or_se(k, n_cls - k, K - k, (N_sym - n_cls) - (K - k))
            rows.append({'classification': vname, 'trait': t, 'is_control': t in CONTROLS,
                         'class': cname, 'n_trait_universe': K, 'n_class': n_cls,
                         'k_overlap': k, 'OR': orr2, 'CI_lo': lo, 'CI_hi': hi,
                         'p_two_sided': p})
    tab = pd.DataFrame(rows)
    tab['q_bh'] = bh_q(tab['p_two_sided'].values)  # BH within classification (36 tests)
    master_rows.append(tab)
master = pd.concat(master_rows, ignore_index=True)
master.to_csv(OUT + 'tier2_anchoring_master_table.csv', index=False)

# ------------------------------------------------------------------ sanity check vs prior
prior = pd.read_csv(OUT + 'tier2_negative_controls.csv')
prior_k = prior[prior['classification'] == 'full_v7'].set_index(['trait', 'class'])['k_overlap']
mine_k = p04a.set_index('trait')['k']
sanity = {'n_prior_rows_compared': int(len(prior_k)), 'k_matches': 0, 'mismatches': []}
# compare all 48 (trait x class) cells using the full table from P0-4c reconstruction
_recon = []
for t in ALL_TRAITS:
    for cname in ['gene_driven', 'gene_driven_all', 'regulation_driven', 'dual_driven']:
        fr = fisher_rd(trait_genes[t], 'cls_full', cname, one_sided=True)
        _recon.append({'trait': t, 'class': cname, 'k': fr['k_overlap']})
_recon = pd.DataFrame(_recon).set_index(['trait', 'class'])
for (t, c), row in prior_k.items():
    if (t, c) in _recon.index:
        if int(row) == int(_recon.loc[(t, c), 'k']):
            sanity['k_matches'] += 1
        else:
            sanity['mismatches'].append({'trait': t, 'class': c, 'prior': int(row),
                                        'mine': int(_recon.loc[(t, c), 'k'])})

# ------------------------------------------------------------------ write master JSON
result = {
    'meta': {'script': 'scripts/phase9_gwas_confound.py',
             'universe_symbols': N_sym, 'n_gwas_rows_scanned': n_rows,
             'n_gws_rows_target_traits': n_gws,
             'gene_span_source': 'GENCODE v47 basic annotation (GTEx v11 download)',
             'genes_with_span': int(len(sub)), 'genes_missing_span': n_missing_len},
    'P0_3a_length_vs_RDS': p03a,
    'P0_3b_logistic': p03b.to_dict(orient='records'),
    'P0_3c_tertiles': {'fisher_by_tertile': p03c_rows, 'cmh': tertile_results},
    'P0_4a_family_test': p04a_summary,
    'P0_4b_permutation': p04b,
    'P0_4c_consistency': p04c,
    'P1_1_mapping_sensitivity': p11.to_dict(orient='records'),
    'P1_2_jaccard_meff': p12,
    'P1_3_power_table': p13.to_dict(orient='records'),
    'P1_4_master_table_stats': {
        'n_rows': int(len(master)),
        'RD_sig_q05_by_classification': {
            v: {'neuro': int(((master.classification == v) & (master['class'] == 'regulation_driven')
                              & (~master.is_control) & (master.q_bh < 0.05)).sum()),
                'control': int(((master.classification == v) & (master['class'] == 'regulation_driven')
                                & (master.is_control) & (master.q_bh < 0.05)).sum())}
            for v in cls_cols}},
    'sanity_vs_prior_negative_controls': sanity,
    'locus_stats': locus_stats,
}
with open(OUT + 'gwas_confound_calibration.json', 'w', encoding='utf-8') as f:
    json.dump(result, f, indent=2, ensure_ascii=False, default=float)

# ------------------------------------------------------------------ console report
lines = []
lines.append('=== P0-3a length vs RDS Spearman ===')
for k, v in p03a.items():
    if k in ('rds_v7', 'rds_doan', 'rds_tau', 'rds_nc', 'rds_brain'):
        lines.append(f'  {k}: rho={v["spearman_rho"]:.4f} p={v["p"]:.2e}')
    elif k in ('RD_vs_rest_len', 'GD_vs_rest_len'):
        lines.append(f'  {k}: {v}')
lines.append(f'  span missing: {n_missing_len}')
lines.append('=== P0-3b logistic (RD ~ trait + log10 len) ===')
for _, r in p03b.iterrows():
    lines.append(f'  {r["trait"]:24s} OR_unadj={r["OR_unadj"]:.3f} (p={r["p_unadj"]:.2e}) -> '
                 f'OR_adj={r["OR_adj"]:.3f} (p={r["p_adj"]:.2e}, q={r["q_adj_bh"]:.3f}) len_p={r["len_p"]:.1e}')
lines.append('=== P0-3c CMH (tertile-stratified) ===')
for t, v in tertile_results.items():
    lines.append(f'  {t}: {v}')
lines.append('=== P0-4a family ===')
lines.append(f'  neuro pooled OR={p04a_summary["neuro_pooled_OR"]:.3f} '
             f'CI=({p04a_summary["neuro_CI"][0]:.3f},{p04a_summary["neuro_CI"][1]:.3f})')
lines.append(f'  control pooled OR={p04a_summary["control_pooled_OR"]:.3f} '
             f'CI=({p04a_summary["control_CI"][0]:.3f},{p04a_summary["control_CI"][1]:.3f})')
lines.append(f'  between z={p04a_summary["z_between"]:.3f} p={p04a_summary["p_between"]:.2e}')
lines.append(f'  len-adjusted between: {p04a_summary["len_adjusted_between"]}')
lines.append('=== P0-4b permutation ===')
lines.append(f'  obs diff={p04b["obs_diff_logOR"]:.3f} perm p={p04b["p_two_sided_empirical"]:.3f}')
lines.append('=== P0-4c consistency ===')
for v, d in p04c.items():
    lines.append(f'  {v}: neuro {d["RD_neuro_sig"]}/8 vs control {d["RD_control_sig"]}/4 '
                 f'(Fisher p={d["prop_test_fisher_p"]:.3f}) sig={d["which_neuro_sig"]}')
lines.append('=== P1-1 mapping sensitivity (RD, full v7) ===')
for _, r in p11.iterrows():
    lines.append(f'  {r["trait"]:24s} {r["convention"]:16s} n={r["n_genes"]:4d} k={r["k_RD"]:3d} '
                 f'OR={r["OR"]:.3f} CI=({r["CI_lo"]:.2f},{r["CI_hi"]:.2f}) p={r["p_one_sided"]:.2e}')
lines.append('=== P1-2 Jaccard/Li&Ji ===')
lines.append(f'  off-diag mean={p12["jaccard_range_offdiag"][2]:.3f} max={p12["pairwise_max_jaccard_traitpair"]}')
lines.append(f'  Meff Li&Ji={meff_liji} Kaiser={meff_kaiser}')
lines.append('=== P1-3 power table (RD, full) ===')
for _, r in p13.iterrows():
    lines.append(f'  {r["trait"]:24s} K={r["n_universe_mapped_genes"]:4d} k={r["k_RD_overlap"]:3d} '
                 f'OR={r["OR"]:.3f} CI=({r["CI_lo"]:.2f},{r["CI_hi"]:.2f}) p={r["fisher_p"]:.3f} q={r["q_bh"]:.3f}')
lines.append('=== P1-4 master RD sig ===')
for v, d in result['P1_4_master_table_stats']['RD_sig_q05_by_classification'].items():
    lines.append(f'  {v}: neuro {d["neuro"]}/8 control {d["control"]}/4')
lines.append('=== sanity (k_overlap vs prior full_v7, all 48 trait x class cells) ===')
lines.append(f'  matches: {sanity["k_matches"]}/{sanity["n_prior_rows_compared"]} '
             f'mismatches: {sanity["mismatches"]}')
with open(OUT + 'log_gwas_confound.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print('DONE')
