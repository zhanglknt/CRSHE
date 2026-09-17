# -*- coding: utf-8 -*-
"""
P0-A Step 5: aggregate summary of LOO Tier1/2 re-runs + GWAS negative controls
==============================================================================
Output: results/phase9_hardening/loo_p0a_summary.json
"""
import json

BASE = 'd:/人类正选择基因项目/'
OUT = BASE + 'results/phase9_hardening/'

with open(OUT + 'loo_variant_check.json', encoding='utf-8') as f:
    check = json.load(f)
with open(OUT + 'loo_tier1_summary.json', encoding='utf-8') as f:
    t1 = json.load(f)
with open(OUT + 'tier2_negative_controls.json', encoding='utf-8') as f:
    neg = json.load(f)

# ---- Tier 2 per-variant (32-test family) RD row extraction ----
t2 = {}
for variant in ['LOO-brain', 'LOO-tau']:
    with open(OUT + f'loo_tier2_{variant}.json', encoding='utf-8') as f:
        d = json.load(f)
    rd = [r for r in d['results'] if r['class'] == 'regulation_driven']
    t2[variant] = {
        'n_tests': d['n_tests'],
        'RD_sig_q05': sum(1 for r in rd if r['q_bh'] < 0.05),
        'RD_total': len(rd),
        'RD_by_trait': {r['trait']: {'OR': r['OR'], 'p': r['p'], 'q': r['q_bh']}
                        for r in rd},
    }

# v7 (full) baseline for Tier2, 32-test family, from original tier2 results
import pandas as pd
t2_v7 = pd.read_csv(BASE + 'results/phase9_hardening/tier2_gwas_enrichment.csv')
rd_v7 = t2_v7[t2_v7['class'] == 'regulation_driven']
t2_baseline = {
    'n_tests': int(len(t2_v7)),
    'RD_sig_q05': int((rd_v7['q_bh'] < 0.05).sum()),
    'RD_total': int(len(rd_v7)),
    'RD_by_trait': {r['trait']: {'OR': r['OR'], 'p': r['p'], 'q': r['q_bh']}
                    for _, r in rd_v7.iterrows()},
}

# ---- negative-control headline ----
neg_head = {}
for variant, d in neg['by_classification'].items():
    rd_ctrl = [r for r in d['results']
               if r['class'] == 'regulation_driven' and r['is_control']]
    neg_head[variant] = {
        'RD_neuro_sig_q05': d['RD_neuro_sig_q05'],
        'RD_neuro_total': d['RD_neuro_total'],
        'RD_control_sig_q05': d['RD_control_sig_q05'],
        'RD_control_total': d['RD_control_total'],
        'RD_control_detail': {r['trait']: {'OR': r['OR'], 'p': r['p'], 'q': r['q_bh']}
                              for r in rd_ctrl},
        'control_trait_gene_counts': {t: neg['trait_gene_counts'][t]
                                      for t in ['Height', 'LDL cholesterol',
                                                'Crohn disease', 'Type 2 diabetes']},
    }

# ---- Tier 1 headline ----
t1_head = {}
for variant in ['LOO-brain', 'LOO-tau']:
    t1_head[variant] = {}
    for cls in ['gene_driven', 'gene_driven_all', 'regulation_driven', 'dual_driven']:
        t1_head[variant][cls] = {l: {'terms_tested': t1[variant][cls][l]['terms_tested'],
                                      'bh_005': t1[variant][cls][l]['bh_005']}
                                  for l in ['GO_BP', 'SynGO']}

summary = {
    'task': 'P0-A: Tier 1/2 re-runs under LOO classification variants + GWAS negative controls',
    'step1_loo_classification': {
        'all_match_loo_full_matrix': check['loo_full_matrix_all_match'],
        'per_gene_agreement_with_saved_columns': check['agreement_with_saved_columns'],
        'counts': check['counts'],
        'loo_weights': check['loo_weights'],
        'loo_weights_note': check['loo_weights_source'],
    },
    'step2_tier1_loo': {
        'v7_baseline': {'gene_driven': {'GO_BP_bh005': 0, 'SynGO_bh005': 0},
                        'regulation_driven': {'GO_BP_bh005': 18, 'SynGO_bh005': 11}},
        'by_variant': t1_head,
        'key_findings': (
            'RD retains SynGO FDR-significant enrichment under both LOO variants '
            '(LOO-brain: SynGO 7 terms, GO BP 15 terms; LOO-tau: SynGO 5 terms, '
            'GO BP 1 term). GD retains 0 FDR-significant terms in every library '
            'under every variant (v7: 0/0; LOO-brain: 0/0; LOO-tau: 0/0).'),
    },
    'step3_tier2_loo': {
        'v7_baseline_32test_family': t2_baseline,
        'by_variant_32test_family': t2,
        'key_findings': (
            'RD enrichment for neuro/cognitive GWAS traits is robust but attenuated '
            'under LOO: v7 8/8 traits FDR-sig; LOO-brain 6/8 (Intelligence q=0.052, '
            'ASD q=0.154); LOO-tau 5/8 (Intelligence q=0.094, ASD q=0.098, '
            'Bipolar q=0.507). GD shows no enrichment in any variant.'),
    },
    'step4_negative_controls': {
        'control_traits': neg['control_traits'],
        'gwas_scan': neg['gwas_scan'],
        'trait_gene_counts': neg['trait_gene_counts'],
        'by_classification': neg_head,
        'key_findings': (
            'Negative-control specificity is PARTIAL: RD is enriched for Crohn '
            '(OR<1, null) and LDL cholesterol (null at v7/LOO-brain, borderline '
            'LOO-tau q=0.052), but IS significantly enriched for Height '
            '(OR 2.2-2.3, all classifications) and Type 2 diabetes (OR 2.1-2.3, '
            'all classifications). RD neuro enrichment persists after joint BH '
            'with controls (v7 8/8; LOO-brain 6/8; LOO-tau 5/8), but RD is NOT '
            'specific to neuro traits among highly polygenic, regulatory-architectured '
            'traits. Height/T2D enrichment likely reflects the regulatory/noncoding '
            'architecture of RD genes rather than cognition-specific biology.'),
    },
    'output_files': {
        'step1': ['loo_variant_per_gene.csv', 'loo_variant_check.json'],
        'step2': ['loo_tier1_summary.json',
                  'loo_tier1_{variant}_{class}_{lib}.csv (16 files)'],
        'step3': ['loo_tier2_LOO-brain.csv/.json', 'loo_tier2_LOO-tau.csv/.json'],
        'step4': ['tier2_negative_controls.csv/.json'],
    },
    'scripts': ['scripts/phase9_loo_tier12_rerun.py (step 1)',
                'scripts/phase9_loo_tier1.py (step 2)',
                'scripts/phase9_loo_tier2_negctl.py (steps 3+4)',
                'scripts/phase9_loo_p0a_summary.py (step 5)'],
}

with open(OUT + 'loo_p0a_summary.json', 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2, ensure_ascii=False, default=float)
print('saved:', OUT + 'loo_p0a_summary.json')
for k in ['step1_loo_classification', 'step2_tier1_loo', 'step3_tier2_loo',
          'step4_negative_controls']:
    print('\n==', k, '==')
    v = summary[k]
    for kk, vv in v.items():
        if kk != 'key_findings' and not isinstance(vv, str):
            print(kk, ':', json.dumps(vv, default=float)[:300])
    print('key_findings:', v.get('key_findings', ''))
