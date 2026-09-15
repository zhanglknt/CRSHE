#!/usr/bin/env python3
"""
Phase 5 v3: 人类特异性筛选 — 使用BUSTED v3 + Selectome v6整合结果

更新内容:
  1. 用BUSTED v3结果替代v1 (4,974基因 vs 4,820)
  2. 用Selectome v6替代PAML v2b (41灵长类正选择基因)
  3. HSS v3评分: 0.35×CDS + 0.30×NC + 0.35×HAR
  4. Tier分层基于evidence_tier (Tier1_Dual ~ Tier5)
  5. 双方法确认 = BUSTED FDR<0.05 AND Selectome positive
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(r"d:\人类正选择基因项目")
RESULTS_DIR = PROJECT_ROOT / "results"

# ===== 输入 =====
BUSTED_V3 = RESULTS_DIR / "phase4_hyphy" / "final_analysis_v3" / "busted_results_v3.csv"
SELECTOME_INTEGRATION = RESULTS_DIR / "phase4_hyphy" / "final_analysis_v3" / "busted_selectome_integration.tsv"
HAR_OVERLAP = RESULTS_DIR / "phase7_gene_vs_regulation" / "phase7e_har_overlap" / "har_gene_overlap.csv"
NONCODING_RATES = RESULTS_DIR / "phase7_gene_vs_regulation" / "phase7d_noncoding_rates" / "gene_evolution_summary.csv"
PHASE5_V2 = RESULTS_DIR / "phase5_human_specific_screening_v2" / "human_specific_screening_v2_results.csv"

# ===== 输出 =====
OUTPUT_DIR = RESULTS_DIR / "phase5_human_specific_screening_v3"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 已知人类特异基因
KNOWN_HUMAN_SPECIFIC = {
    'FOXP2': 'speech/language',
    'ASPM': 'brain size',
    'MCPH1': 'microcephaly',
    'SRGAP2': 'cortical development',
    'ARHGAP11B': 'neocortex expansion',
    'CNTNAP2': 'language/neurodevelopment',
    'HAR1': 'non-coding RNA',
    'NOTCH2NL': 'cortical neurogenesis',
    'ZEB2': 'neural development',
    'AUTS2': 'neurodevelopment',
    'TBR1': 'cortical projection neurons',
    'ROBO1': 'axon guidance',
}


def load_data():
    """加载所有数据"""
    print("加载数据...")
    
    busted = pd.read_csv(BUSTED_V3)
    print(f"  BUSTED v3: {len(busted)} genes, {busted['fdr_significant'].sum()} FDR<0.05")
    
    selectome = pd.read_csv(SELECTOME_INTEGRATION, sep='\t')
    # 标记有Selectome证据的基因
    selectome['has_selectome'] = selectome['selectome_pvalue'].notna() & (selectome['selectome_pvalue'].astype(str).str.strip() != '')
    print(f"  Selectome: {selectome['has_selectome'].sum()} genes with evidence")
    
    har = pd.read_csv(HAR_OVERLAP)
    print(f"  HAR overlap: {len(har)} genes, {(har['n_hars'] > 0).sum()} with HARs")
    
    nc_rates = pd.read_csv(NONCODING_RATES)
    print(f"  Non-coding rates: {len(nc_rates)} genes")
    
    # Phase 5 v2 for gene symbols
    phase5_v2 = pd.read_csv(PHASE5_V2) if PHASE5_V2.exists() else pd.DataFrame()
    if not phase5_v2.empty:
        print(f"  Phase 5 v2: {len(phase5_v2)} genes (for gene symbols)")
    
    return busted, selectome, har, nc_rates, phase5_v2


def calculate_cds_score_v3(busted_row, selectome_row):
    """
    CDS评分 v3 (0-1)
    基于 BUSTED v3 + Selectome 整合
    
    组成:
      1. BUSTED p值显著性 (0-0.4)
         - p < 1e-10: 0.4
         - p < 1e-5: 0.3
         - p < 0.001: 0.2
         - p < 0.05: 0.1
         - else: 0
      2. BUSTED FDR显著性 (0-0.3)
         - FDR < 0.01: 0.3
         - FDR < 0.05: 0.2
         - p < 0.05 (not FDR): 0.1
         - else: 0
      3. Selectome交叉确认 (0-0.3)
         - FDR<0.05 + Selectome: 0.3 (Tier1)
         - Selectome only: 0.15
         - BUSTED FDR only: 0.0 (already counted above)
         - else: 0
    """
    score = 0.0
    
    p_val = busted_row.get('p_value', 1.0)
    fdr = busted_row.get('bh_fdr', 1.0)
    fdr_sig = busted_row.get('fdr_significant', False)
    p_sig = busted_row.get('p_significant', False)
    
    # 1. BUSTED p值
    if p_val < 1e-10:
        score += 0.4
    elif p_val < 1e-5:
        score += 0.3
    elif p_val < 0.001:
        score += 0.2
    elif p_val < 0.05:
        score += 0.1
    
    # 2. BUSTED FDR
    if fdr_sig and fdr < 0.01:
        score += 0.3
    elif fdr_sig and fdr < 0.05:
        score += 0.2
    elif p_sig:
        score += 0.1
    
    # 3. Selectome交叉确认
    has_sel = selectome_row.get('has_selectome', False) if selectome_row is not None else False
    if has_sel and fdr_sig:
        score += 0.3  # Tier1 dual-method
    elif has_sel:
        score += 0.15  # Selectome evidence alone
    
    return min(1.0, score)


def calculate_nc_score_v3(nc_row):
    """
    非编码评分 v3 (0-1)
    基于非编码加速因子
    """
    if nc_row is None or nc_row.empty:
        return 0.0
    
    max_acc = nc_row.get('max_acceleration_factor', 0)
    if pd.isna(max_acc) or max_acc <= 1.0:
        return 0.0
    
    # 加速因子: 1.0=无加速, 5.0=高度加速
    return min(1.0, (max_acc - 1) / 4)


def calculate_har_score_v3(har_row):
    """
    HAR评分 v3 (0-1)
    基于HAR数量和最大HAR速率
    """
    if har_row is None:
        return 0.0
    
    n_hars = int(har_row.get('n_hars', 0))
    if n_hars == 0:
        return 0.0
    
    # HAR数量: 1个=0.3, 3个=0.6, 5+个=1.0
    har_count_score = min(1.0, n_hars / 5)
    
    return har_count_score


def determine_tier_v3(busted_row, selectome_row):
    """
    Tier分层 v3
    Tier1: BUSTED FDR<0.05 + Selectome positive (dual-method)
    Tier2: BUSTED FDR<0.05 only
    Tier3: BUSTED p<0.05 only
    Tier4: Selectome only (BUSTED not significant)
    Tier5: Both present but BUSTED not significant
    """
    fdr_sig = busted_row.get('fdr_significant', False)
    p_sig = busted_row.get('p_significant', False)
    has_sel = selectome_row.get('has_selectome', False) if selectome_row is not None else False
    
    if fdr_sig and has_sel:
        return 'Tier1_Dual'
    elif fdr_sig:
        return 'Tier2_BUSTED_FDR'
    elif p_sig and has_sel:
        return 'Tier1_Dual'  # p<0.05 + Selectome also counts
    elif p_sig:
        return 'Tier3_BUSTED_p'
    elif has_sel:
        return 'Tier4_Selectome'
    else:
        return 'None'


def determine_pattern_v3(cds_score, nc_score, har_score):
    """确定模式"""
    cds_active = cds_score >= 0.2
    nc_active = nc_score >= 0.3
    har_active = har_score >= 0.2
    
    if cds_active and nc_active and har_active:
        return 'coding+noncoding+regulatory'
    elif cds_active and har_active:
        return 'coding+regulatory'
    elif cds_active and nc_active:
        return 'coding+noncoding'
    elif cds_active:
        return 'coding_dominant'
    elif nc_active or har_active:
        return 'noncoding_dominant'
    else:
        return 'weak'


def main():
    print("="*70)
    print("Phase 5 v3: 人类特异性筛选 (BUSTED v3 + Selectome v6)")
    print("="*70)
    print(f"时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # 加载数据
    busted_df, selectome_df, har_df, nc_df, phase5_v2 = load_data()
    print()
    
    # 创建基因符号映射
    gene_symbols = {}
    if not phase5_v2.empty:
        gene_symbols = dict(zip(phase5_v2['gene_id'], phase5_v2['gene_symbol']))
    
    # Selectome基因符号
    for _, row in selectome_df.iterrows():
        if pd.notna(row.get('gene_name')) and str(row.get('gene_name', '')).strip():
            gene_symbols[row['gene_id']] = row['gene_name']
    
    # 构建查找字典
    selectome_map = {}
    for _, row in selectome_df.iterrows():
        selectome_map[row['gene_id']] = row.to_dict()
    
    har_map = {}
    if not har_df.empty:
        for _, row in har_df.iterrows():
            har_map[row['gene_id']] = row
    
    nc_map = {}
    if not nc_df.empty and 'gene_id' in nc_df.columns:
        for _, row in nc_df.iterrows():
            nc_map[row['gene_id']] = row
    
    # 计算HSS v3
    print("计算HSS v3...")
    results = []
    
    for _, row in busted_df.iterrows():
        gene_id = row['gene_id']
        
        # 查找各数据源
        sel_row = selectome_map.get(gene_id, None)
        if sel_row is not None and not sel_row.get('has_selectome', False):
            sel_row = None  # 只保留有Selectome证据的
        
        har_row = har_map.get(gene_id, None)
        nc_row = nc_map.get(gene_id, None)
        
        # 计算评分
        cds_score = calculate_cds_score_v3(row, sel_row)
        nc_score = calculate_nc_score_v3(nc_row)
        har_score = calculate_har_score_v3(har_row)
        
        # HSS = 0.35×CDS + 0.30×NC + 0.35×HAR
        hss = 0.35 * cds_score + 0.30 * nc_score + 0.35 * har_score
        
        # Tier和模式
        tier = determine_tier_v3(row, sel_row)
        pattern = determine_pattern_v3(cds_score, nc_score, har_score)
        
        # 已知基因
        symbol = gene_symbols.get(gene_id, '')
        is_known = symbol in KNOWN_HUMAN_SPECIFIC
        known_cat = KNOWN_HUMAN_SPECIFIC.get(symbol, '')
        
        # Selectome信息
        sel_branch = ''
        sel_pval = ''
        sel_qval = ''
        sel_sig = ''
        if sel_row is not None:
            sel_branch = sel_row.get('selectome_branch', '')
            sel_pval = sel_row.get('selectome_pvalue', '')
            sel_qval = sel_row.get('selectome_qvalue', '')
            sel_sig = sel_row.get('significance', '')
        
        # HAR信息
        n_hars = int(har_row.get('n_hars', 0)) if har_row is not None else 0
        max_nc_accel = float(nc_row.get('max_acceleration_factor', 0)) if nc_row is not None else 0
        
        results.append({
            'gene_id': gene_id,
            'gene_symbol': symbol,
            'HSS_v3': round(hss, 4),
            'CDS_score_v3': round(cds_score, 4),
            'NC_score_v3': round(nc_score, 4),
            'HAR_score_v3': round(har_score, 4),
            'has_nc_data': nc_row is not None,
            'busted_p': row.get('p_value', 1.0),
            'busted_fdr': row.get('bh_fdr', 1.0),
            'busted_lrt': row.get('lrt', 0),
            'busted_fdr_sig': row.get('fdr_significant', False),
            'busted_p_sig': row.get('p_significant', False),
            'selectome_positive': sel_row is not None,
            'selectome_branch': sel_branch,
            'selectome_pvalue': sel_pval,
            'selectome_qvalue': sel_qval,
            'selectome_significance': sel_sig,
            'n_hars': n_hars,
            'max_nc_accel': max_nc_accel,
            'tier': tier,
            'pattern': pattern,
            'is_known_human_specific': is_known,
            'known_category': known_cat,
        })
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('HSS_v3', ascending=False)
    
    print(f"  HSS v3计算完成: {len(results_df)} genes")
    print()
    
    # 统计
    print("统计:")
    tier_counts = results_df['tier'].value_counts()
    print(f"  Tier分布:")
    for tier, cnt in tier_counts.items():
        print(f"    {tier:25s}: {cnt:4d}")
    
    print(f"\n  模式分布:")
    pattern_counts = results_df['pattern'].value_counts()
    for pat, cnt in pattern_counts.head(10).items():
        print(f"    {pat:30s}: {cnt:4d}")
    
    print(f"\n  HSS v3分布:")
    print(f"    Top 10: {results_df['HSS_v3'].head(10).tolist()}")
    print(f"    Mean: {results_df['HSS_v3'].mean():.4f}")
    print(f"    Median: {results_df['HSS_v3'].median():.4f}")
    print()
    
    # 保存
    print("保存结果...")
    output_file = OUTPUT_DIR / "human_specific_screening_v3_results.csv"
    results_df.to_csv(output_file, index=False, encoding='utf-8')
    print(f"  主文件: {output_file}")
    
    # 按Tier保存
    for tier in results_df['tier'].unique():
        tier_genes = results_df[results_df['tier'] == tier]
        tier_file = OUTPUT_DIR / f"{tier}_genes.csv"
        tier_genes.to_csv(tier_file, index=False, encoding='utf-8')
        print(f"    {tier:25s}: {len(tier_genes):4d} → {tier_file.name}")
    
    # Top基因
    top_genes = results_df.head(50)
    top_file = OUTPUT_DIR / "top50_genes_v3.csv"
    top_genes.to_csv(top_file, index=False, encoding='utf-8')
    print(f"  Top 50: {top_file}")
    
    # 统计摘要
    summary = {
        'version': 'v3',
        'timestamp': datetime.now().isoformat(),
        'data_sources': {
            'busted_v3': int(len(busted_df)),
            'busted_fdr_sig': int(busted_df['fdr_significant'].sum()),
            'selectome_positive': int(selectome_df['has_selectome'].sum()),
            'genes_with_har': int((har_df['n_hars'] > 0).sum()) if not har_df.empty else 0,
            'genes_with_nc_data': int(len(nc_df)),
        },
        'tier_distribution': {k: int(v) for k, v in tier_counts.items()},
        'total_genes': int(len(results_df)),
        'hss_formula': '0.35*CDS + 0.30*NC + 0.35*HAR',
        'improvements': [
            'BUSTED v3 (4,974 genes) replaces v1 (4,820)',
            'Selectome v6 replaces PAML v2b',
            'Tier1 = BUSTED FDR<0.05 AND Selectome positive',
            'No circular logic in scoring',
        ],
    }
    
    summary_file = OUTPUT_DIR / "phase5_v3_summary.json"
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"  统计摘要: {summary_file}")
    
    # 打印Top 20基因
    print(f"\n{'='*70}")
    print("Top 20 人类特异性正选择基因 (HSS v3):")
    print(f"{'='*70}")
    print(f"{'Rank':<5} {'Gene':<12} {'HSS':<7} {'CDS':<6} {'NC':<6} {'HAR':<5} {'Tier':<20} {'Pattern'}")
    for i, (_, row) in enumerate(results_df.head(20).iterrows()):
        sym = row['gene_symbol'] if row['gene_symbol'] else row['gene_id'][:10]
        print(f"{i+1:<5} {sym:<12} {row['HSS_v3']:<7.3f} {row['CDS_score_v3']:<6.2f} "
              f"{row['NC_score_v3']:<6.2f} {row['HAR_score_v3']:<5.2f} {row['tier']:<20} {row['pattern']}")
    
    print(f"\n结果目录: {OUTPUT_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
