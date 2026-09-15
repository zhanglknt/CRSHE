# -*- coding: utf-8 -*-
"""
Phase 7A: Parse PAML v4 candidate gene results
Extract positive selection signals from branch-site model outputs.

For each gene:
  - H0 (null): foreground omega fixed at 1.0
  - H1 (alternative): foreground omega > 1.0 allowed
  - LRT = 2 * (lnL_H1 - lnL_H0), df = 1
  - BEB sites with Prob(w>1) > 0.95

Output: CSV with all genes + positively selected gene list
"""

import os
import re
import json
import logging
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
from datetime import datetime

# ===================== Configuration =====================
PROJECT_ROOT = r"d:\人类正选择基因项目"
PAML_DIR = os.path.join(PROJECT_ROOT, "results", "phase7_gene_vs_regulation", 
                         "phase7a_paml_results_v4_candidates")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "results", "phase7_gene_vs_regulation",
                          "phase7a_paml_parsed")
HAR_OVERLAP = os.path.join(PROJECT_ROOT, "results", "phase7_gene_vs_regulation",
                           "phase7e_har_overlap", "har_gene_overlap.csv")
SHARED_GENES = os.path.join(PROJECT_ROOT, "results", "phase2_final", "shared_genes_all10.csv")

# PAML significance thresholds
PVALUE_THRESHOLD = 0.05
BEB_THRESHOLD = 0.95
OMEGA_THRESHOLD = 1.0

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ===================== Logging =====================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)


def parse_lnl(filepath):
    """Parse log-likelihood and number of parameters from PAML output."""
    if not os.path.exists(filepath):
        return None, None
    
    lnl = None
    np_count = None
    
    with open(filepath, 'r', errors='replace') as f:
        for line in f:
            # lnL(ntime: 18  np: 23): -27832.058334      +0.000000
            m = re.search(r'lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*([\-\d.]+)', line)
            if m:
                np_count = int(m.group(1))
                lnl = float(m.group(2))
                break
    
    return lnl, np_count


def parse_omega_and_proportions(filepath):
    """Parse foreground omega, proportions, and BEB sites from H1 output."""
    if not os.path.exists(filepath):
        return None
    
    result = {
        'omega_fg': None,
        'p0': None, 'p1': None, 'p2a': None, 'p2b': None,
        'bg_w0': None, 'bg_w1': None, 'bg_w2a': None, 'bg_w2b': None,
        'fg_w0': None, 'fg_w1': None, 'fg_w2a': None, 'fg_w2b': None,
        'kappa': None,
        'beb_sites': [],
        'beb_sites_high': [],
    }
    
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines):
        # kappa (ts/tv) = 0.43506
        m = re.search(r'kappa\s*\(ts/tv\)\s*=\s*([\d.]+)', line)
        if m:
            result['kappa'] = float(m.group(1))
        
        # MLEs of dN/dS (w) for site classes (K=4)
        # site class             0        1       2a       2b
        # proportion       0.10506  0.89494  0.00000  0.00000
        # background w     0.08658  1.00000  0.08658  1.00000
        # foreground w     0.08658  1.00000 11.72060 11.72060
        if 'proportion' in line and 'site class' not in line:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    result['p0'] = float(parts[1])
                    result['p1'] = float(parts[2])
                    result['p2a'] = float(parts[3])
                    result['p2b'] = float(parts[4])
                except (ValueError, IndexError):
                    pass
        
        if 'background w' in line:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    result['bg_w0'] = float(parts[2])
                    result['bg_w1'] = float(parts[3])
                    result['bg_w2a'] = float(parts[4])
                    result['bg_w2b'] = float(parts[5])
                except (ValueError, IndexError):
                    pass
        
        if 'foreground w' in line:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    result['fg_w0'] = float(parts[2])
                    result['fg_w1'] = float(parts[3])
                    result['fg_w2a'] = float(parts[4])
                    result['fg_w2b'] = float(parts[5])
                    # foreground omega for positive selection class
                    result['omega_fg'] = max(float(parts[4]), float(parts[5]))
                except (ValueError, IndexError):
                    pass
        
        # BEB analysis
        # Positive sites for foreground lineages Prob(w>1):
        # Followed by lines like: 1 K 0.987*  5 R 0.963*
        if 'Positive sites for foreground lineages' in line:
            # Read next lines for BEB sites
            j = i + 1
            while j < len(lines) and j < i + 50:
                beb_line = lines[j].strip()
                if not beb_line or beb_line.startswith('The grid') or beb_line.startswith('Posterior'):
                    break
                # Parse BEB sites: "1 K 0.987*" or "13 A 0.956*  42 G 0.973*"
                sites = re.findall(r'(\d+)\s+([A-Z\*])\s+([\d.]+)\*?', beb_line)
                for site_num, aa, prob in sites:
                    prob_val = float(prob)
                    result['beb_sites'].append({
                        'site': int(site_num),
                        'aa': aa,
                        'prob': prob_val
                    })
                    if prob_val >= BEB_THRESHOLD:
                        result['beb_sites_high'].append({
                            'site': int(site_num),
                            'aa': aa,
                            'prob': prob_val
                        })
                j += 1
    
    return result


def parse_gene_name(filename):
    """Extract gene ID from filename like ENSG00000000971_H1.out"""
    return filename.replace('_H1.out', '').replace('_H0.out', '')


def load_har_overlap():
    """Load HAR overlap data."""
    if not os.path.exists(HAR_OVERLAP):
        logger.warning(f"HAR overlap file not found: {HAR_OVERLAP}")
        return {}
    
    df = pd.read_csv(HAR_OVERLAP)
    har_dict = {}
    for _, row in df.iterrows():
        gene_id = str(row.get('gene_id', ''))
        n_hars = int(row.get('n_hars', 0))
        har_dict[gene_id] = {
            'har_count': n_hars,
            'has_har': n_hars > 0,
            'gene_name': str(row.get('gene_name', '')),
        }
    return har_dict


def load_gene_symbols():
    """Load gene ID to symbol mapping from shared genes file and HAR overlap."""
    gene_map = {}
    
    # From HAR overlap CSV (has gene_name column)
    if os.path.exists(HAR_OVERLAP):
        df_har = pd.read_csv(HAR_OVERLAP)
        for _, row in df_har.iterrows():
            gid = str(row.get('gene_id', ''))
            gname = str(row.get('gene_name', ''))
            if gid and gname and gname != 'nan':
                gene_map[gid] = gname
    
    # From shared genes CSV (has human_gene_name column)
    if os.path.exists(SHARED_GENES):
        df_shared = pd.read_csv(SHARED_GENES, encoding='utf-8-sig')
        id_col = None
        sym_col = None
        for col in df_shared.columns:
            if col.lower() == 'human_gene_id':
                id_col = col
            if col.lower() == 'human_gene_name':
                sym_col = col
        if id_col and sym_col:
            for _, row in df_shared.iterrows():
                gid = str(row[id_col])
                gname = str(row[sym_col])
                if gid and gname and gname != 'nan' and gid not in gene_map:
                    gene_map[gid] = gname
    
    return gene_map


def main():
    logger.info("=" * 70)
    logger.info("Phase 7A: Parse PAML v4 Candidate Gene Results")
    logger.info("=" * 70)
    
    # Find all gene IDs
    h0_files = [f for f in os.listdir(PAML_DIR) if f.endswith('_H0.out')]
    h1_files = [f for f in os.listdir(PAML_DIR) if f.endswith('_H1.out')]
    
    gene_ids = set()
    for f in h0_files:
        gene_ids.add(parse_gene_name(f))
    for f in h1_files:
        gene_ids.add(parse_gene_name(f))
    
    logger.info(f"Found {len(gene_ids)} unique genes")
    logger.info(f"H0 files: {len(h0_files)}, H1 files: {len(h1_files)}")
    
    # Load HAR overlap and gene symbols
    har_data = load_har_overlap()
    gene_symbols = load_gene_symbols()
    logger.info(f"Loaded HAR overlap for {len(har_data)} genes")
    logger.info(f"Loaded gene symbols for {len(gene_symbols)} genes")
    
    # Parse each gene
    results = []
    parse_errors = 0
    
    for i, gene_id in enumerate(sorted(gene_ids)):
        if (i + 1) % 100 == 0:
            logger.info(f"  Processing gene {i+1}/{len(gene_ids)}...")
        
        h0_path = os.path.join(PAML_DIR, f"{gene_id}_H0.out")
        h1_path = os.path.join(PAML_DIR, f"{gene_id}_H1.out")
        
        # Parse lnL
        lnl_h0, np_h0 = parse_lnl(h0_path)
        lnl_h1, np_h1 = parse_lnl(h1_path)
        
        if lnl_h0 is None or lnl_h1 is None:
            parse_errors += 1
            results.append({
                'gene_id': gene_id,
                'gene_symbol': gene_symbols.get(gene_id, ''),
                'status': 'parse_failed',
            })
            continue
        
        # Parse omega and BEB from H1
        h1_data = parse_omega_and_proportions(h1_path)
        
        if h1_data is None:
            parse_errors += 1
            results.append({
                'gene_id': gene_id,
                'gene_symbol': gene_symbols.get(gene_id, ''),
                'status': 'h1_parse_failed',
            })
            continue
        
        # Calculate LRT
        lrt = 2 * (lnl_h1 - lnl_h0)
        df = max(1, (np_h1 or 23) - (np_h0 or 22))
        
        # Guard against negative LRT (numerical issues)
        if lrt < 0:
            lrt = 0
            p_value = 1.0
        else:
            p_value = stats.chi2.sf(lrt, df)
        
        # Extract key values
        omega_fg = h1_data.get('omega_fg')
        p2a = h1_data.get('p2a', 0) or 0
        p2b = h1_data.get('p2b', 0) or 0
        p_positive = p2a + p2b
        
        beb_sites = h1_data.get('beb_sites', [])
        beb_sites_high = h1_data.get('beb_sites_high', [])
        
        # Format BEB sites as string
        beb_str = '; '.join([f"{s['site']}{s['aa']}({s['prob']:.3f})" for s in beb_sites])
        beb_high_str = '; '.join([f"{s['site']}{s['aa']}({s['prob']:.3f})" for s in beb_sites_high])
        
        # HAR overlap
        har_info = har_data.get(gene_id, {})
        
        # Determine positive selection status
        is_lrt_significant = p_value < PVALUE_THRESHOLD
        has_omega_gt_1 = omega_fg is not None and omega_fg > OMEGA_THRESHOLD
        has_beb_sites = len(beb_sites_high) > 0
        has_positive_proportion = p_positive > 0.001
        
        # Classification
        if is_lrt_significant and has_omega_gt_1 and has_beb_sites:
            ps_status = "positive_selection"
        elif is_lrt_significant and has_omega_gt_1 and has_positive_proportion:
            ps_status = "likely_positive"
        elif is_lrt_significant:
            ps_status = "lrt_significant_no_beb"
        else:
            ps_status = "neutral"
        
        results.append({
            'gene_id': gene_id,
            'gene_symbol': gene_symbols.get(gene_id, ''),
            'status': 'ok',
            'lnl_h0': lnl_h0,
            'lnl_h1': lnl_h1,
            'lrt_statistic': lrt,
            'lrt_df': df,
            'p_value': p_value,
            'omega_foreground': omega_fg,
            'p0': p2a,  # proportion of positive selection sites
            'p_positive': p_positive,
            'kappa': h1_data.get('kappa'),
            'n_beb_sites': len(beb_sites),
            'n_beb_sites_high': len(beb_sites_high),
            'beb_sites': beb_str,
            'beb_sites_high': beb_high_str,
            'has_har_overlap': har_info.get('has_har', False),
            'har_count': har_info.get('har_count', 0),
            'ps_status': ps_status,
        })
    
    logger.info(f"Parsing complete. Errors: {parse_errors}")
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Multiple testing correction (Benjamini-Hochberg FDR)
    valid_mask = df['status'] == 'ok'
    if valid_mask.sum() > 0:
        valid = df[valid_mask].copy()
        # Sort p-values for BH correction
        n = len(valid)
        valid = valid.sort_values('p_value')
        valid['rank'] = range(1, n + 1)
        valid['fdr_bh'] = valid['p_value'] * n / valid['rank']
        # Enforce monotonicity (from the end)
        valid['fdr_bh'] = valid['fdr_bh'].iloc[::-1].cummin().iloc[::-1]
        valid['fdr_bh'] = valid['fdr_bh'].clip(upper=1.0)
        valid = valid.drop(columns=['rank'])
        df = pd.concat([valid, df[~valid_mask]], ignore_index=True)
    
    # Save full results
    full_csv = os.path.join(OUTPUT_DIR, "paml_results_all_genes.csv")
    df.to_csv(full_csv, index=False)
    logger.info(f"Full results saved: {full_csv} ({len(df)} genes)")
    
    # Filter positively selected genes
    ps_genes = df[df['ps_status'] == 'positive_selection'].sort_values('p_value')
    ps_csv = os.path.join(OUTPUT_DIR, "positively_selected_genes.csv")
    ps_genes.to_csv(ps_csv, index=False)
    logger.info(f"Positive selection genes: {len(ps_genes)} -> {ps_csv}")
    
    # Likely positive
    likely_ps = df[df['ps_status'].isin(['positive_selection', 'likely_positive'])].sort_values('p_value')
    likely_csv = os.path.join(OUTPUT_DIR, "likely_positive_genes.csv")
    likely_ps.to_csv(likely_csv, index=False)
    logger.info(f"Likely positive genes (incl. no BEB): {len(likely_ps)} -> {likely_csv}")
    
    # LRT significant (broader)
    lrt_sig = df[(df['p_value'] < PVALUE_THRESHOLD) & (valid_mask)].sort_values('p_value')
    lrt_csv = os.path.join(OUTPUT_DIR, "lrt_significant_genes.csv")
    lrt_sig.to_csv(lrt_csv, index=False)
    logger.info(f"LRT significant genes (p<{PVALUE_THRESHOLD}): {len(lrt_sig)} -> {lrt_csv}")
    
    # FDR significant
    fdr_sig = df[(df['fdr_bh'] < 0.10) & (valid_mask)].sort_values('p_value')
    fdr_csv = os.path.join(OUTPUT_DIR, "fdr_significant_genes.csv")
    fdr_sig.to_csv(fdr_csv, index=False)
    logger.info(f"FDR significant genes (FDR<0.10): {len(fdr_sig)} -> {fdr_csv}")
    
    # Summary statistics
    summary = {
        'total_genes': len(df),
        'parsed_ok': int(valid_mask.sum()),
        'parse_errors': parse_errors,
        'lrt_significant_p005': int(lrt_sig.shape[0]),
        'fdr_significant_010': int(fdr_sig.shape[0]),
        'omega_gt_1': int(df[(df['omega_foreground'] > 1) & valid_mask].shape[0]),
        'has_beb_sites': int(df[(df['n_beb_sites_high'] > 0) & valid_mask].shape[0]),
        'positive_selection': int(ps_genes.shape[0]),
        'likely_positive': int(likely_ps.shape[0]),
        'has_har_overlap': int(df[df['has_har_overlap'] == True].shape[0]),
        'ps_with_har': int(ps_genes[ps_genes['has_har_overlap'] == True].shape[0]) if len(ps_genes) > 0 else 0,
        'timestamp': datetime.now().isoformat(),
    }
    
    summary_path = os.path.join(OUTPUT_DIR, "parsing_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info("=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    for k, v in summary.items():
        logger.info(f"  {k}: {v}")
    
    # Print top positively selected genes
    if len(ps_genes) > 0:
        logger.info("\nTop positively selected genes:")
        for _, row in ps_genes.head(20).iterrows():
            logger.info(f"  {row['gene_id']} ({row['gene_symbol']}): "
                       f"p={row['p_value']:.2e}, omega={row['omega_foreground']:.2f}, "
                       f"BEB={row['n_beb_sites_high']}, HAR={row['har_count']}")
    
    return df


if __name__ == "__main__":
    df = main()
