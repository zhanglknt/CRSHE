# -*- coding: utf-8 -*-
"""Build ENSG->gene_symbol map from GENCODE v47 GTF, save to data/ for reuse."""
import gzip
import re

import pandas as pd

GTF = 'd:/人类正选择基因项目/data/downloads/gtex_v11/gencode.v47.basic.annotation.gtf.gz'
OUT = 'd:/人类正选择基因项目/data/gene_id_to_symbol_gencode47.csv'

ensg_re = re.compile(r'gene_id "([^"]+)"')
sym_re = re.compile(r'gene_name "([^"]+)"')

rows = []
with gzip.open(GTF, 'rt', encoding='utf-8') as f:
    for line in f:
        if line.startswith('#'):
            continue
        if '\tgene\t' not in line:
            continue
        m1 = ensg_re.search(line)
        m2 = sym_re.search(line)
        if m1 and m2:
            gid = m1.group(1).split('.')[0]  # strip version
            rows.append((gid, m2.group(1)))

mp = pd.DataFrame(rows, columns=['gene_id', 'gene_symbol']).drop_duplicates('gene_id')
mp.to_csv(OUT, index=False)

# how many of the 929 unnamed can we recover?
df = pd.read_csv('d:/人类正选择基因项目/results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv')
sym_ok = df['gene_symbol'].notna() & (df['gene_symbol'].astype(str).str.strip() != '')
unnamed = df[~sym_ok]
merged = unnamed[['gene_id']].merge(mp, on='gene_id', how='left')
recovered = merged['gene_symbol'].notna() & (merged['gene_symbol'].astype(str).str.strip() != '') & (merged['gene_symbol'] != 'nan')
with open('d:/人类正选择基因项目/tmp_map.txt', 'w', encoding='utf-8') as f:
    f.write(f'GTF genes: {len(mp)}\nunnamed in v7: {len(unnamed)}\n'
            f'recovered: {recovered.sum()}\nstill missing: {(~recovered).sum()}')
