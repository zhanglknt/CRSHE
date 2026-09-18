# -*- coding: utf-8 -*-
"""Extract longest-span gene coordinates for shared genes from GENCODE GTF
(replicates phase7e_har_overlap.extract_gene_coordinates). Saves to JSON cache."""
import json, csv, sys

GTF = r'd:\人类正选择基因项目\data\human\Homo_sapiens_genes.gtf'
SHARED = r'd:\人类正选择基因项目\results\phase2_final\shared_genes_all10.csv'
OUT = r'd:\人类正选择基因项目\results\phase9_hardening\_gtf_gene_coords.json'

gene_ids = set()
with open(SHARED, encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        gene_ids.add(row['human_gene_id'].split('.')[0])

coords = {}
with open(GTF, encoding='utf-8') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        attrs = parts[8]
        gid = None
        for attr in attrs.split(';'):
            attr = attr.strip()
            if attr.startswith('gene_id'):
                gid = attr.split('"')[1] if '"' in attr else attr.split()[-1].strip('"')
                if '.' in gid:
                    gid = gid.split('.')[0]
                break
        if gid and gid in gene_ids:
            chrom, start, end = parts[0], int(parts[3]), int(parts[4])
            if gid not in coords or end - start > coords[gid][2] - coords[gid][1]:
                coords[gid] = [chrom, start, end]

with open(OUT, 'w', encoding='utf-8') as f:
    json.dump(coords, f)
print('saved', len(coords), 'gene coords')
