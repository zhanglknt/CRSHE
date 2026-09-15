"""
Phase 1: 提取基因序列（仅10个灵长类物种）
"""

import os
import gzip
import csv
from collections import defaultdict
import sys

# 项目路径
PROJECT_DIR = 'C:\\Users\\admin\\Desktop\\人类正选择基因项目'
DATA_DIR = os.path.join(PROJECT_DIR, 'data')
RESULTS_DIR = os.path.join(PROJECT_DIR, 'results', 'phase1_gene_extraction')

# 创建结果目录
os.makedirs(RESULTS_DIR, exist_ok=True)

# 10个灵长类物种
PRIMATE_SPECIES = {
    'human': 'Homo sapiens',
    'chimpanzee': 'Pan troglodytes',
    'gorilla': 'Gorilla gorilla',
    'bonobo': 'Pan paniscus',
    'gibbon': 'Nomascus leucogenys',
    'marmoset': 'Callithrix jacchus',
    'vervet': 'Chlorocebus sabaeus',
    'rhesus_monkey': 'Macaca mulatta',
    'baboon': 'Papio anubis',
    'squirrel_monkey': 'Saimiri boliviensis'
}

def extract_genes_from_gtf(gtf_file, fasta_file, species_key):
    """从GTF和FASTA中提取基因CDS序列"""
    
    print(f"\nProcessing {species_key}...")
    print(f"  GTF: {gtf_file}")
    print(f"  FASTA: {fasta_file}")
    
    # 检查文件是否存在
    if not os.path.exists(gtf_file):
        print(f"  [X] GTF file not found")
        return 0, []
    
    if not os.path.exists(fasta_file):
        print(f"  [X] FASTA file not found")
        return 0, []
    
    # 1. 解析GTF文件，收集基因信息
    gene_info = []
    gene_cds_coords = defaultdict(list)  # gene_id -> [(chr, start, end, strand), ...]
    
    print(f"  Parsing GTF...")
    try:
        with open(gtf_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature_type = parts[2]
                if feature_type != 'CDS':
                    continue
                
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                
                # 解析属性
                attributes = parts[8]
                gene_id = None
                
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1]
                        break
                
                if gene_id:
                    gene_cds_coords[gene_id].append((chrom, start, end, strand))
        
        print(f"  Found {len(gene_cds_coords)} genes with CDS")
    except Exception as e:
        print(f"  [X] Error parsing GTF: {e}")
        return 0, []
    
    # 2. 读取FASTA文件
    print(f"  Reading genome FASTA...")
    genome_seqs = {}
    
    try:
        with open(fasta_file, 'r', encoding='utf-8', errors='ignore') as f:
            current_chrom = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_chrom:
                        genome_seqs[current_chrom] = ''.join(current_seq)
                    current_chrom = line[1:].split()[0]  # 取第一个词作为染色体名
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_chrom:
                genome_seqs[current_chrom] = ''.join(current_seq)
        
        print(f"  Loaded {len(genome_seqs)} chromosomes")
    except Exception as e:
        print(f"  [X] Error reading FASTA: {e}")
        return 0, []
    
    # 3. 提取每个基因的CDS序列
    print(f"  Extracting gene CDS sequences...")
    genes_cds = []
    gene_sequences = []
    
    for gene_id, coords in gene_cds_coords.items():
        # 排序坐标（按染色体和位置）
        coords_sorted = sorted(coords, key=lambda x: (x[0], x[1]))
        
        # 获取染色体
        chrom = coords_sorted[0][0]
        
        # 检查染色体是否存在
        if chrom not in genome_seqs:
            continue
        
        genome_seq = genome_seqs[chrom]
        strand = coords_sorted[0][3]
        
        # 合并所有CDS
        cdss = []
        for chr, start, end, strnd in coords_sorted:
            # 转换为0-based索引
            cds_seq = genome_seq[start-1:end]
            cdss.append(cds_seq)
        
        # 根据链方向合并
        if strand == '-':
            cdss = cdss[::-1]
            cds_sequence = ''.join(cdss)
            # 反向互补
            cds_sequence = reverse_complement(cds_sequence)
        else:
            cds_sequence = ''.join(cdss)
        
        # 过滤：只保留有CDS序列的基因
        if len(cds_sequence) > 0 and len(cds_sequence) % 3 == 0:
            genes_cds.append({
                'gene_id': gene_id,
                'chromosome': chrom,
                'strand': strand,
                'cds_length': len(cds_sequence)
            })
            gene_sequences.append((gene_id, cds_sequence))
    
    print(f"  Extracted {len(genes_cds)} genes with valid CDS")
    
    # 4. 保存结果
    # 保存CDS序列
    cds_file = os.path.join(RESULTS_DIR, f"{species_key}_genes_cds.fasta")
    with open(cds_file, 'w', encoding='utf-8') as f:
        for gene_id, seq in gene_sequences:
            f.write(f">{gene_id}\n{seq}\n")
    
    # 保存基因信息
    info_file = os.path.join(RESULTS_DIR, f"{species_key}_genes_info.csv")
    with open(info_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['gene_id', 'chromosome', 'strand', 'cds_length'])
        for gene in genes_cds:
            writer.writerow([gene['gene_id'], gene['chromosome'],
                           gene['strand'], gene['cds_length']])
    
    print(f"  Saved to {cds_file}")
    print(f"  Saved to {info_file}")
    
    return len(genes_cds), genes_cds

def reverse_complement(seq):
    """反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join([complement.get(base, base) for base in reversed(seq)])

def main():
    print("="*70)
    print("Phase 1: Extracting Gene Sequences from 10 Primate Species")
    print("="*70)
    
    # 查找文件
    species_results = {}
    
    for species_key, species_name in PRIMATE_SPECIES.items():
        species_dir = os.path.join(DATA_DIR, species_key)
        
        # 查找FASTA和GTF文件
        genome_files = []
        gtf_files = []
        
        if os.path.exists(species_dir):
            for filename in os.listdir(species_dir):
                if filename.endswith('.fa') and 'genome' in filename.lower():
                    genome_files.append(os.path.join(species_dir, filename))
                elif filename.endswith('.gtf') and ('gene' in filename.lower() or 'annotation' in filename.lower()):
                    gtf_files.append(os.path.join(species_dir, filename))
        
        if genome_files and gtf_files:
            print(f"\n{'='*70}")
            print(f"Species: {species_name} ({species_key})")
            print(f"{'='*70}")
            
            # 使用第一个找到的文件
            gene_count, genes = extract_genes_from_gtf(
                gtf_files[0],
                genome_files[0],
                species_key
            )
            
            species_results[species_key] = {
                'name': species_name,
                'gene_count': gene_count,
                'genes': genes,
                'genome_file': genome_files[0],
                'gtf_file': gtf_files[0]
            }
        else:
            print(f"\n[!] {species_key}: Missing files")
            if not os.path.exists(species_dir):
                print(f"    Directory not found")
            else:
                if not genome_files:
                    print(f"    No genome FASTA found")
                if not gtf_files:
                    print(f"    No GTF annotation found")
    
    # 保存汇总
    summary_file = os.path.join(RESULTS_DIR, 'phase1_summary.csv')
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['species', 'gene_count', 'genome_file', 'gtf_file'])
        
        for species_key, result in species_results.items():
            writer.writerow([
                species_key,
                result['gene_count'],
                result['genome_file'],
                result['gtf_file']
            ])
    
    # 打印总结
    print(f"\n{'='*70}")
    print("Summary")
    print(f"{'='*70}")
    print(f"\nTotal species processed: {len(species_results)}")
    print(f"\nSpecies list:")
    
    total_genes = 0
    for species_key, result in species_results.items():
        print(f"  {result['name']:30s} ({species_key:15s}): {result['gene_count']:6,} genes")
        total_genes += result['gene_count']
    
    print(f"\nTotal genes: {total_genes:,}")
    print(f"\nResults saved to: {RESULTS_DIR}")
    print(f"Summary file: {summary_file}")
    
    print(f"\n{'='*70}")
    print("Phase 1 Complete!")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
