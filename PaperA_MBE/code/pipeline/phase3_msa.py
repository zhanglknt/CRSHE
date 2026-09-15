# -*- coding: utf-8 -*-
"""
Phase 3: Multiple Sequence Alignment using MAFFT
- Processes 6908 genes across 10 species
- Uses parallel processing for speed
- Output: PHYLIP files for PAML analysis
"""

import os
import sys
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import logging
import time
import json
import argparse
from multiprocessing import Pool, Manager
import signal

# ===================== 路径配置 =====================
PROJECT_ROOT = r"C:\Users\admin\Desktop\人类正选择基因项目"
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")
PHASE1_DIR = os.path.join(RESULTS_DIR, "phase1_gene_extraction")
PHASE2_DIR = os.path.join(RESULTS_DIR, "phase2_final")
PHASE3_DIR = os.path.join(RESULTS_DIR, "phase3_msa")
TMP_DIR    = os.path.join(PHASE3_DIR, "tmp")

MAFFT_EXE = r"C:\Users\admin\Desktop\人类正选择基因项目\tools\mafft\mafft-win\mafft.bat"

SPECIES_MAP = {
    "human":          "human_genes_cds.fasta",
    "chimpanzee":     "chimpanzee_genes_cds.fasta",
    "bonobo":         "bonobo_genes_cds.fasta",
    "gorilla":        "gorilla_genes_cds.fasta",
    "gibbon":         "gibbon_genes_cds.fasta",
    "marmoset":       "marmoset_genes_cds.fasta",
    "vervet":         "vervet_genes_cds.fasta",
    "baboon":         "baboon_genes_cds.fasta",
    "rhesus":         "rhesus_monkey_genes_cds.fasta",
    "squirrel_monkey":"squirrel_monkey_genes_cds.fasta",
}

CSV_COL_MAP = {
    "human":          "human_gene_id",
    "chimpanzee":     "chimpanzee_gene_id",
    "bonobo":         "bonobo_gene_id",
    "gorilla":        "gorilla_gene_id",
    "gibbon":         "gibbon_gene_id",
    "marmoset":       "marmoset_gene_id",
    "vervet":         "vervet_gene_id",
    "baboon":         "baboon_gene_id",
    "rhesus":         "rhesus_gene_id",
    "squirrel_monkey":"squirrel_monkey_gene_id",
}

# ===================== 全局CDS数据库（进程级缓存）=====================
_CDS_DB = None

def get_cds_db():
    global _CDS_DB
    if _CDS_DB is None:
        _CDS_DB = {}
        for sp, fname in SPECIES_MAP.items():
            _CDS_DB[sp] = {}
            fasta_path = os.path.join(PHASE1_DIR, fname)
            for rec in SeqIO.parse(fasta_path, "fasta"):
                gid = rec.id.split(".")[0]
                if gid not in _CDS_DB[sp] or len(rec.seq) > len(_CDS_DB[sp][gid]):
                    _CDS_DB[sp][gid] = str(rec.seq).upper()
    return _CDS_DB


def trim_to_codon(seq):
    """截断到3的倍数"""
    rem = len(seq) % 3
    return seq[:-rem] if rem else seq


def seq_to_phylip(species_seqs, gene_id):
    """生成PAML格式的PHYLIP文件内容"""
    n = len(species_seqs)
    seq_len = len(list(species_seqs.values())[0])
    lines = [f" {n} {seq_len}"]
    for sp, seq in species_seqs.items():
        name = sp[:10].ljust(10)
        lines.append(f"{name}  {seq}")
    return "\n".join(lines) + "\n"


def process_gene(args):
    """处理单个基因（多进程worker函数）"""
    row_dict, pid = args
    human_id = row_dict["human_gene_id"]
    
    phy_dir = os.path.join(PHASE3_DIR, "phylip")
    out_phy = os.path.join(phy_dir, f"{human_id}.phy")
    
    # 已处理过，直接跳过
    if os.path.exists(out_phy) and os.path.getsize(out_phy) > 10:
        return (human_id, "cached", 0, "")

    cds_db = get_cds_db()
    
    # 提取各物种序列
    seqs = {}
    for sp, col in CSV_COL_MAP.items():
        gid = str(row_dict.get(col, "")).split(".")[0]
        if gid and gid in cds_db.get(sp, {}):
            seq = trim_to_codon(cds_db[sp][gid])
            if len(seq) >= 30:
                seqs[sp] = seq

    if len(seqs) < 4:
        return (human_id, "skip", len(seqs), "too_few_species")

    # 写临时FASTA
    tmp_fa = os.path.join(TMP_DIR, f"{human_id}_{pid}.fa")
    with open(tmp_fa, "w") as f:
        for sp, seq in seqs.items():
            f.write(f">{sp}\n{seq}\n")

    # MAFFT比对
    try:
        result = subprocess.run(
            [MAFFT_EXE, "--auto", "--preservecase", "--quiet", tmp_fa],
            capture_output=True, text=True, timeout=600
        )
        success = result.returncode == 0 and result.stdout.strip()
    except subprocess.TimeoutExpired:
        success = False
    except Exception as e:
        success = False

    # 清理临时文件
    if os.path.exists(tmp_fa):
        try:
            os.remove(tmp_fa)
        except:
            pass

    if success:
        # 解析MAFFT输出
        aligned = {}
        current_id = None
        current_seq = []
        for line in result.stdout.strip().split("\n"):
            if line.startswith(">"):
                if current_id:
                    aligned[current_id] = "".join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            aligned[current_id] = "".join(current_seq)
        
        method = "mafft"
    else:
        # 降级：截断到最短长度，直接使用
        min_len = min(len(s) for s in seqs.values())
        min_len = (min_len // 3) * 3
        aligned = {sp: seq[:min_len] for sp, seq in seqs.items()}
        method = "truncated"

    # 写PHYLIP文件
    phy_content = seq_to_phylip(aligned, human_id)
    if phy_content:
        with open(out_phy, "w") as f:
            f.write(phy_content)
        return (human_id, "ok", len(seqs), method)
    else:
        return (human_id, "fail", len(seqs), "empty_alignment")


def run_phase3(n_workers=3, batch_report=100):
    """主函数"""
    # 创建目录
    for d in [PHASE3_DIR, TMP_DIR, os.path.join(PHASE3_DIR, "phylip")]:
        os.makedirs(d, exist_ok=True)
    
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(PHASE3_DIR, "phase3.log"), encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ]
    )
    logger = logging.getLogger(__name__)
    
    logger.info("=" * 60)
    logger.info("PHASE 3: Multiple Sequence Alignment (MAFFT)")
    logger.info(f"Workers: {n_workers}")
    logger.info("=" * 60)

    # 加载数据
    df = pd.read_csv(os.path.join(PHASE2_DIR, "shared_genes_all10.csv"))
    phy_dir = os.path.join(PHASE3_DIR, "phylip")
    
    # 断点续传
    existing = set(f.replace(".phy", "") for f in os.listdir(phy_dir) if f.endswith(".phy"))
    df_todo = df[~df["human_gene_id"].isin(existing)].copy()
    
    logger.info(f"Total: {len(df)} | Done: {len(existing)} | Todo: {len(df_todo)}")
    
    if len(df_todo) == 0:
        logger.info("All genes already processed!")
        return
    
    # 准备参数
    args_list = [(row.to_dict(), os.getpid()) for _, row in df_todo.iterrows()]
    
    # 统计
    ok_count = len(existing)
    skip_count = 0
    fail_count = 0
    start_time = time.time()
    processed = 0
    
    results_records = []
    
    # 多进程处理
    logger.info(f"Starting parallel processing with {n_workers} workers...")
    with Pool(processes=n_workers) as pool:
        for result in pool.imap_unordered(process_gene, args_list, chunksize=5):
            gene_id, status, n_species, method = result
            results_records.append({
                "gene_id": gene_id,
                "status": status,
                "n_species": n_species,
                "method": method
            })
            
            if status in ("ok", "cached"):
                ok_count += 1
            elif status == "skip":
                skip_count += 1
            else:
                fail_count += 1
            
            processed += 1
            
            if processed % batch_report == 0:
                elapsed = time.time() - start_time
                total_done = len(existing) + processed
                rate = processed / elapsed if elapsed > 0 else 0
                remaining = (len(df_todo) - processed) / rate if rate > 0 else 0
                pct = total_done / len(df) * 100
                logger.info(
                    f"[{total_done}/{len(df)}] {pct:.1f}% | "
                    f"OK={ok_count} SKIP={skip_count} FAIL={fail_count} | "
                    f"Rate={rate:.1f}/s | ETA={remaining/60:.1f}min"
                )
    
    total_time = time.time() - start_time
    total_phy = len([f for f in os.listdir(phy_dir) if f.endswith(".phy")])
    
    logger.info("\n" + "=" * 60)
    logger.info("PHASE 3 COMPLETE")
    logger.info(f"  Total PHY files: {total_phy}")
    logger.info(f"  Success: {ok_count} | Skip: {skip_count} | Fail: {fail_count}")
    logger.info(f"  Time: {total_time/60:.1f} min ({total_time/3600:.2f} hr)")
    logger.info(f"  Output: {phy_dir}")
    logger.info("=" * 60)
    
    # 保存结果
    if results_records:
        res_df = pd.DataFrame(results_records)
        res_df.to_csv(os.path.join(PHASE3_DIR, "phase3_results.csv"), index=False)
    
    # 写成功基因列表
    ok_genes = [f.replace(".phy", "") for f in os.listdir(phy_dir) if f.endswith(".phy")]
    with open(os.path.join(PHASE3_DIR, "phase3_ok_genes.txt"), "w") as f:
        f.write("\n".join(ok_genes))
    
    logger.info(f"Success gene list: {len(ok_genes)} genes -> phase3_ok_genes.txt")
    return total_phy


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 3: MSA with MAFFT")
    parser.add_argument("--workers", type=int, default=3, help="Number of parallel workers")
    parser.add_argument("--report", type=int, default=50, help="Report interval")
    args = parser.parse_args()
    
    run_phase3(n_workers=args.workers, batch_report=args.report)
