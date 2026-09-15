# -*- coding: utf-8 -*-
"""
Phase 7A (v3): PAML Branch-Site Model Batch Analysis (Short Genes Only)
只分析≤1000密码子的基因，增加timeout确保完成
"""
import os
import sys
import csv
import json
import time
import shutil
import subprocess
import multiprocessing
from pathlib import Path
from datetime import datetime

# ===== 路径配置 =====
PROJECT_ROOT = r"d:\人类正选择基因项目"
CODEML      = os.path.join(PROJECT_ROOT, r"tools\paml\paml-4.10.10-win-x86_64\bin\codeml.exe")
PHY_DIR     = os.path.join(PROJECT_ROOT, r"results\phase3_msa\phylip_paml_cleaned_v3_short")
TREE_FILE   = os.path.join(PROJECT_ROOT, r"results\species_tree\species_tree_binary.nwk")
OUT_DIR     = os.path.join(PROJECT_ROOT, r"results\phase7_gene_vs_regulation\phase7a_paml_results_v3")
LOG_DIR     = os.path.join(OUT_DIR, "logs")
STATUS_FILE = os.path.join(OUT_DIR, "run_status.json")
TEMP_BASE   = r"C:\paml_tmp_v4"

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR,  exist_ok=True)
os.makedirs(TEMP_BASE, exist_ok=True)

N_WORKERS = min(8, multiprocessing.cpu_count())

# ===== 动态超时（v3: 增加超时确保完成）=====
# 所有基因≤1000密码子，timeout = 30-60分钟
def calc_timeout(n_codons):
    # 1000密码子 → 60分钟(3600秒) 超时
    # 500密码子 → 30分钟(1800秒) 超时
    t = int(n_codons / 1000 * 3600 + 600)  # 基础600秒 + 按比例增加
    return min(t, 7200)  # 最大2小时

# ===== 解析lnL =====
def parse_lnL(out_file):
    try:
        with open(out_file, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
        for line in content.split('\n'):
            if 'lnL' in line:
                if ':' in line:
                    try:
                        val = line.split(':')[1].strip().split()[0]
                        return float(val)
                    except:
                        pass
                if '=' in line:
                    try:
                        val = line.split('=')[1].strip().split()[0]
                        return float(val)
                    except:
                        pass
    except:
        pass
    return None

# ===== 读取基因列表 =====
def load_gene_list():
    genes = []
    for fname in sorted(os.listdir(PHY_DIR)):
        if fname.endswith('.phy'):
            gene_id = fname.replace('.phy', '')
            phy_path = os.path.join(PHY_DIR, fname)
            try:
                with open(phy_path) as f:
                    header = f.readline().strip().split()
                    n_sites = int(header[1])
                    n_codons = n_sites // 3
                    if n_codons >= 30:
                        genes.append((gene_id, n_codons, phy_path))
            except:
                pass
    return genes

# ===== 读取已完成基因（检查所有历史目录）=====
OLD_DIRS = [
    os.path.join(PROJECT_ROOT, r"results\phase7_gene_vs_regulation\phase7a_paml_results"),
    os.path.join(PROJECT_ROOT, r"results\phase7_gene_vs_regulation\phase7a_paml_results_v2"),
]

def load_done():
    done = {}
    for out_dir in OLD_DIRS:
        if os.path.exists(out_dir):
            for fname in os.listdir(out_dir):
                if fname.endswith('_H1.out'):
                    gene_id = fname.replace('_H1.out', '')
                    h0_file = os.path.join(out_dir, f"{gene_id}_H0.out")
                    h1_file = os.path.join(out_dir, f"{gene_id}_H1.out")
                    if os.path.exists(h0_file) and os.path.exists(h1_file):
                        h0_lnL = parse_lnL(h0_file)
                        h1_lnL = parse_lnL(h1_file)
                        if h0_lnL is not None and h1_lnL is not None:
                            done[gene_id] = (h0_lnL, h1_lnL)
    # 检查v3目录
    if os.path.exists(OUT_DIR):
        for fname in os.listdir(OUT_DIR):
            if fname.endswith('_H1.out'):
                gene_id = fname.replace('_H1.out', '')
                h0_file = os.path.join(OUT_DIR, f"{gene_id}_H0.out")
                h1_file = os.path.join(OUT_DIR, f"{gene_id}_H1.out")
                if os.path.exists(h0_file) and os.path.exists(h1_file):
                    h0_lnL = parse_lnL(h0_file)
                    h1_lnL = parse_lnL(h1_file)
                    if h0_lnL is not None and h1_lnL is not None:
                        done[gene_id] = (h0_lnL, h1_lnL)
    return done

# ===== 写状态文件 =====
def write_status(done_count, failed_count, total, start_ts):
    elapsed = max(time.time() - start_ts, 0.1)
    speed = done_count / elapsed * 3600
    remain = (total - done_count - failed_count) / (done_count / elapsed) if done_count > 0 else 0
    data = {
        "total": total,
        "done": done_count,
        "failed": failed_count,
        "pct": round(done_count / max(total, 1) * 100, 2),
        "elapsed_h": round(elapsed / 3600, 2),
        "speed_per_h": round(speed, 1),
        "eta_h": round(remain / 3600, 2),
        "updated": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    tmp = STATUS_FILE + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    shutil.move(tmp, STATUS_FILE)

# ===== 单基因分析（子进程） =====
def run_one_gene(args):
    gene_id, n_codons, phy_path = args
    timeout = calc_timeout(n_codons)
    
    results = {}
    lnLs = {}
    
    for model in ("H0", "H1"):
        out_dst = os.path.join(OUT_DIR, f"{gene_id}_{model}.out")
        if os.path.exists(out_dst) and os.path.getsize(out_dst) > 100:
            existing_lnL = parse_lnL(out_dst)
            if existing_lnL is not None:
                results[model] = "skip"
                lnLs[model] = existing_lnL
                continue
        
        tmp_id = f"{gene_id}_{model}_{os.getpid()}"
        tmp_dir = os.path.join(TEMP_BASE, tmp_id)
        os.makedirs(tmp_dir, exist_ok=True)
        
        phy_name = "seq.phy"
        tree_name = "tree.nwk"
        out_name = "result.out"
        ctl_name = "run.ctl"
        
        try:
            shutil.copy2(phy_path, os.path.join(tmp_dir, phy_name))
            shutil.copy2(TREE_FILE, os.path.join(tmp_dir, tree_name))
            
            ctl = (
                f"seqfile  = {phy_name}\n"
                f"outfile  = {out_name}\n"
                f"treefile = {tree_name}\n"
                "noisy    = 0\n"
                "verbose  = 0\n"
                "runmode  = 0\n"
                "seqtype  = 1\n"
                "CodonFreq = 2\n"
                "clock    = 0\n"
                "aaDist   = 0\n"
                "model    = 2\n"
                "NSsites  = 2\n"
                "icode    = 0\n"
                "Mgene    = 0\n"
                "fix_kappa = 0\n"
                "kappa    = 2\n"
                + ("fix_omega = 1\n" if model == "H0" else "fix_omega = 0\n") +
                "omega    = 1\n"
                "fix_alpha = 1\n"
                "alpha    = 0\n"
                "Malpha   = 0\n"
                "ncatG    = 10\n"
                "getSE    = 0\n"
                "RateAncestor = 0\n"
                "Small_Diff = .5e-6\n"
                "cleandata = 0\n"
            )
            with open(os.path.join(tmp_dir, ctl_name), "w", encoding="ascii") as f:
                f.write(ctl)
            
            proc = subprocess.run(
                [CODEML, ctl_name],
                cwd=tmp_dir, timeout=timeout,
                capture_output=True, text=True,
                encoding="utf-8", errors="replace"
            )
            
            actual = os.path.join(tmp_dir, out_name)
            if os.path.exists(actual) and os.path.getsize(actual) > 100:
                shutil.copy2(actual, out_dst)
                lnL = parse_lnL(out_dst)
                if lnL is not None:
                    results[model] = "ok"
                    lnLs[model] = lnL
                else:
                    results[model] = "no_lnL"
            else:
                results[model] = "fail"
                log = os.path.join(LOG_DIR, f"{gene_id}_{model}_fail.txt")
                with open(log, "w", encoding="utf-8", errors="replace") as lf:
                    lf.write(proc.stdout + "\n" + proc.stderr)
        except subprocess.TimeoutExpired:
            results[model] = "timeout"
            with open(os.path.join(LOG_DIR, f"{gene_id}_{model}_timeout.txt"), "w") as lf:
                lf.write(f"timeout {timeout}s\n")
        except Exception as e:
            results[model] = f"err:{e}"
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)
    
    ok = all(v in ("ok", "skip") for v in results.values())
    return gene_id, ok, results, lnLs

# ===== 主程序 =====
def main():
    print("[Phase 7A v3] PAML Branch-Site 批量分析 (短基因优化版)")
    print(f"  codeml  : {CODEML}")
    print(f"  PHY目录 : {PHY_DIR}")
    print(f"  树文件  : {TREE_FILE}")
    print(f"  输出目录: {OUT_DIR}")
    print(f"  临时目录: {TEMP_BASE}")
    print(f"  并行进程: {N_WORKERS}")
    
    genes = load_gene_list()
    already = load_done()
    todo = [g for g in genes if g[0] not in already]
    total = len(genes)
    done = len(already)
    failed = 0
    
    print(f"  有效基因: {total}")
    print(f"  已完成:   {done}")
    print(f"  待分析:   {len(todo)}")
    
    if not todo:
        print("[DONE] 全部分析已完成！")
        return
    
    start_ts = time.time()
    write_status(done, failed, total, start_ts)
    
    print(f"\n  [开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]")
    print("  格式: [时间] 进度% (完成/总数) 速度/h ETA")
    
    with multiprocessing.Pool(N_WORKERS) as pool:
        for i, (gene_id, ok, res, lnLs) in enumerate(pool.imap_unordered(run_one_gene, todo), 1):
            if ok:
                done += 1
            else:
                failed += 1
            
            if i % 10 == 0 or i == len(todo):
                elapsed = time.time() - start_ts
                speed = done / elapsed * 3600
                remain = (len(todo) - i) / (i / elapsed) if i > 0 else 0
                pct = done / total * 100
                print(f"  [{datetime.now().strftime('%H:%M:%S')}] {pct:.1f}% ({done}/{total}) "
                      f"速度:{speed:.1f}/h ETA:{remain/3600:.1f}h 失败:{failed}")
                write_status(done, failed, total, start_ts)
    
    print(f"\n[DONE] 完成: {done}/{total}, 失败: {failed}")
    print(f"  总耗时: {(time.time()-start_ts)/3600:.2f}小时")
    write_status(done, failed, total, start_ts)

if __name__ == '__main__':
    main()
