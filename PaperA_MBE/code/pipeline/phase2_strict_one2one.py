# -*- coding: utf-8 -*-
"""
Phase 2 v2: Strict one-to-one ortholog identification
- Query Ensembl REST API with type=orthologues
- Client-side filter: ONLY ortholog_one2one (no fallback to one2many/many2many)
- Multi-threaded for speed (10 concurrent workers)
- Resume capability (saves progress every 500 genes)
- Output: results/phase2_final_v2/
"""

import requests
import time
import json
import csv
import os
import sys
import logging
from pathlib import Path
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Encoding
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# ===================== Paths =====================
PROJECT_ROOT = Path(r"d:\人类正选择基因项目")
RESULTS_DIR = PROJECT_ROOT / "results" / "phase2_final_v2"
PHASE1_DIR = PROJECT_ROOT / "results" / "phase1_gene_extraction"
OLD_RESULTS = PROJECT_ROOT / "results" / "phase2_final"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ===================== Ensembl REST API =====================
ENSEMBL_REST = "https://rest.ensembl.org"

# Target species (9 non-human primates)
TARGET_SPECIES = [
    "pan_troglodytes",
    "pan_paniscus",
    "gorilla_gorilla",
    "nomascus_leucogenys",
    "callithrix_jacchus",
    "chlorocebus_sabaeus",
    "papio_anubis",
    "macaca_mulatta",
    "saimiri_boliviensis_boliviensis",
]

SP_DISPLAY = {
    "pan_troglodytes": "chimpanzee",
    "pan_paniscus": "bonobo",
    "gorilla_gorilla": "gorilla",
    "nomascus_leucogenys": "gibbon",
    "callithrix_jacchus": "marmoset",
    "chlorocebus_sabaeus": "vervet",
    "papio_anubis": "baboon",
    "macaca_mulatta": "rhesus",
    "saimiri_boliviensis_boliviensis": "squirrel_monkey",
}

TARGET_SET = set(TARGET_SPECIES)

# Progress file for resume
PROGRESS_FILE = RESULTS_DIR / "progress.json"
PROGRESS_LOCK = threading.Lock()

# Logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(RESULTS_DIR / "phase2_v2.log", encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)


def load_human_genes():
    """Load human gene IDs from Phase 1 FASTA"""
    fasta_file = PHASE1_DIR / "human_genes_cds.fasta"
    genes = []
    for rec in SeqIO.parse(fasta_file, "fasta"):
        gene_id = rec.id.split("|")[0].split(".")[0]
        gene_name = rec.id.split("|")[1] if "|" in rec.id else ""
        genes.append({"gene_id": gene_id, "gene_name": gene_name})
    return genes


def load_progress():
    """Load progress for resume"""
    if PROGRESS_FILE.exists():
        with open(PROGRESS_FILE, 'r') as f:
            return json.load(f)
    return {"completed": {}, "all10": [], "failed": []}


def save_progress(progress):
    """Save progress atomically"""
    tmp = str(PROGRESS_FILE) + ".tmp"
    with open(tmp, 'w') as f:
        json.dump(progress, f)
    os.replace(tmp, PROGRESS_FILE)


def query_one2one_orthologs(gene_id, max_retries=3):
    """
    Query Ensembl REST API for orthologs.
    Returns: (gene_id, {species: ortholog_id}, status)
    Only returns ortholog_one2one entries - NO fallback to other types.
    """
    url = f"{ENSEMBL_REST}/homology/id/homo_sapiens/{gene_id}"
    params = {"type": "orthologues"}  # Only orthologs, not paralogs

    for attempt in range(max_retries):
        try:
            r = requests.get(
                url,
                headers={"Content-Type": "application/json"},
                params=params,
                timeout=30
            )

            if r.status_code == 429:
                # Rate limited - exponential backoff
                wait = min(2 ** attempt, 10)
                time.sleep(wait)
                continue

            if r.status_code == 404:
                # Gene not found in Ensembl
                return gene_id, {}, "not_found"

            if r.status_code != 200:
                return gene_id, {}, f"http_{r.status_code}"

            d = r.json()

            # Parse - ONLY ortholog_one2one
            one2one = {}
            for entry in d.get("data", []):
                for h in entry.get("homologies", []):
                    sp = h.get("target", {}).get("species", "")
                    htype = h.get("type", "")
                    tid = h["target"]["id"]

                    # Strict filter: only ortholog_one2one in target species
                    if sp in TARGET_SET and htype == "ortholog_one2one":
                        # one2one means exactly one per species, but API might
                        # return duplicates - keep first
                        if sp not in one2one:
                            one2one[sp] = tid

            return gene_id, one2one, "ok"

        except requests.exceptions.Timeout:
            if attempt < max_retries - 1:
                time.sleep(2)
                continue
            return gene_id, {}, "timeout"
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(1)
                continue
            return gene_id, {}, f"error:{str(e)[:50]}"

    return gene_id, {}, "max_retries"


def process_gene(gene_info):
    """Process a single gene - wrapper for threading"""
    gene_id = gene_info["gene_id"]
    gene_name = gene_info["gene_name"]
    gid, one2one, status = query_one2one_orthologs(gene_id)

    result = {
        "gene_id": gene_id,
        "gene_name": gene_name,
        "status": status,
        "n_one2one": len(one2one),
        "all9": len(one2one) == 9,
        "orthologs": one2one,
    }
    return result


def main():
    logger.info("=" * 60)
    logger.info("Phase 2 v2: Strict one-to-one ortholog identification")
    logger.info("Only ortholog_one2one - NO fallback to one2many/many2many")
    logger.info("=" * 60)

    # Load human genes
    human_genes = load_human_genes()
    total = len(human_genes)
    logger.info(f"Loaded {total} human genes from Phase 1")

    # Load progress for resume
    progress = load_progress()
    completed_ids = set(progress["completed"].keys())
    logger.info(f"Already completed: {len(completed_ids)} | Remaining: {total - len(completed_ids)}")

    # Filter to remaining genes
    todo = [g for g in human_genes if g["gene_id"] not in completed_ids]

    if not todo:
        logger.info("All genes already processed! Generating output...")
    else:
        logger.info(f"Processing {len(todo)} genes with 10 concurrent workers...")

        # Multi-threaded processing
        n_workers = 10
        batch_save_interval = 500
        processed_since_save = 0

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(process_gene, g): g for g in todo}

            for future in as_completed(futures):
                result = future.result()
                gid = result["gene_id"]

                with PROGRESS_LOCK:
                    progress["completed"][gid] = {
                        "status": result["status"],
                        "n_one2one": result["n_one2one"],
                        "all9": result["all9"],
                        "orthologs": result["orthologs"],
                        "gene_name": result["gene_name"],
                    }

                    if result["all9"]:
                        progress["all10"].append(gid)
                    elif result["status"] != "ok":
                        progress["failed"].append(gid)

                    processed_since_save += 1

                    if processed_since_save >= batch_save_interval:
                        save_progress(progress)
                        processed_since_save = 0
                        n_done = len(progress["completed"])
                        n_all10 = len(progress["all10"])
                        n_failed = len(progress["failed"])
                        logger.info(
                            f"  [{n_done}/{total}] {n_done/total*100:.1f}% | "
                            f"all9={n_all10} failed={n_failed}"
                        )

        # Final save
        save_progress(progress)

    # Generate output files
    logger.info("\nGenerating output files...")

    # CSV
    out_csv = RESULTS_DIR / "shared_genes_one2one_all10.csv"
    fieldnames = ["human_gene_id", "human_gene_name"] + \
                 [f"{SP_DISPLAY[sp]}_gene_id" for sp in TARGET_SPECIES]

    all10_genes = []
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()

        for gid, info in progress["completed"].items():
            if info["all9"] and info["status"] == "ok":
                row = {"human_gene_id": gid, "human_gene_name": info.get("gene_name", "")}
                orthologs = info["orthologs"]
                for sp in TARGET_SPECIES:
                    row[f"{SP_DISPLAY[sp]}_gene_id"] = orthologs.get(sp, "")
                writer.writerow(row)
                all10_genes.append({
                    "gid": gid,
                    "gname": info.get("gene_name", ""),
                    "homologs": orthologs,
                })

    # JSON
    json_out = RESULTS_DIR / "shared_genes_one2one_all10.json"
    with open(json_out, "w", encoding="utf-8") as f:
        json.dump(all10_genes, f, indent=2, ensure_ascii=False)

    # Summary
    n_ok = sum(1 for v in progress["completed"].values() if v["status"] == "ok")
    n_failed = sum(1 for v in progress["completed"].values() if v["status"] != "ok")
    n_all10 = len(all10_genes)

    # Species coverage stats
    sp_counts = {sp: 0 for sp in TARGET_SPECIES}
    for g in all10_genes:
        for sp in g["homologs"]:
            sp_counts[sp] += 1

    # Compare with old results
    old_csv = OLD_RESULTS / "shared_genes_all10.csv"
    old_count = 0
    old_ids = set()
    if old_csv.exists():
        with open(old_csv, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                old_count += 1
                old_ids.add(row["human_gene_id"])

    new_ids = set(g["gid"] for g in all10_genes)
    overlap = old_ids & new_ids
    only_old = old_ids - new_ids
    only_new = new_ids - old_ids

    summary = []
    summary.append("=" * 60)
    summary.append("PHASE 2 v2 SUMMARY: Strict ortholog_one2one")
    summary.append("=" * 60)
    summary.append("")
    summary.append(f"Total human genes queried  : {total}")
    summary.append(f"Successful queries         : {n_ok}")
    summary.append(f"Failed queries             : {n_failed}")
    summary.append(f"Genes with one2one in ALL 9: {n_all10}")
    summary.append(f"Success rate (all9)        : {n_all10/total*100:.1f}%")
    summary.append("")
    summary.append("Species coverage (one2one orthologs):")
    for sp in TARGET_SPECIES:
        cnt = sp_counts[sp]
        summary.append(f"  {SP_DISPLAY[sp]:20s}: {cnt} ({cnt/n_all10*100:.1f}%)")
    summary.append("")
    summary.append("=" * 60)
    summary.append("COMPARISON WITH OLD RESULTS (v1)")
    summary.append("=" * 60)
    summary.append(f"Old v1 (with fallback)     : {old_count} genes")
    summary.append(f"New v2 (strict one2one)    : {n_all10} genes")
    summary.append(f"Overlap (in both)          : {len(overlap)} genes")
    summary.append(f"Only in v1 (removed)       : {len(only_old)} genes")
    summary.append(f"Only in v2 (newly found)   : {len(only_new)} genes")
    summary.append(f"Reduction                  : {old_count - n_all10} genes ({(old_count - n_all10)/old_count*100:.1f}%)")
    summary.append("")
    summary.append("Status distribution of removed genes:")
    removed_statuses = {}
    for gid in only_old:
        if gid in progress["completed"]:
            s = progress["completed"][gid]["status"]
            n = progress["completed"][gid]["n_one2one"]
            removed_statuses.setdefault(f"{s} (n_one2one={n})", 0)
            removed_statuses[f"{s} (n_one2one={n})"] += 1
        else:
            removed_statuses.setdefault("not_in_v2", 0)
            removed_statuses["not_in_v2"] += 1
    for status, count in sorted(removed_statuses.items(), key=lambda x: -x[1])[:20]:
        summary.append(f"  {count:5d}  {status}")

    summary_text = "\n".join(summary)
    summary_out = RESULTS_DIR / "phase2_v2_summary.txt"
    with open(summary_out, "w", encoding="utf-8") as f:
        f.write(summary_text)

    logger.info("\n" + summary_text)
    logger.info(f"\nOutput: {RESULTS_DIR}")


if __name__ == "__main__":
    main()
