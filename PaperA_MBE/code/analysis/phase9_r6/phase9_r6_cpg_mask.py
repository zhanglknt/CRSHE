#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate CpG-masked alignments for the 400-gene stratified subset.

Identifies CpG dinucleotides in the human reference sequence (including
CpG sites spanning codon boundaries) and masks the affected codon columns
(replaced with ???) across ALL species in the alignment.

Usage (WSL): python3 phase9_r6_cpg_mask.py <genes_csv> <fasta_dir> <out_dir> <stats_csv>
"""
import csv
import json
import sys
from pathlib import Path

genes_csv, fasta_dir, out_dir, stats_csv = sys.argv[1:5]
fasta_dir = Path(fasta_dir)
out_dir = Path(out_dir)
out_dir.mkdir(parents=True, exist_ok=True)


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if name:
        seqs[name] = "".join(chunks).upper()
    return seqs


def write_fasta(path, seqs, width=60):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")


rows = []
genes = [r["gene_id"] for r in csv.DictReader(open(genes_csv, encoding="utf-8"))]
for gid in genes:
    src = fasta_dir / f"{gid}.fasta"
    if not src.exists():
        rows.append({"gene_id": gid, "status": "missing_fasta"})
        continue
    seqs = read_fasta(src)
    if "human" not in seqs:
        rows.append({"gene_id": gid, "status": "no_human"})
        continue
    L = {len(s) for s in seqs.values()}
    if len(L) != 1:
        rows.append({"gene_id": gid, "status": "ragged"})
        continue
    L = L.pop()
    if L % 3 != 0:
        rows.append({"gene_id": gid, "status": "not_multiple_of_3"})
        continue
    ncod = L // 3
    h = seqs["human"]
    # find CpG positions in human (skip ambiguous chars)
    mask_codons = set()
    n_cpg = 0
    for i in range(L - 1):
        if h[i] == "C" and h[i + 1] == "G":
            n_cpg += 1
            mask_codons.add(i // 3)
            mask_codons.add((i + 1) // 3)
    masked = dict(seqs)
    for sp, s in seqs.items():
        chars = list(s)
        for c in mask_codons:
            chars[c * 3:c * 3 + 3] = "???"
        masked[sp] = "".join(chars)
    write_fasta(out_dir / f"{gid}.fasta", masked)
    rows.append({
        "gene_id": gid, "status": "ok", "n_codons": ncod,
        "n_cpg_sites_human": n_cpg, "n_masked_codons": len(mask_codons),
        "frac_codons_masked": round(len(mask_codons) / ncod, 5),
    })

with open(stats_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

ok = [r for r in rows if r["status"] == "ok"]
summary = {
    "n_genes": len(rows),
    "n_ok": len(ok),
    "mean_frac_codons_masked": sum(r["frac_codons_masked"] for r in ok) / len(ok) if ok else None,
    "median_frac_codons_masked": sorted(r["frac_codons_masked"] for r in ok)[len(ok) // 2] if ok else None,
    "mean_n_cpg_sites": sum(r["n_cpg_sites_human"] for r in ok) / len(ok) if ok else None,
    "status_counts": {s: sum(1 for r in rows if r["status"] == s) for s in {r["status"] for r in rows}},
}
print(json.dumps(summary, indent=2))
Path(stats_csv).with_suffix(".summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
