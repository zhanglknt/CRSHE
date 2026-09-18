#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Neutral simulation with PAML evolver (codon model, omega=1) for 400 sampled genes.

For each sampled real gene:
  1. read its alignment (length + per-species gap mask + species set)
  2. simulate a codon alignment on the SAME tree (branch lengths, {Test} stripped)
     with omega=1, kappa=2, F3x4 codon frequencies estimated from the real data
  3. apply the real gene's gap mask and species-missingness pattern
  4. write FASTA for BUSTED

Usage (WSL): python3 phase9_r6_evolver.py <sim_genes_csv> <fasta_dir> <tree_dir> <freq_json> <out_dir> <work_dir> <evolver_bin>
"""
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sim_csv, fasta_dir, tree_dir, freq_json, out_dir, work_dir, evolver_bin = sys.argv[1:8]
fasta_dir, tree_dir = Path(fasta_dir), Path(tree_dir)
out_dir, work_dir = Path(out_dir), Path(work_dir)
out_dir.mkdir(parents=True, exist_ok=True)
work_dir.mkdir(parents=True, exist_ok=True)

freqs = json.load(open(freq_json, encoding="utf-8"))["codon_freqs_f3x4_stops_zeroed_renorm"]
# renormalize to exactly 1.0 at 10 dp
freqs = [round(f / sum(freqs), 10) for f in freqs]
freqs[-1] = round(freqs[-1] + (1 - sum(freqs)), 10)


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


def parse_mc(path):
    """Parse evolver mc.txt (paml format, possibly interleaved)."""
    lines = [l.rstrip("\n") for l in open(path)]
    # header: two ints on one line
    hdr_idx = None
    for i, l in enumerate(lines):
        t = l.split()
        if len(t) == 2 and all(x.isdigit() for x in t):
            hdr_idx = i
            break
    nseq, nchar = (int(x) for x in lines[hdr_idx].split())
    seqs = {}
    order = []
    cur = None
    for l in lines[hdr_idx + 1:]:
        if not l.strip():
            continue
        t = l.split()
        if len(t) >= 2 and not t[0].isdigit() and set(t[0]) <= set(
                "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789"):
            name = t[0]
            if name not in seqs:
                seqs[name] = ""
                order.append(name)
            cur = name
            seqs[name] += "".join(t[1:])
        elif t and cur is not None:
            seqs[cur] += "".join(t)
    return seqs


rows = []
genes = list(csv.DictReader(open(sim_csv, encoding="utf-8")))
for i, r in enumerate(genes):
    gid = r["gene_id"]
    src_fa = fasta_dir / f"{gid}.fasta"
    src_tr = tree_dir / f"{gid}.nwk"
    if not src_fa.exists() or not src_tr.exists():
        rows.append({"gene_id": gid, "status": "missing_input"})
        continue
    seqs = read_fasta(src_fa)
    L = {len(s) for s in seqs.values()}
    if len(L) != 1 or L.pop() % 3 != 0:
        rows.append({"gene_id": gid, "status": "bad_alignment"})
        continue
    ncod = r["aln_codons"]
    species = list(seqs)
    # gap mask per species: codon masked if any gap/N/? in real codon
    masks = {sp: [("?" if ("-" in s[j*3:j*3+3] or "N" in s[j*3:j*3+3] or "?" in s[j*3:j*3+3])
                   else ".") for j in range(int(ncod))] for sp, s in seqs.items()}
    # tree: strip {Test} annotations for evolver
    tree = src_tr.read_text().strip().replace("{Test}", "")
    # count leaves in tree
    nleaf = tree.count(",") + 1
    seed = 1000003 + 2 * i

    with tempfile.TemporaryDirectory(dir=work_dir) as td:
        td = Path(td)
        freq_block = "\n".join("  ".join(f"{v:.10f}" for v in freqs[k*4:k*4+4]) for k in range(16))
        ctl = (
            "0\n"
            f"{seed}\n"
            f"{nleaf} {int(ncod)} 1\n"
            "-1\n"
            f"{tree}\n"
            "1.0\n"
            "2.0\n"
            f"{freq_block}\n"
            "0\n"
        )
        (td / "MCcodon.dat").write_text(ctl)
        p = subprocess.run([evolver_bin, "6", "MCcodon.dat"], cwd=td,
                           capture_output=True, text=True, timeout=120)
        mc = td / "mc.txt"
        if not mc.exists() or p.returncode != 0:
            rows.append({"gene_id": gid, "status": f"evolver_fail_rc{p.returncode}"})
            continue
        sim = parse_mc(mc)
        # sanity
        if len(sim) != nleaf or any(len(s) != 3 * int(ncod) for s in sim.values()):
            rows.append({"gene_id": gid, "status": "mc_parse_mismatch",
                         "nseq": len(sim), "lens": sorted({len(s) for s in sim.values()})[:3]})
            continue
        # apply species set + gap mask of the real gene
        out_seqs = {}
        for sp in species:
            if sp not in sim:
                rows.append({"gene_id": gid, "status": f"species_missing_in_sim:{sp}"})
                out_seqs = None
                break
            s = sim[sp]
            chars = []
            for j in range(int(ncod)):
                cod = s[j*3:j*3+3]
                if masks[sp][j] == "?":
                    chars.append("---")
                else:
                    chars.append(cod)
            out_seqs[sp] = "".join(chars)
        if out_seqs is None:
            continue
        write_fasta(out_dir / f"{gid}.fasta", out_seqs)
        rows.append({"gene_id": gid, "status": "ok", "n_codons": int(ncod), "n_species": len(species)})
    if (i + 1) % 50 == 0:
        print(f"{i+1}/{len(genes)}", flush=True)

ok = [r for r in rows if r["status"] == "ok"]
report = {
    "n_genes": len(rows), "n_ok": len(ok),
    "status_counts": {s: sum(1 for r in rows if r["status"] == s) for s in {r["status"] for r in rows}},
    "params": {"model": "PAML evolver 4.10.10 codon (option 6), M0", "omega": 1.0, "kappa": 2.0,
               "codon_freqs": "F3x4 estimated from all 4,974 real alignments (stops zeroed, renormalized)",
               "tree": "per-gene species tree (same as BUSTED run, {Test} stripped), absolute branch lengths",
               "missingness": "real gene per-species per-codon gap mask applied to simulated sequences"},
}
print(json.dumps(report, indent=2))
# write sim manifest csv
with open(out_dir / "_sim_manifest.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["gene_id", "status", "n_codons", "n_species"])
    w.writeheader()
    for r in rows:
        w.writerow({"gene_id": r["gene_id"], "status": r["status"],
                    "n_codons": r.get("n_codons", ""), "n_species": r.get("n_species", "")})
json.dump(report, open(out_dir / "_sim_report.json", "w"), indent=2)
