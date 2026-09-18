# -*- coding: utf-8 -*-
"""Phase9 R6 setup: stratified sampling + alignment stats + F3x4 codon frequencies.

Outputs (results/phase9_hardening/):
  - stratified_400_genes.csv : 400-gene stratified subset for rerun groups (Tasks 3-4)
  - sim_400_genes.csv        : 400 genes sampled by alignment-length deciles for neutral sim (Task 2)
  - alignment_quality_stats.json : pairwise identity / gap fraction on 500 sampled genes (Task 3d)
  - f3x4_codon_freqs.json    : codon frequencies (F3x4) estimated from all alignments
  - setup_report.log
"""
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd

random.seed(20260917)
np.random.seed(20260917)

BASE = Path(r"d:\人类正选择基因项目")
MASTER = BASE / "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv"
FASTA_DIR = BASE / "results/phase4_hyphy/fasta_v3_cleaned"
OUT_DIR = BASE / "results/phase9_hardening"
OUT_DIR.mkdir(parents=True, exist_ok=True)

log_lines = []
def log(s):
    print(s)
    log_lines.append(str(s))

def read_fasta(path):
    """Return dict name -> sequence (whitespace-stripped, uppercase)."""
    seqs = {}
    name = None
    chunks = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs

# ------------------------------------------------------------------
# 1. Load master table
# ------------------------------------------------------------------
df = pd.read_csv(MASTER)
log(f"master table: {len(df)} genes")
log(f"bh_fdr_sig: {int(df.bh_fdr_sig.sum())} significant, {int((~df.bh_fdr_sig.astype(bool)).sum())} non-significant")

# ------------------------------------------------------------------
# 2. Read all alignments once: length, gap stats, base composition
# ------------------------------------------------------------------
aln = {}
pos_counts = np.zeros((3, 4))  # codon position x (T,C,A,G)
IDX = {"T": 0, "C": 1, "A": 2, "G": 3}

for gid in df.gene_id:
    fp = FASTA_DIR / f"{gid}.fasta"
    if not fp.exists():
        log(f"WARN missing fasta {gid}")
        continue
    seqs = read_fasta(fp)
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        log(f"WARN ragged alignment {gid}: {lengths}")
    L = max(lengths)
    aln[gid] = {"L_nt": L, "L_codon": L // 3, "n_seq": len(seqs),
                "seqs": seqs}
    # base composition at codon positions (ungapped codons, all sequences)
    for s in seqs.values():
        for i in range(0, L - L % 3, 3):
            codon = s[i:i+3]
            if "-" in codon or "N" in codon or "?" in codon:
                continue
            for p, b in enumerate(codon):
                if b in IDX:
                    pos_counts[p, IDX[b]] += 1

log(f"read {len(aln)} alignments")

# F3x4 codon frequencies
pos_freq = pos_counts / pos_counts.sum(axis=1, keepdims=True)
log("pos1 T,C,A,G: " + " ".join(f"{x:.4f}" for x in pos_freq[0]))
log("pos2 T,C,A,G: " + " ".join(f"{x:.4f}" for x in pos_freq[1]))
log("pos3 T,C,A,G: " + " ".join(f"{x:.4f}" for x in pos_freq[2]))

BASES = "TCAG"
codon_freqs = []
for b1 in BASES:
    for b2 in BASES:
        for b3 in BASES:
            codon = b1 + b2 + b3
            if codon in ("TAA", "TAG", "TGA"):
                codon_freqs.append(0.0)
            else:
                codon_freqs.append(pos_freq[0, IDX[b1]] * pos_freq[1, IDX[b2]] * pos_freq[2, IDX[b3]])
tot = sum(codon_freqs)
codon_freqs = [x / tot for x in codon_freqs]
with open(OUT_DIR / "f3x4_codon_freqs.json", "w", encoding="utf-8") as f:
    json.dump({
        "codon_order": "TTT TTC TTA TTG TCT ... GGG (PAML MCcodon.dat fixed order)",
        "position_base_freqs": {"pos1": pos_freq[0].tolist(), "pos2": pos_freq[1].tolist(), "pos3": pos_freq[2].tolist()},
        "codon_freqs_f3x4_stops_zeroed_renorm": codon_freqs,
        "gc_content_overall": float((pos_counts[:, 1].sum() + pos_counts[:, 3].sum()) / pos_counts.sum()),
    }, f, indent=2)

# alignment length distribution
Lc = np.array([aln[g]["L_codon"] for g in aln])
log(f"alignment codon lengths: n={len(Lc)}, min={Lc.min()}, q25={np.percentile(Lc,25):.0f}, "
    f"median={np.median(Lc):.0f}, q75={np.percentile(Lc,75):.0f}, max={Lc.max()}")

# ------------------------------------------------------------------
# 3. Stratified 400-gene subset (Task 3): sig/nonsig x length tercile x LRT tercile
# ------------------------------------------------------------------
df = df[df.gene_id.isin(aln)].copy()
df["aln_codons"] = df.gene_id.map(lambda g: aln[g]["L_codon"])
sig = df[df.bh_fdr_sig.astype(bool)].copy()
nonsig = df[~df.bh_fdr_sig.astype(bool)].copy()
log(f"sig={len(sig)} nonsig={len(nonsig)}")

def tercile_labels(s):
    q1, q2 = np.percentile(s, [33.33, 66.67])
    return np.digitize(s, [q1, q2])  # 0,1,2

for d in (sig, nonsig):
    d["len_terc"] = tercile_labels(d.aln_codons.values)
    d["lrt_terc"] = tercile_labels(d.busted_lrt.fillna(0).values)

TARGET_PER_STRATUM = 200
subset_rows = []
rng = np.random.default_rng(20260917)
for d, label in ((sig, "sig"), (nonsig, "nonsig")):
    cells = d.groupby(["len_terc", "lrt_terc"])
    n_cells = len(cells)
    per_cell = TARGET_PER_STRATUM // n_cells
    remainder = TARGET_PER_STRATUM - per_cell * n_cells
    cell_sizes = cells.size().sort_values(ascending=False)
    alloc = {}
    for i, (key, size) in enumerate(cell_sizes.items()):
        alloc[key] = per_cell + (1 if i < remainder else 0)
    leftovers = []
    count = 0
    for key, idx in cells.groups.items():
        sub = d.loc[idx]
        take = min(alloc[key], len(sub))
        chosen = rng.choice(len(sub), size=take, replace=False)
        rest = np.setdiff1d(np.arange(len(sub)), chosen)
        subset_rows.append(sub.iloc[chosen])
        count += take
        if len(rest):
            leftovers.append(sub.iloc[rest])
    pool = pd.concat(leftovers) if leftovers else d.iloc[0:0]
    while count < TARGET_PER_STRATUM and len(pool):
        pick = rng.choice(len(pool), size=min(TARGET_PER_STRATUM - count, len(pool)), replace=False)
        subset_rows.append(pool.iloc[pick])
        count += len(pick)
        pool = pool.drop(pool.index[pick])
    log(f"stratum {label}: selected {count}")

subset = pd.concat(subset_rows)
assert len(subset) == 400, f"got {len(subset)}"
subset = subset.sample(frac=1, random_state=20260917).reset_index(drop=True)  # shuffle
subset_out = subset[["gene_id", "gene_symbol", "cds_length", "aln_codons", "busted_p",
                     "busted_lrt", "bh_fdr", "bh_fdr_sig", "len_terc", "lrt_terc", "classification_v7"]]
subset_out.to_csv(OUT_DIR / "stratified_400_genes.csv", index=False)
log(f"stratified_400_genes.csv: {len(subset)} genes; sig={int(subset.bh_fdr_sig.astype(bool).sum())}")
log("length tercile x significance counts:\n" +
    str(subset.groupby([subset.bh_fdr_sig.astype(bool), "len_terc"]).size()))

# ------------------------------------------------------------------
# 4. Simulation sample: 400 genes by alignment-length deciles (Task 2)
# ------------------------------------------------------------------
# exclude genes already in the rerun subset? No -- independent sample OK, but prefer
# disjoint to maximize information; sample from remaining genes.
remaining = df[~df.gene_id.isin(subset.gene_id)].copy()
remaining["len_decile"] = pd.qcut(remaining.aln_codons, 10, labels=False, duplicates="drop")
sim_rows = []
for dec, idx in remaining.groupby("len_decile").groups.items():
    take = min(40, len(idx))
    chosen = rng.choice(len(idx), size=take, replace=False)
    sim_rows.append(remaining.loc[idx[chosen]])
sim_sample = pd.concat(sim_rows)
while len(sim_sample) < 400 and len(remaining) > len(sim_sample):
    rest = remaining[~remaining.gene_id.isin(sim_sample.gene_id)]
    pick = rng.choice(len(rest), size=min(400 - len(sim_sample), len(rest)), replace=False)
    sim_sample = pd.concat([sim_sample, rest.iloc[pick]])
log(f"sim sample: {len(sim_sample)} genes")
sim_out = sim_sample[["gene_id", "gene_symbol", "aln_codons", "cds_length"]].copy()
sim_out["n_species"] = sim_sample.gene_id.map(lambda g: aln[g]["n_seq"])
sim_out.to_csv(OUT_DIR / "sim_400_genes.csv", index=False)

# ------------------------------------------------------------------
# 5. Alignment quality stats on 500 sampled genes (Task 3d)
# ------------------------------------------------------------------
qual_pool = df.sample(n=min(500, len(df)), random_state=20260917)
stats_rows = []
for gid in qual_pool.gene_id:
    seqs = aln[gid]["seqs"]
    names = list(seqs)
    L = aln[gid]["L_nt"]
    gap_fracs, idents = [], []
    for i in range(len(names)):
        si = seqs[names[i]]
        gap_fracs.append(si.count("-") / len(si))
        for j in range(i + 1, len(names)):
            sj = seqs[names[j]]
            match = tot = 0
            for a, b in zip(si, sj):
                if a in "-N??" or b in "-N??":
                    continue
                tot += 1
                if a == b:
                    match += 1
            if tot > 0:
                idents.append(match / tot)
    stats_rows.append({
        "gene_id": gid,
        "aln_codons": aln[gid]["L_codon"],
        "n_seq": len(names),
        "mean_gap_frac": float(np.mean(gap_fracs)),
        "mean_pairwise_identity": float(np.mean(idents)) if idents else np.nan,
    })
qdf = pd.DataFrame(stats_rows)
qdf.to_csv(OUT_DIR / "alignment_quality_500.csv", index=False)
qual_summary = {
    "n_genes_sampled": int(len(qdf)),
    "mean_pairwise_identity": {
        "mean": float(qdf.mean_pairwise_identity.mean()),
        "median": float(qdf.mean_pairwise_identity.median()),
        "p5": float(qdf.mean_pairwise_identity.quantile(0.05)),
        "p95": float(qdf.mean_pairwise_identity.quantile(0.95)),
    },
    "gap_fraction_per_sequence": {
        "mean": float(qdf.mean_gap_frac.mean()),
        "median": float(qdf.mean_gap_frac.median()),
        "p5": float(qdf.mean_gap_frac.quantile(0.05)),
        "p95": float(qdf.mean_gap_frac.quantile(0.95)),
    },
    "n_species_per_alignment": {
        "mean": float(qdf.n_seq.mean()),
        "min": int(qdf.n_seq.min()),
        "max": int(qdf.n_seq.max()),
    },
}
with open(OUT_DIR / "alignment_quality_stats.json", "w", encoding="utf-8") as f:
    json.dump(qual_summary, f, indent=2, ensure_ascii=False)
log("alignment quality: " + json.dumps(qual_summary, indent=1))

with open(OUT_DIR / "setup_report.log", "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines))
log("DONE")
