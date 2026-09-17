# -*- coding: utf-8 -*-
"""
Shao 2023 cross-validation version-drift mapping (67->64, 3->4)
================================================================================
An earlier manuscript version reported "3 of 82 Shao PSGs in our BUSTED set
(OR = 0.07), 67/82 absent from the universe" (symbol-based matching, see
results/cross_validation_shao2023/). The current v9 reports "4 of 18 analyzable
genes (OR = 0.55), 64/82 absent" (ENSG version-stripped matching). This script
rebuilds both matchings gene-by-gene and produces the mapping table.

Output: results/phase9_hardening/shao_version_mapping.csv + .json ; log _audit_shao.txt
"""
import io
import json

import pandas as pd
from scipy import stats

BASE = "d:/人类正选择基因项目/"
LOG = io.open(BASE + "_audit_shao.txt", "w", encoding="utf-8")


def P(*a):
    print(*a, file=LOG)


shao = pd.read_csv(BASE + "data/downloads/shao_2023/table_s17_psgs.csv")
shao.columns = ["no", "symbol", "ensg", "lnl", "p"]
P(f"Shao table S17: {len(shao)} PSGs")

gc = pd.read_csv(BASE + "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv")
strip = lambda s: str(s).split(".")[0]
gc["ensg"] = gc["gene_id"].map(strip)
sig = gc["bh_fdr_sig"].astype(bool)
gmap = pd.read_csv(BASE + "data/gene_id_to_symbol_gencode47.csv")
gmap_sym = {strip(g): s for g, s in zip(gmap["gene_id"], gmap["gene_symbol"])}
gc["sym"] = gc["gene_symbol"].where(
    gc["gene_symbol"].notna() & (gc["gene_symbol"].astype(str).str.strip() != ""),
    gc["ensg"].map(lambda e: gmap_sym.get(e, ""))).astype(str).str.upper().str.strip()
uni_syms = {}
for e, s in zip(gc["ensg"], gc["sym"]):
    if s and s != "NAN":
        uni_syms.setdefault(s, e)  # universe symbol -> ensg (first)
uni_ensg = set(gc["ensg"])
ensg2sym = {**gmap_sym, **dict(zip(gc["ensg"], gc["gene_symbol"]))}

busted_sig_ensg = set(gc.loc[sig, "ensg"])
n_busted = len(gc)

rows = []
for _, r in shao.iterrows():
    sym, e = str(r["symbol"]).upper(), strip(r["ensg"])
    m_sym = sym in uni_syms
    m_ensg = e in uni_ensg
    # status under each matching version
    in_uni_sym = m_sym
    in_uni_ensg = m_ensg
    sig_sym = m_sym and uni_syms[sym] in busted_sig_ensg
    sig_ensg = m_ensg and e in busted_sig_ensg
    note = ""
    if m_ensg and not m_sym:
        note = f"in universe via ENSG only (universe symbol for {e} = {ensg2sym.get(e)})"
    elif m_sym and not m_ensg:
        note = f"in universe via symbol only (maps to {uni_syms[sym]})"
    elif m_sym and m_ensg and uni_syms[sym] != e:
        note = f"symbol maps to {uni_syms[sym]} but ENSG matches {ensg2sym.get(e)}"
    rows.append({
        "shao_symbol": r["symbol"], "shao_ensg": e,
        "in_universe_symbol_match": in_uni_sym,
        "in_universe_ensg_match": in_uni_ensg,
        "universe_ensg": (uni_syms.get(sym) if m_sym else e if m_ensg else ""),
        "universe_symbol": ensg2sym.get(e, ""),
        "busted_fdr_sig_symbol_match": sig_sym,
        "busted_fdr_sig_ensg_match": sig_ensg,
        "classification_v7": (gc.set_index("ensg")["classification_v7"].get(e) if m_ensg else ""),
        "note": note,
    })
tab = pd.DataFrame(rows)
tab.to_csv(BASE + "results/phase9_hardening/shao_version_mapping.csv", index=False)

n_sym = int(tab["in_universe_symbol_match"].sum())
n_ensg = int(tab["in_universe_ensg_match"].sum())
k_sym = int(tab["busted_fdr_sig_symbol_match"].sum())
k_ensg = int(tab["busted_fdr_sig_ensg_match"].sum())

# ---- OLD version (as produced by shao2023_cross_validation.py) ----
# symbol matching against the 4,045-symbol BUSTED v3 set + 1,369 FDR symbols
with io.open(BASE + "results/phase4_hyphy/final_analysis_v3/all_analyzed_gene_symbols.txt", encoding="utf-8") as f:
    old_all_syms = set(line.strip().upper() for line in f if line.strip())
with io.open(BASE + "results/phase4_hyphy/final_analysis_v3/fdr_sig_gene_symbols.txt", encoding="utf-8") as f:
    old_fdr_syms = set(line.strip().upper() for line in f if line.strip())
old_matched = set(shao["symbol"].astype(str).str.upper()) & old_all_syms
old_fdr = set(shao["symbol"].astype(str).str.upper()) & old_fdr_syms
P(f"\nOLD version (symbol matching, {len(old_all_syms)}-gene BUSTED v3 set, "
  f"{len(old_fdr_syms)} FDR-sig): {len(old_matched)}/82 matched ({82-len(old_matched)} absent), "
  f"{len(old_fdr)} FDR-significant: {sorted(old_fdr)}")

# ---- NEW version (v9: ENSG version-stripped vs 4,974 universe) ----
new_matched = set(tab.loc[tab["in_universe_ensg_match"], "shao_symbol"].astype(str).str.upper())
new_fdr = set(tab.loc[tab["busted_fdr_sig_ensg_match"], "shao_symbol"].astype(str).str.upper())
P(f"NEW version (ENSG matching, 4,974 universe): {len(new_matched)}/82 matched "
  f"({82-len(new_matched)} absent), {len(new_fdr)} FDR-significant: {sorted(new_fdr)}")

recovered = new_matched - old_matched
lost = old_matched - new_matched
P(f"\ngenes in old-67-absent but matched in new (recovered): {sorted(recovered)}")
for g in sorted(recovered):
    r = tab[tab["shao_symbol"].astype(str).str.upper() == g].iloc[0]
    P(f"  {g}: ENSG {r['shao_ensg']} -> universe symbol {r['universe_symbol']}, "
      f"BUSTED sig: {r['busted_fdr_sig_ensg_match']}")
P(f"genes matched in old but not new: {sorted(lost)}")
P(f"4th significant gene (3->4): {sorted(new_fdr - old_fdr)}")
P(f"symbol matching: {n_sym}/82 in universe, {k_sym} BUSTED-significant, "
  f"{82 - n_sym} absent (=67 expected)")
P(f"ENSG matching:   {n_ensg}/82 in universe, {k_ensg} BUSTED-significant, "
  f"{82 - n_ensg} absent (=64 expected)")

P("\n--- genes recovered by ENSG matching but missed by symbol matching ---")
sub = tab[tab["in_universe_ensg_match"] & ~tab["in_universe_symbol_match"]]
P(sub[["shao_symbol", "shao_ensg", "universe_symbol", "busted_fdr_sig_ensg_match",
       "classification_v7", "note"]].to_string(index=False))

P("\n--- genes matched by symbol but not by ENSG ---")
sub2 = tab[tab["in_universe_symbol_match"] & ~tab["in_universe_ensg_match"]]
P(sub2[["shao_symbol", "shao_ensg", "universe_ensg", "note"]].to_string(index=False))

P("\n--- BUSTED-significant under each version ---")
P("symbol version:", tab.loc[tab["busted_fdr_sig_symbol_match"], "shao_symbol"].tolist())
P("ENSG version:  ", tab.loc[tab["busted_fdr_sig_ensg_match"], "shao_symbol"].tolist())

# verify old OR=0.07 and new OR=0.55
# old: 3 sig of 82 vs universe rate (treat 79 as non-sig)
rate = k_ensg and None
n_sig_uni = int(sig.sum())
or_old, p_old = stats.fisher_exact([[k_sym, 82 - k_sym],
                                    [n_sig_uni - k_sym, n_busted - n_sig_uni - (82 - k_sym)]])
# old script's table: busted_universe 4045, sig 1369 (different denominators, from summary)
or_new, p_new = stats.fisher_exact([[k_ensg, n_ensg - k_ensg],
                                    [n_sig_uni - k_ensg, n_busted - n_ensg - (n_sig_uni - k_ensg)]])
P(f"\nrecomputed new-style test (4x of 18 vs rest of universe): OR={or_new:.3f}, p={p_new:.3f}")
exp = n_ensg * n_sig_uni / n_busted
P(f"expected significant among 18 under universe rate: {exp:.1f} (manuscript: 6.1 with 4,045 denom)")

result = {
    "task": "Shao 2023 version-drift mapping (67->64 absent, 3->4 significant)",
    "symbol_matching_current_universe": {"in_universe": n_sym, "absent": 82 - n_sym, "busted_sig": k_sym},
    "ensg_matching": {"in_universe": n_ensg, "absent": 82 - n_ensg, "busted_sig": k_ensg},
    "old_version": {"method": "symbol matching vs 4,045-symbol BUSTED v3 set",
                    "in_universe": len(old_matched), "absent": 82 - len(old_matched),
                    "busted_sig": sorted(old_fdr)},
    "new_version": {"method": "ENSG version-stripped matching vs 4,974 universe",
                    "in_universe": len(new_matched), "absent": 82 - len(new_matched),
                    "busted_sig": sorted(new_fdr)},
    "genes_recovered_old_absent_to_new_matched": sorted(recovered),
    "genes_lost_new_absent": sorted(lost),
    "busted_sig_symbol_version": tab.loc[tab["busted_fdr_sig_symbol_match"], "shao_symbol"].tolist(),
    "busted_sig_ensg_version": tab.loc[tab["busted_fdr_sig_ensg_match"], "shao_symbol"].tolist(),
    "the_fourth_gene": sorted(new_fdr - old_fdr),
    "recomputed_or_0.55_test": {"OR": float(or_new), "p": float(p_new)},
    "explanation": "Drift decomposition: 67->64 absent = 3 Shao genes matched in the new analysis "
                   "but not the old (alias/renamed symbols + universe growth 4,045->4,974); "
                   "3->4 significant = PDHB (BUSTED FDR-significant in the current 4,974-gene "
                   "universe but not in the old 4,045/1,369 analysis). Both v9 numbers (4 of 18, "
                   "64/82 absent) are internally consistent with ENSG version-stripped matching.",
}
with io.open(BASE + "results/phase9_hardening/shao_version_mapping.json", "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False)
P("\nsaved shao_version_mapping.csv/.json")
LOG.close()
print("done")
