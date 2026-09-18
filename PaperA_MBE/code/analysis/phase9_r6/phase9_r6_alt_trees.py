#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate alternative-topology trees ((human,gorilla) sister) for the 400-gene subset.

Transforms the local arrangement ((chimp,bonobo),human),gorilla into
((human,gorilla),(chimp,bonobo)), i.e. swaps gorilla with the human-sister clade,
keeping human's {Test} annotation and approximate branch lengths.

Usage (WSL): python3 phase9_r6_alt_trees.py <genes_csv> <tree_dir> <out_dir> <report_json>
"""
import csv
import json
import sys
from pathlib import Path

genes_csv, tree_dir, out_dir, report_json = sys.argv[1:5]
tree_dir = Path(tree_dir)
out_dir = Path(out_dir)
out_dir.mkdir(parents=True, exist_ok=True)


class Node:
    def __init__(self):
        self.children = []
        self.label = ""      # taxon label incl. {} annotation, or "" for internal
        self.length = None   # float or None

    def is_leaf(self):
        return not self.children

    def leaf_labels(self):
        if self.is_leaf():
            return {self.label.split("{")[0].strip()}
        s = set()
        for c in self.children:
            s |= c.leaf_labels()
        return s


def parse_newick(s):
    s = s.strip().rstrip(";").strip()
    pos = [0]

    def parse_node():
        node = Node()
        if s[pos[0]] == "(":
            pos[0] += 1
            while True:
                node.children.append(parse_node())
                if s[pos[0]] == ",":
                    pos[0] += 1
                    continue
                if s[pos[0]] == ")":
                    pos[0] += 1
                    break
                raise ValueError(f"bad newick at {pos[0]}: {s[pos[0]-20:pos[0]+20]}")
        # label
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ",():;":
            pos[0] += 1
        node.label = s[start:pos[0]].strip()
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            start = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ",()":
                pos[0] += 1
            node.length = float(s[start:pos[0]])
        return node

    tree = parse_node()
    return tree


def fmt(node):
    if node.is_leaf():
        out = node.label
    else:
        out = "(" + ",".join(fmt(c) for c in node.children) + ")" + node.label
    if node.length is not None:
        out += f":{node.length:g}"
    return out


def find_lca(node, a, b):
    """Return the deepest node whose leaves contain both a and b, plus path info."""
    leaves = node.leaf_labels()
    if a in leaves and b in leaves:
        best = None
        for c in node.children:
            r = find_lca(c, a, b)
            if r is not None:
                best = r
                break
        return ("child", best) if best else ("self", node)
    return None


def lca_node(root, a, b):
    res = find_lca(root, a, b)
    if res is None:
        return None
    while res[0] == "child":
        res = res[1]
    return res[1]


def transform(tree_str):
    """Swap gorilla with the human-sister clade. Returns (new_tree_str, note)."""
    root = parse_newick(tree_str)
    n = lca_node(root, "human", "gorilla")
    if n is None:
        return None, "human/gorilla not both present"
    hum_side = gor_side = None
    for c in n.children:
        lv = c.leaf_labels()
        if "human" in lv:
            hum_side = c
        elif "gorilla" in lv:
            gor_side = c
    if hum_side is None or gor_side is None:
        return None, "LCA children unexpected"
    # human must be a direct child of hum_side
    direct = [c for c in hum_side.children if c.is_leaf() and c.label.split("{")[0] == "human"]
    if not direct:
        return None, "human not direct child of its LCA-side (nested)"
    human = direct[0]
    sibs = [c for c in hum_side.children if c is not human]
    if gor_side.is_leaf() and gor_side.label.split("{")[0] != "gorilla":
        return None, "gorilla side is a non-gorilla leaf"
    # build new clade: ((human, gorilla_side), siblings...)
    new_pair = Node()
    new_pair.children = [human, gor_side]
    new_pair.length = hum_side.length  # keep A's length for the new pair
    new_n = Node()
    new_n.children = [new_pair] + sibs
    new_n.length = n.length
    # splice into root in place of n
    replaced = [False]

    def splice(node):
        if replaced[0]:
            return
        for i, c in enumerate(node.children):
            if c is n:
                node.children[i] = new_n
                replaced[0] = True
                return
            splice(c)

    if n is root:
        root = new_n
    else:
        splice(root)
        if not replaced[0]:
            return None, "splice failed"
    return fmt(root) + ";", "ok"


rows = []
for r in csv.DictReader(open(genes_csv, encoding="utf-8")):
    gid = r["gene_id"]
    src = tree_dir / f"{gid}.nwk"
    if not src.exists():
        rows.append({"gene_id": gid, "status": "missing_tree"})
        continue
    ts = src.read_text().strip()
    try:
        new_ts, note = transform(ts)
    except Exception as e:
        new_ts, note = None, f"error: {e}"
    if new_ts is None:
        rows.append({"gene_id": gid, "status": note, "orig": ts})
    else:
        (out_dir / f"{gid}.nwk").write_text(new_ts + "\n")
        rows.append({"gene_id": gid, "status": "ok", "orig": ts, "alt": new_ts})

ok = [r for r in rows if r["status"] == "ok"]
report = {
    "n_genes": len(rows),
    "n_ok": len(ok),
    "status_counts": {s: sum(1 for r in rows if r["status"] == s) for s in {r["status"] for r in rows}},
    "example_original": rows[0].get("orig", ""),
    "example_alt": next((r["alt"] for r in rows if r["status"] == "ok"), ""),
}
Path(report_json).write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps({k: report[k] for k in ("n_genes", "n_ok", "status_counts")}, indent=2))
for r in rows:
    if r["status"] != "ok":
        print("NON-OK:", r["gene_id"], r["status"])
