"""R8 optional Task 2: does the codon-aware realignment + filtering inflate
detection under the NULL? Realign + BUSTED of the 400 omega=1 neutral sims.

Compares BH FDR<0.05 rates: original sim alignments (R7 Task 2 parametric
null, busted_rerun_sim2) vs codon-aware realigned sim alignments
(results_sim_realign). Under proper calibration both should be near zero;
a large positive delta under null would mean the +14.8pp observed on real
data is a realignment artifact.

Output: neutral_realign_calibration_r8.json
"""
import glob
import json
import os
import statistics as st

PROJ = r"d:\人类正选择基因项目"
HARD = os.path.join(PROJ, r"results\phase9_hardening")


def read_p(json_dir):
    out = {}
    for f in glob.glob(os.path.join(json_dir, "*.json")):
        gene = os.path.basename(f)[:-5]
        try:
            if os.path.getsize(f) < 50:
                continue
            d = json.load(open(f))
            tr = d.get("test results")
            if not tr:
                continue
            out[gene] = float(tr["p-value"])
        except Exception:
            continue
    return out


def bh_q(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q_sorted = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        q = min(prev, pvals[i] * n / rank)
        q_sorted[rank - 1] = q
        prev = q
    q = [0.0] * n
    for rank in range(n):
        q[order[rank]] = q_sorted[rank]
    return q


def main():
    # NOTE: the realigned BUSTED master runs over ~800 alignments because the
    # realign output dir accumulated BOTH sim sets (R6 fixed-param fasta_sim and
    # R7 parametric fasta_sim2). Matched-subset comparisons are therefore done
    # per set: primary = R7 parametric (pre-registered), replication = R6.
    sim2_p = read_p(os.path.join(HARD, "busted_rerun_sim2"))   # R7 parametric null
    sim6_p = read_p(os.path.join(HARD, "busted_rerun_sim"))    # R6 fixed-param null
    re_p = read_p(os.path.join(HARD, "busted_rerun_sim_realign"))
    if not re_p:
        print("no sim_realign results yet")
        return

    def stats(pmap, alpha=0.05):
        genes = sorted(pmap)
        ps = [pmap[g] for g in genes]
        qs = bh_q(ps)
        n = len(ps)
        return {
            "n": n,
            "p_lt_0.05": sum(1 for p in ps if p < 0.05) / n,
            "p_eq_0.5": sum(1 for p in ps if p == 0.5) / n,
            "bh_fdr": sum(1 for q in qs if q < alpha) / n,
            "n_rejections": sum(1 for q in qs if q < alpha),
        }

    def matched(base_p, tag):
        # restrict realigned results to genes present in the given baseline
        sub = {g: p for g, p in re_p.items() if g in base_p}
        common = sorted(set(base_p) & set(sub))
        s_b = stats({g: base_p[g] for g in common})
        s_r = stats({g: sub[g] for g in common})
        return {
            "set": tag,
            "n_common": len(common),
            "baseline": s_b,
            "realigned": s_r,
            "delta_p_lt_0.05_pp": round((s_r["p_lt_0.05"] - s_b["p_lt_0.05"]) * 100, 2),
            "delta_bh_fdr_pp": round((s_r["bh_fdr"] - s_b["bh_fdr"]) * 100, 2),
        }

    primary = matched(sim2_p, "R7_parametric_null_400")
    replication = matched(sim6_p, "R6_fixed_param_null_400")
    comb_base = {g: p for g, p in {**sim6_p, **sim2_p}.items()}
    combined = matched(comb_base, "combined_800")

    out = {
        "description": "R8 optional Task 2: matched neutral calibration for the +14.8pp "
                       "detection-rate increase observed after codon-aware realignment of the "
                       "stratified 400 genes. The same realignment+filtering pipeline applied "
                       "to omega=1 neutral simulations; BUSTED srv=No with per-gene trees; all "
                       "alignments null by construction. Primary comparison uses the R7 Task 2 "
                       "parametric-bootstrap set (per-gene null MLE params); the R6 fixed-param "
                       "set serves as an independent replication (the realign output directory "
                       "accumulated both sets, and each is compared only to its OWN matched "
                       "baseline). If realignment inflated detection under null, the artifact "
                       "interpretation would hold; if rates stay at the conservative null "
                       "level, the +14.8pp reflects genuine power gain.",
        "primary_r7_parametric": primary,
        "replication_r6_fixed": replication,
        "combined": combined,
        "n_realigned_total": len(re_p),
        "reference_real_data": "stratified 400: baseline 49.4% -> realigned 64.1% (+14.8pp)",
    }
    with open(os.path.join(HARD, "neutral_realign_calibration_r8.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
