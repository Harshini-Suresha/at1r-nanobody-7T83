#!/usr/bin/env python3
"""
03_analyze.py -- compare the 7T83 (WT) and 7T84 (G26D T57I) replicates.

Questions answered (with replicate-level uncertainty):
  1. Loop mobility: per-residue RMSF, ΔRMSF(7T84 - 7T83) with a Welch t stat,
     flagged where |t|>2 and |Δ|>0.05 nm.
  2. Replicate spread: across-replicate mean/SD/SE of RMSD, Rg, H-bonds and
     per-residue RMSF for each system (the "sampling error" of the run).
  3. Domain view: mean RMSF per CDR loop / framework, mapped 7T83 <-> 7T84
     by global sequence alignment (their residue numbers are not comparable
     directly -- 7T84 chain A is 117 aa vs 120 in 7T83).

Inputs : $WORK/<7T83|7T84>/rep*_seed*/{rmsd,rmsf_bb,gyrate,hbond}.xvg
Outputs: $RESULTS/*.csv and *.png
Run    : python3 simulation/03_analyze.py
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xvg_utils as U  # noqa: E402

SYSTEMS = ("7T83", "7T84")
CDR_MOTIFS = {"CDR1": "GYIYRRYR", "CDR2": "ISGGSSTN", "CDR3": "RIVSDPRVY"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="simulation/work")
    ap.add_argument("--results", default="simulation/results")
    ap.add_argument("--data", default="simulation/data")
    ap.add_argument("--dmin", type=float, default=0.05,
                    help="min |dRMSF| (nm) to flag a residue")
    ap.add_argument("--tmin", type=float, default=2.0,
                    help="min |t| to flag a residue")
    args = ap.parse_args()
    os.makedirs(args.results, exist_ok=True)

    # ---- gather per-replicate data ----------------------------------------
    rmsf = {s: {} for s in SYSTEMS}    # system -> {resid: np.array over reps}
    reps = {s: U.find_reps(args.work, s) for s in SYSTEMS}
    rmsd_s, rg_s, hb_s = {s: [] for s in SYSTEMS}, {s: [] for s in SYSTEMS}, \
                         {s: [] for s in SYSTEMS}
    for s in SYSTEMS:
        if not reps[s]:
            print(f"!! {s}: no replicate dirs under {args.work}/{s} -- run "
                  f"02_run_replicates.sh first (or check $WORK).")
        for d in reps[s]:
            t, rmsd = U.load_xvg(os.path.join(d, "rmsd.xvg"))
            resids, rmsf_vals = U.load_rmsf(os.path.join(d, "rmsf_bb.xvg"))
            _, rg = U.load_xvg(os.path.join(d, "gyrate.xvg"))
            _, hb = U.load_xvg(os.path.join(d, "hbond.xvg"))
            for res, val in zip(resids, rmsf_vals):
                rmsf[s].setdefault(int(round(res)), []).append(float(val))
            rmsd_s[s].append(float(rmsd.mean()))
            rg_s[s].append(float(rg.mean()))
            hb_s[s].append(float(hb.mean()))

    # ---- 7T83 <-> 7T84 residue map (alignment) -----------------------------
    aa83 = U.aa_from_pdb(os.path.join(args.data, "7T83_nb.pdb"))
    aa84 = U.aa_from_pdb(os.path.join(args.data, "7T84_nb.pdb"))
    r83s = sorted(aa83); r84s = sorted(aa84)
    pairs = U.needleman_wunsch([aa83[r] for r in r83s],
                               [aa84[r] for r in r84s])
    align = [(r83s[i], r84s[j]) for i, j in pairs if i is not None and j is not None]

    # locates CDR residues on the 7T83 reference numbering
    cdr83 = {name: set() for name in CDR_MOTIFS}
    seq83 = "".join(aa83[r] for r in r83s)
    for name, motif in CDR_MOTIFS.items():
        k = seq83.find(motif)
        if k >= 0:
            cdr83[name] = set(r83s[k:k + len(motif)])

    # ---- per-residue statistics --------------------------------------------
    rows = []
    for r83, r84 in align:
        v83 = rmsf["7T83"].get(r83)
        v84 = rmsf["7T84"].get(r84)
        if not v83 or not v84:
            continue
        m83, _, se83 = U.runstats(v83)
        m84, _, se84 = U.runstats(v84)
        if not np.isfinite(m83) or not np.isfinite(m84):
            continue
        delta = m84 - m83
        t = U.welch_t(delta, se83, se84)
        region = next((name for name, s in cdr83.items() if r83 in s), "FR/other")
        rows.append((r83, aa83[r83], r84, aa84[r84],
                     m83, se83, m84, se84, delta, t,
                     abs(t) >= args.tmin and abs(delta) >= args.dmin, region))
    df = np.asarray([(r[0], r[2], r[4], r[5], r[6], r[7], r[8], r[9], int(r[10])) for r in rows])

    # ---- loop-level mobility (CDR means per variant) ------------------------
    loop_rows = []
    for name, res83 in cdr83.items():
        res84 = {r84 for r83, r84 in align if r83 in res83}
        v83 = [np.mean(rmsf["7T83"][r]) for r in sorted(res83) if rmsf["7T83"].get(r)]
        v84 = [np.mean(rmsf["7T84"][r]) for r in sorted(res84) if rmsf["7T84"].get(r)]
        if v83 and v84:
            loop_rows.append((name, len(res83) & len(res84), np.mean(v83), np.mean(v84),
                              np.mean(v84) - np.mean(v83)))

    # ---- write outputs ------------------------------------------------------
    with open(os.path.join(args.results, "rmsf_delta.csv"), "w") as fh:
        fh.write("resid_7T83,aa_7T83,resid_7T84,aa_7T84,rmsf83_mean_nm,"
                 "rmsf83_se_nm,rmsf84_mean_nm,rmsf84_se_nm,delta_nm,t_stat,"
                 "significant,region\n")
        for r in rows:
            fh.write(",".join(map(str, r)) + "\n")

    with open(os.path.join(args.results, "replicate_summary.csv"), "w") as fh:
        fh.write("system,replicate,rmsd_mean_nm,rg_mean_nm,hbond_mean\n")
        for s in SYSTEMS:
            for i, d in enumerate(reps[s]):
                fh.write(f"{s},rep{i+1},{rmsd_s[s][i]:.4f},{rg_s[s][i]:.4f},{hb_s[s][i]:.1f}\n")

    with open(os.path.join(args.results, "loop_summary.csv"), "w") as fh:
        fh.write("loop,n_common_residues,rmsf83_mean_nm,rmsf84_mean_nm,"
                 "delta_nm\n")
        for row in loop_rows:
            fh.write(",".join(map(str, row)) + "\n")

    # ---- console report ----------------------------------------------------
    print("\n=== replicate spread (mean +/- SE over replicates) ===")
    print(f"{'system':<6}{'RMSD (nm)':>16}{'Rg (nm)':>12}{'H-bonds':>12}")
    for s in SYSTEMS:
        if rmsd_s[s]:
            _, _, se_r = U.runstats(rmsd_s[s]); _, _, se_g = U.runstats(rg_s[s])
            _, _, se_h = U.runstats(hb_s[s])
            print(f"{s:<6}{np.mean(rmsd_s[s]):>10.4f}+/-{float(se_r):.4f}"
                  f"{np.mean(rg_s[s]):>9.3f}+/-{float(se_g):.3f}"
                  f"{np.mean(hb_s[s]):>9.1f}+/-{float(se_h):.1f}")
    print("\n=== loop mobility ===")
    for row in loop_rows:
        print(f"{row[0]}: 7T83 {row[2]:.3f} nm | 7T84 {row[3]:.3f} nm | "
              f"delta {row[4]:+.3f} nm")
    sig = [r for r in rows if r[10]]
    print(f"\n=== significant CDR/loop residues (|t|>={args.tmin}, "
          f"|delta|>={args.dmin} nm) ===")
    for r in sorted([r for r in sig], key=lambda x: -abs(x[8])):
        print(f"  {r[0]}({r[1]}) -> {r[2]}({r[3]})  d={r[8]:+.3f} nm  t={r[9]:+.1f}  [{r[11]}]")
    if not sig:
        print("  none")

    # ---- figures ------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("(matplotlib not available; skipping figures)")
        return

    r83a = [r[0] for r in rows]
    delta = np.array([r[8] for r in rows])
    t = np.array([r[9] for r in rows])
    sig_flag = np.array([r[10] for r in rows])
    fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    ax[0].bar(r83a, delta, color=["#c0392b" if f else "#bdc3c7" for f in sig_flag])
    ax[0].axhline(0, color="k", lw=0.5)
    ax[0].set_ylabel("dRMSF 7T84-7T83 (nm)")
    ax[0].set_title("Per-residue mobility change, G26D + T57I vs WT "
                    "(red = |t|>2 and |d|>0.05 nm)")
    ax[1].plot(r83a, t, "o-", ms=3)
    ax[1].axhline(args.tmin, c="gray", ls="--"); ax[1].axhline(-args.tmin, c="gray", ls="--")
    ax[1].set_ylabel("Welch t"); ax[1].set_xlabel("7T83 residue number")
    for name, s in cdr83.items():
        if s:
            x0, x1 = min(s), max(s)
            ax[0].axvspan(x0 - 0.5, x1 + 0.5, alpha=0.08, color="#2980b9")
    fig.tight_layout()
    fig.savefig(os.path.join(args.results, "rmsf_delta.png"), dpi=150)
    print("wrote", os.path.join(args.results, "rmsf_delta.png"))

    fig, ax = plt.subplots(figsize=(8, 4))
    for s, c in zip(SYSTEMS, ("#2980b9", "#27ae60")):
        xs = rmsf[s]
        if not xs:
            continue
        r = sorted(xs)
        m = np.array([np.mean(xs[res]) for res in r])
        sd = np.array([np.std(xs[res]) if len(xs[res]) > 1 else 0 for res in r])
        ax.plot(r, m, label=s, color=c)
        ax.fill_between(r, m - sd, m + sd, color=c, alpha=0.18)
    ax.set_xlabel("residue"); ax.set_ylabel("RMSF (nm)"); ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.results, "replicate_spread.png"), dpi=150)
    print("wrote", os.path.join(args.results, "replicate_spread.png"))


if __name__ == "__main__":
    main()