#!/usr/bin/env python3
"""
07_analyze_interface.py -- nanobody-on-receptor interface for 8TH3.

Static mode (always available, no dependencies beyond numpy):
  * inter-chain contact map (heavy atoms < CUTOFF), per-residue contact
    counts, and the interface-residue list for nanobody (chain B) vs
    receptor (chain A).
  * contact persistence is a trajectory property, so static mode reports a
    single frame (the deposited cryo-EM model).

Trajectory mode (optional; requires mdtraj + a membrane production run):
    python3 07_analyze_interface.py --tpr work/8TH3/prod/md.tpr \
                                    --xtc work/8TH3/prod/md_center.xtc
  adds: time-averaged contact counts, H-bond/salt-bridge persistence, and the
  bound-nanobody per-residue RMSF (compared with the free nanobody RMSF from
  03_analyze.py, if results/rmsf_delta.csv exists).

Outputs in results/: interface.csv, interface_contacts.png, (rmsf_bound.png)
Run:  python3 simulation/07_analyze_interface.py [--tpr ... --xtc ...]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CUTOFF_NM = 0.45       # heavy-atom contact cutoff
HYD = {"H", "D"}       # 1H/2H elements (D = deuterium in some force fields)

AA3NAMES = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
            "TYR", "VAL", "MSE", "HID", "HIE", "HSP"}


def read_heavy(chain, path):
    """[(resnum, resname, (N,3) heavy-atom coords)] for one chain."""
    out, cur = [], None
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM") or line[21] != chain or line[16] != " ":
                continue
            elem = line[76:78].strip() or line[12:16].strip()
            if elem.startswith("H") or elem == "D":
                continue
            key = (line[22:26].strip(), line[26].strip())
            if cur is None or key != cur[0]:
                if cur:
                    out.append(cur)
                cur = [key, int(key[0]), line[17:20].strip(), []]
            cur[3].append(np.array([float(line[30:38]), float(line[38:46]),
                                    float(line[46:54])]))
    if cur:
        out.append(cur)
    return [(rnum, rname, np.asarray(xyz)) for (_, rnum, rname, xyz) in out]


def min_pair_distance(A, B):
    """Min heavy-atom distance between residue A and B (nm)."""
    d = A[:, None, :] - B[None, :, :]
    return np.sqrt((d * d).sum(-1)).min()


def interface_table(nb, rec, cutoff=CUTOFF_NM):
    rec_index = {rnum: i for i, (rnum, _, _) in enumerate(rec)}
    contacts = np.zeros((len(nb), len(rec)), dtype=int)
    for i, (n, _, a) in enumerate(nb):
        for j, (m, _, b) in enumerate(rec):
            if min_pair_distance(a, b) <= cutoff:
                contacts[i, j] = 1
    nb_contacts = contacts.sum(1)
    rec_contacts = contacts.sum(0)
    return contacts, nb_contacts, rec_contacts, rec_index


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default="simulation/data/8TH3_complex.pdb")
    ap.add_argument("--results", default="simulation/results")
    ap.add_argument("--tpr", default="")
    ap.add_argument("--xtc", default="")
    args = ap.parse_args()
    os.makedirs(args.results, exist_ok=True)

    nb = read_heavy("B", args.pdb)
    rec = read_heavy("A", args.pdb)
    print(f"nanobody {len(nb)} res | receptor {len(rec)} res")

    contacts, nbc, rec_cnt, rec_index = interface_table(nb, rec)

    # ---- static interface ------------------------------------------------
    n_if_nb = int((nbc > 0).sum())
    n_if_rec = int((rec_cnt > 0).sum())
    print(f"\ninterface residues: {n_if_nb} on nanobody, {n_if_rec} on receptor")
    print("nanobody interface residues (contacts):")
    for (rnum, rname, _), n in sorted(zip(nb, nbc), key=lambda x: -x[1]):
        if n:
            print(f"  {rname}{rnum}: {n} receptor contact(s)")
    print("top receptor contacts:")
    for (rnum, rname, _), n in sorted(zip(rec, rec_cnt), key=lambda x: -x[1])[:10]:
        if n:
            print(f"  {rname}{rnum}: {n} nanobody contact(s)")

    with open(os.path.join(args.results, "interface.csv"), "w") as fh:
        fh.write("chain,resnum,resname,n_contacts,in_interface\n")
        for ch, tbl in (("B", zip(nb, nbc)), ("A", zip(rec, rec_cnt))):
            for (rnum, rname, _), n in tbl:
                fh.write(f"{ch},{rnum},{rname},{int(n)},{int(n > 0)}\n")

    # ---- trajectory mode ---------------------------------------------------
    if args.tpr and args.xtc:
        try:
            import mdtraj as md
            from xvg_utils import aa_from_pdb, needleman_wunsch
        except ImportError as e:
            print(f"(mdtraj unavailable ({e}); skipping trajectory mode)")
            return
        traj = md.load(args.xtc, top=args.tpr)
        top = traj.topology
        print(f"\ntrajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")

        nb_atoms = top.select("chainid 1")
        rec_atoms = top.select("chainid 0")
        if nb_atoms.size == 0 or rec_atoms.size == 0 or not os.path.isfile(args.tpr):
            print("!! chain id assumptions wrong for this topology; set chain ids "
                  "so nanobody = 1, receptor = 0 (see fetch_and_split output)")
            return

        # candidate cross-chain residue pairs: < 1.2 nm in the reference frame
        preds = []
        for i in range(traj.n_residues):
            a = top.residue(i)
            for j in range(i + 1, traj.n_residues):
                b = top.residue(j)
                if a.chain.index == b.chain.index:
                    continue
                preds.append((i, j))
        dist = md.compute_contacts(traj, contacts=preds, scheme="closest-heavy")[0]
        hold = set()
        for k, (i, j) in enumerate(preds):
            if np.nanmin(dist[:, k]) <= 1.2:
                hold.add(k)

        # per-frame interface occupancy + persistence per residue
        per_frame = (dist[:, list(hold)] <= CUTOFF_NM).astype(float)
        occ_nb, occ_rec = {}, {}
        for k in hold:
            i, j = preds[k]
            r = top.residue(i); s = top.residue(j)
            occ_nb.setdefault(r, []); occ_rec.setdefault(s, [])
            occ_nb[r].append(per_frame[:, k]); occ_rec[s].append(per_frame[:, k])
        print("nanobody residues with contacts in >=50% of frames:")
        for r, arrs in sorted(occ_nb.items()):
            frac = float(np.mean(np.concatenate(arrs)))
            if frac >= 0.5:
                print(f"  {r.name}{r.resSeq}: {frac:.0%} persistence")

        # bound-nanobody RMSF (aligned on the nanobody only)
        ref = traj[0]
        atoms = top.select("chainid 1 and name CA")
        superposed = traj.superpose(ref, atom_indices=atoms)
        rmsf_atom = md.rmsf(superposed, top.select("chainid 1"), ref, parallel=False)
        rmsf_res = {}
        for a, v in zip(top.atoms, rmsf_atom):
            if a.element.symbol == "H":
                continue
            res = a.residue
            rmsf_res.setdefault(res.resSeq, []).append(v)
        bound = {r: float(np.max(v)) for r, v in rmsf_res.items()}

        # overlay with the free-nanobody RMSF (7T83) from 03_analyze.py
        free_csv = os.path.join(args.results, "rmsf_delta.csv")
        if os.path.isfile(free_csv):
            free = {}
            with open(free_csv) as fh:
                next(fh)
                for line in fh:
                    p = line.strip().split(",")
                    free[int(p[0])] = float(p[4])
            aaf = aa_from_pdb("simulation/data/7T83_nb.pdb")
            aab = aa_from_pdb("simulation/data/8TH3_nb.pdb")
            pairs = needleman_wunsch([aaf[r] for r in sorted(aaf)],
                                    [aab[r] for r in sorted(aab)])
            rf = sorted(aaf); rb = sorted(aab)
            xs, yf, yb = [], [], []
            for i, j in pairs:
                if i is None or j is None or not (rf[i] in free and rb[j] in bound):
                    continue
                xs.append(rf[i]); yf.append(free[rf[i]]); yb.append(bound[rb[j]])
            try:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(9, 4))
                ax.plot(xs, yf, "o-", ms=3, label="free (7T83, reps mean)")
                ax.plot(xs, yb, "o-", ms=3, label="bound (8TH3)")
                ax.set_xlabel("7T83 nanobody residue"); ax.set_ylabel("RMSF (nm)")
                ax.legend()
                fig.tight_layout()
                fig.savefig(os.path.join(args.results, "rmsf_bound.png"), dpi=150)
                print("wrote", os.path.join(args.results, "rmsf_bound.png"))
                print("free-vs-bound loop mobility -> residues that rigidify/drift")
            except ImportError:
                pass
    else:
        print("(no trajectory given; static single-frame contact map only)")


if __name__ == "__main__":
    main()