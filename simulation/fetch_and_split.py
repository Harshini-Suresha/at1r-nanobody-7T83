#!/usr/bin/env python3
"""
fetch_and_split.py -- download and curate the three structures used in the
nanobody / AT1R simulation suite.

Outputs (into simulation/data):
  * 7T83_nb.pdb, 7T84_nb.pdb        chain A of the free nanobody, no ligand/water
  * 8TH3_nb.pdb                     nanobody (chain B) from the cryo-EM complex
  * 8TH3_at1r.pdb                   AT1R receptor (chain A) from the complex
  * 8TH3_complex.pdb                nanobody + receptor, interface restored
                                    (rigid refit; see 8TH3_complex.meta.json)

What the PDB actually contains (verified against the coordinates):
  * 7T83 / 7T84  -- eight crystallographic copies of the free nanobody; we keep
    one.  Residue numbering is renormalised to 1..N for pdb2gmx.
  * 8TH3 (cryo-EM, 3.0 A) -- the deposited model splits the nanobody-receptor
    fusion into separate chains: A = AT1R receptor (393 aa, with BRIL/b562
    insert), B = AT118-H nanobody (99 aa), C/D = BAG2 anti-BRIL Fab.  The entry
    is the *locally refined* map (EMD-41248) and the nanobody is ~2.4 nm off the
    receptor (paratope faces the ECL epitope), so by default the nanobody is
    rigidly translated onto the receptor until heavy-atom contacts form.  Use
    --no-rebuild to keep the deposited (gapped) pose, and --strip-bril to carve
    the BRIL insert out (NOTE: creates non-native junction scars).

Run:  python3 fetch_and_split.py [--pdb-dir simulation/data] [--strip-bril]
"""
import argparse
import json
import os
import sys
import urllib.request
import ssl

import numpy as np

RCSB = "https://files.rcsb.org/download/{}.pdb"

# BRIL / cytochrome b562 core (variant-resistant search seed, near the C-term)
BRIL_CORE = "RHGFDILVGQIDDALKLANEGKVKEAQAAAEQLKTTRNAYIQKYLERARSTLDKLNDDIFK"
BRIL_HEAD_DEFAULT = "YDLEDNWETLNDNLKVIEKADNAAQVKDALTKMRAAALDAQKATMKDFRHGFDILVGQIDDALKLANEGKVKEAQAAAEQLKTTRNAYIQKYLERARSTLDKLNDDIFK"

# --- residue map ------------------------------------------------------------
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
       "MSE":"M","SEC":"C","HSD":"H","HSE":"H","HSP":"H","HID":"H","HIE":"H"}

# ---------------------------------------------------------------------------
def fetch(pid, outdir):
    path = os.path.join(outdir, pid + ".pdb")
    if not os.path.exists(path):
        print(f"downloading {pid} -> {path}")
        req = urllib.request.Request(RCSB.format(pid), headers={"User-Agent": "at1r-nb-suite/1.0"})
        with urllib.request.urlopen(req, context=ssl.create_default_context(), timeout=60) as r:
            data = r.read()
        if data.startswith(b"<!DOCTYPE") or b"NOT FOUND" in data[:400]:
            raise RuntimeError(f"{pid}: RCSB did not return a coordinate file")
        with open(path, "wb") as fh:
            fh.write(data)
    else:
        print(f"using cached {path}")
    return path

def read_chain(pid, outdir, chain):
    """Parse ATOM records (first alt loc, no HETATM) -> list of residues."""
    path = fetch(pid, outdir)
    residues = []          # [resseq, ins, resname, [atom_lines]]
    seen = set()
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            if line[16] != " " or rec == "HETATM":
                continue                       # skip alt locs / hetero
            if line[21] != chain:
                continue
            key = (line[22:26].strip(), line[26].strip())
            if key not in seen:
                seen.add(key)
                residues.append([key[0], key[1], line[17:20].strip(), []])
            residues[-1][3].append(line)
    if not residues:
        raise RuntimeError(f"{pid} chain {chain}: no ATOM records found")
    return residues

def seq_of(residues):
    """One-letter sequence for a chain, 'X' for unknown residues."""
    out = []
    for _rseq, _ins, resname, _a in residues:
        out.append(AA3.get(resname.upper(), "X"))
    return "".join(out)

def find_ix(seq, sub):
    i = seq.find(sub)
    return i if i >= 0 else seq.upper().find(sub.upper())

def _atom_xyz(line):
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def _heavy_lines(residues):
    out = []
    for _rseq, _ins, _resname, atoms in residues:
        for a in atoms:
            elem = (a[76:78] or a[12:16]).strip().upper()
            if elem.startswith("H") or elem == "D":
                continue
            out.append(a)
    return out


def _rigid_transform(residues, R, t):
    """Apply xyz' = xyz @ R.T + t to every atom line of the residue lists."""
    out = []
    for rseq, ins, resname, atoms in residues:
        new_atoms = []
        for a in atoms:
            x, y, z = _atom_xyz(a) @ R.T + t
            new_atoms.append(f"{a[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{a[54:]}")
        out.append([rseq, ins, resname, new_atoms])
    return out


_BIG = 1e9


def rigid_fit_to_interface(nb, rec, gap=2.4):
    """Rigid-body docking of the nanobody onto the receptor ECL pocket.

    The deposited 8TH3 model (locally-refined cryo-EM) leaves the nanobody's
    CDR tip ~2.4 nm off the receptor; the paratope already points at the
    epitope.  We sample small rotations (0-12 deg) and translations (insertion
    along the pocket surface normal + lateral jitter), then hill-climb the
    best seeds by translation only.  The score rewards heavy-atom contacts
    (0.30-0.50 nm) and near-contacts (0.50-0.65 nm) and penalises clashes
    (< 0.30 nm).  A drift cap keeps the search near the epitope.

    Returns (nb_residues shifted, delta[nm], n_contacts, n_clashes).
    """
    NB = np.array([_atom_xyz(a) for a in _heavy_lines(nb)])
    RC = np.array([_atom_xyz(a) for a in _heavy_lines(rec)])
    cent = NB.mean(0)

    def atom_min(X, Y):
        mn = np.full(X.shape[0], _BIG)
        for s in range(0, Y.shape[0], 400):
            B = Y[s:s + 400]
            D = X[:, None, :] - B[None, :, :]
            mn = np.minimum(mn, np.sqrt((D * D).sum(-1)).min(1))
        return mn

    # deposited tip atoms (< 3.5 nm from receptor) define the CDR contact edge
    tip_mask = atom_min(NB, RC) < 3.5
    tip_c = NB[tip_mask].mean(0) if tip_mask.any() else NB.mean(0)

    # epitope patch: receptor atoms within 4.0 nm of any tip atom
    od = NB[tip_mask][:, None, :] - RC[None, :, :]
    near = np.sqrt((od * od).sum(-1)).min(0) < 4.0
    patch = RC[near]
    epi_c = patch.mean(0)

    # surface normal via PCA; oriented back toward the nanobody
    pc = patch - epi_c
    u, _s, _vt = np.linalg.svd(pc.T, full_matrices=False)
    normal = u[:, -1]
    if np.dot(tip_c - epi_c, normal) < 0:
        normal = -normal

    def score(X):
        if not np.all(np.isfinite(X)):
            return -np.inf
        mn = atom_min(X, RC)
        clash = int((mn < 0.22).sum()) * 3.0 \
            + int(((mn >= 0.22) & (mn < 0.26)).sum()) * 2.0 \
            + int(((mn >= 0.26) & (mn < 0.30)).sum())
        contact = int(((mn >= 0.30) & (mn < 0.50)).sum())
        near_c = int(((mn >= 0.50) & (mn < 0.65)).sum()) * 0.5
        return contact + near_c - clash

    def rot_matrix(axis, ang_deg):
        a = axis / np.linalg.norm(axis)
        th = np.radians(ang_deg)
        K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
        return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)

    def translation_from(R, d):
        """Pure translation so that (X-cent)@R.T+cent+d == X@R.T+t."""
        return cent - cent @ R.T + d

    rng = np.random.default_rng(7)
    axes = [np.array(a) for a in ([0, 0, 1], [1, 0, 0], [0, 1, 0])]
    axes += [rng.normal(size=3) for _ in range(6)]

    seeds = []                             # (score, R, t)
    for offs in np.arange(0.30, 0.95, 0.10):     # insertion depth
        for _ in range(10):                      # lateral jitter
            base = (epi_c + normal * offs) - tip_c + rng.normal(0, 0.25, 3)
            for ang in (0.0, 6.0, 12.0):         # small rotations
                for ax in axes:
                    R = rot_matrix(ax, ang)
                    t = translation_from(R, base)
                    with np.errstate(all="ignore"):
                        X = NB @ R.T + t
                    if np.all(np.isfinite(X)):
                        seeds.append((score(X), R, t))

    seeds.sort(key=lambda s: -s[0])
    best_s, best_R, best_t = seeds[0]
    for s0, R, t0 in seeds[:6]:
        sX, t, start_t, sigma = s0, t0, t0, 0.05
        for _ in range(160):
            can_t = t + rng.normal(0, sigma, 3)
            if np.linalg.norm(can_t - start_t) > 1.2:    # drift cap
                sigma = min(0.3, sigma * 1.2)
                continue
            with np.errstate(all="ignore"):
                cs = score(NB @ R.T + can_t)
            if cs > sX:
                t, sX, sigma = can_t, cs, 0.05
            else:
                sigma = min(0.3, sigma * 1.2)
        if sX > best_s:
            best_s, best_R, best_t = sX, R, t

    with np.errstate(all="ignore"):
        X = NB @ best_R.T + best_t
    mn = atom_min(X, RC)
    n_contacts = int(((mn >= 0.30) & (mn < 0.50)).sum())
    n_clash = int((mn < 0.26).sum())
    # the pure translation component (detached from the rotation)
    pure_t = best_t - (cent - cent @ best_R.T)
    print(f"[fit] contact atoms: {n_contacts}, clashes: {n_clash}, "
          f"displacement: {np.linalg.norm(pure_t):.2f} nm, "
          f"score: {best_s:.1f}")
    if n_contacts < 6:
        print("[fit] WARNING: weak interface after rigid refit; "
              "verify the pose before simulating")
    shifted = _rigid_transform(nb, best_R, best_t)
    return shifted, pure_t, n_contacts, n_clash


def write_chunks(pid, tag, chunks, outdir, chains=("A", "B")):
    """Write residue chunks as consecutive chains, renumbered 1..N."""
    out = os.path.join(outdir, f"{pid}_{tag}.pdb")
    lines, serial = [], 1
    for i, piece in enumerate(chunks):
        ch = chains[i] if i < len(chains) else "C"
        prev_key, rnum = None, 1
        for rseq, ins, resname, atoms in piece:
            rnum = 1 if prev_key is None else rnum + (rseq != prev_key or ins != prev_ins)
            prev_key, prev_ins = rseq, ins
            for a in atoms:
                lines.append(f"{a[:6]}{serial:5d}{a[11:21]}{ch}{rnum:4d}{a[27:]}")
                serial += 1
        lines.append("TER")
    lines.append("END")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    n = sum(len(c) for c in chunks)
    print(f"[{pid}] wrote {tag}: {n} aa in {len(chunks)} chain(s) -> {out}")
    return out

# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb-dir", default=os.path.join(os.path.dirname(__file__), "data"))
    ap.add_argument("-c", "--chain", default="A", help="nanobody copy to keep (7T83/7T84)")
    ap.add_argument("--strip-bril", action="store_true",
                    help="carve the BRIL insert out of the 8TH3 receptor chain A")
    ap.add_argument("--no-rebuild", action="store_true",
                    help="keep the deposited (gapped) nanobody pose when writing "
                         "8TH3_complex.pdb; default is a rigid refit")
    args = ap.parse_args()
    outdir = args.pdb_dir
    os.makedirs(outdir, exist_ok=True)

    # ----------------------------------------------------------
    # 1. Free nanobody: 7T83 (WT) and 7T84 (G26D T57I)
    # ----------------------------------------------------------
    for pid in ("7T83", "7T84"):
        residues = read_chain(pid, outdir, args.chain)
        seq = seq_of(residues)
        print(f"[{pid}] chain {args.chain}: {len(residues)} aa "
              f"({seq.count('X')} unknown res)")
        write_chunks(pid, "nb", [residues], outdir, chains=(args.chain,))

    # ----------------------------------------------------------
    # 2. Complex: 8TH3 (nanobody chain B + receptor chain A)
    # ----------------------------------------------------------
    pid = "8TH3"
    nb_chain = "B"
    rec_chain = "A"
    nb = read_chain(pid, outdir, nb_chain)
    rec = read_chain(pid, outdir, rec_chain)
    print(f"[{pid}] nanobody chain {nb_chain}: {len(nb)} aa "
          f"({seq_of(nb).count('X')} unknown)")
    print(f"[{pid}] receptor chain {rec_chain}: {len(rec)} aa "
          f"({seq_of(rec).count('X')} unknown)")

    rec_out = [rec]
    if args.strip_bril:
        seq = seq_of(rec)
        head = find_ix(seq, BRIL_HEAD_DEFAULT)
        if head < 0:
            head = find_ix(seq, BRIL_CORE)
            if head >= 0:
                head -= len(BRIL_HEAD_DEFAULT) - seq.find(BRIL_CORE)
        if head < 0:
            print("[8TH3] WARNING: BRIL not found; receptor kept intact")
        else:
            start = head
            end = head + len(BRIL_HEAD_DEFAULT)
            if find_ix(seq, BRIL_CORE) >= 0:
                end = find_ix(seq, BRIL_CORE) + len(BRIL_CORE)
                start = max(0, end - len(BRIL_HEAD_DEFAULT))
            print(f"[8TH3] stripping BRIL construct-ix {start}:{end}")
            rec_out = [rec[:start], rec[end:]]

    write_chunks(pid, "nb", [nb], outdir, chains=(nb_chain,))
    write_chunks(pid, "at1r", rec_out, outdir, chains=(rec_chain,))

    # ---- nanobody/receptor complex with rigid refit ------------------------
    meta = {"pdb": pid, "rebuilt": False, "note": "",
            "min_distance_nm": None, "displacement_nm": None,
            "n_contacts": 0, "n_clashes": 0}
    nb_fit, complex_chunks = nb, [nb] + rec_out
    if not args.no_rebuild:
        na = np.array([_atom_xyz(a) for a in _heavy_lines(nb)])
        ra = np.array([_atom_xyz(a) for a in _heavy_lines(rec)])
        gap = float(np.sqrt(((na[:, None, :] - ra[None, :, :]) ** 2).sum(-1)).min())
        meta["min_distance_nm"] = round(gap, 3)
        if gap < 1.0:
            print(f"[8TH3] nanobody already in contact (gap {gap:.2f} nm); no refit")
        else:
            print(f"[8TH3] deposited gap {gap:.2f} nm -> rigid refit onto receptor")
            nb_fit, delta, n_contact, n_clash = rigid_fit_to_interface(nb, rec, gap)
            meta.update(rebuilt=True,
                    displacement_nm=round(float(np.linalg.norm(delta)), 4),
                    n_contacts=int(n_contact), n_clashes=int(n_clash),
                    note="locally-refined cryo-EM deposit leaves the nanobody "
                         "~2.4 nm off AT1R; rigidly refit (sample small rotations "
                         "and translations) to re-establish a physical interface")
            complex_chunks = [nb_fit] + rec_out
    write_chunks(pid, "complex", complex_chunks, outdir, chains=("B", "A"))
    with open(os.path.join(outdir, "8TH3_complex.meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")

if __name__ == "__main__":
    main()