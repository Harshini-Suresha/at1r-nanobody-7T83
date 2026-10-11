"""
xvg_utils.py -- helpers shared by 03_analyze.py and the Colab notebook.

* load_xvg            : parse GROMACS .xvg (rmsd / gyrate / hbond / rmsf)
* needleman_wunsch    : residue-residue map between variants (7T83 vs 7T84)
* aa_from_pdb         : residue number -> one-letter code from a _nb.pdb
* runstats            : replicate mean / std / SE / t-stat helpers
"""
import os
import numpy as np

BLOSUM_RES = None  # simple +1/-3 scoring is sufficient for 100-aa constructs


def load_xvg(path, tcol=0, vcol=1):
    """Columns of a GROMACS .xvg file as two numpy arrays (skips # and @)."""
    t, v = [], []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(("#", "@", "&")):
                continue
            tok = s.split()
            try:
                t.append(float(tok[tcol])); v.append(float(tok[vcol]))
            except (IndexError, ValueError):
                continue
    return np.asarray(t), np.asarray(v)


def load_rmsf(path):
    """Per-residue RMSF (.xvg from `gmx rmsf -res`): (resids, rmsf)."""
    return load_xvg(path, 0, 1)


def aa_from_pdb(pdb_path):
    """Return {residue_number: one_letter_aa} for the first chain of a PDB."""
    aa3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
           "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
           "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
           "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
           "MSE": "M"}
    out, prev = {}, None
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            rn = int(line[22:26].strip())
            if rn != prev:
                prev = rn
                out[rn] = aa3.get(line[17:20].strip(), "X")
    return out


def needleman_wunsch(a, b, match=2, mismatch=-1, gap=-2):
    """Global alignment -> list of (i_a, i_b) aligned indices.

    i_a / i_b index operator[] on the *aa lists*; None = gap.
    """
    n, m = len(a), len(b)
    S = np.zeros((n + 1, m + 1), dtype=float)
    T = np.zeros((n + 1, m + 1), dtype=np.int8)
    for i in range(1, n + 1):
        S[i, 0], T[i, 0] = i * gap, 1           # up
    for j in range(1, m + 1):
        S[0, j], T[0, j] = j * gap, 2           # left
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = S[i - 1, j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up = S[i - 1, j] + gap
            left = S[i, j - 1] + gap
            best = max(diag, up, left)
            S[i, j] = best
            T[i, j] = 0 if best == diag else (1 if best == up else 2)
    # backtrace: pairs of (i_a, i_b) ORIGIN indices (None = gap)
    i, j = n, m
    ia, ib = [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and T[i, j] == 0:      # diagonal: both consumed
            ia.append(i - 1); ib.append(j - 1); i -= 1; j -= 1
        elif i > 0 and T[i, j] == 1:              # "up": a[i-1] hits a gap in b
            ia.append(i - 1); ib.append(None); i -= 1
        else:                                     # "left": b[j-1] hits a gap in a
            ia.append(None); ib.append(j - 1); j -= 1
    ia.reverse(); ib.reverse()
    return list(zip(ia, ib))


def runstats(values):
    """values: (nrep, ...) -> (mean, sd_across_reps, se)."""
    arr = np.asarray(values, dtype=float)
    mean = arr.mean(axis=0)
    sd = arr.std(axis=0, ddof=1) if arr.shape[0] > 1 else np.zeros_like(mean)
    se = sd / np.sqrt(arr.shape[0])
    return mean, sd, se


def welch_t(delta, se_a, se_b):
    return delta / np.sqrt(se_a ** 2 + se_b ** 2 + 1e-12)


def natural_sort_key(k):
    import re
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", k)]


def find_reps(root, system):
    """Rep dirs <root>/<system>/rep*_seed* containing rmsd.xvg."""
    base = os.path.join(root, system)
    if not os.path.isdir(base):
        return []
    reps = []
    for d in sorted(os.listdir(base)):
        if d.startswith("rep") and os.path.isfile(os.path.join(base, d, "rmsd.xvg")):
            reps.append(os.path.join(base, d))
    reps.sort(key=natural_sort_key)
    return reps