"""Streamlit companion to the AT118i4h32 / PDB 7T83 project.

Highlights the same analyses as the static site and the Colab notebook:
crystallographic checks, sequence anatomy, a box-size / run planner, example
MD observables (upload your own .xvg), and the 8TH3 nanobody-receptor
interface contact map.

Run:  streamlit run app.py
"""

import os
import math
import pathlib

import numpy as np
import pandas as pd
import streamlit as st

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = pathlib.Path(__file__).resolve().parent

# --------------------------------------------------------------------------
# Constants (deposited values for 7T83, also returned by the RCSB API)
# --------------------------------------------------------------------------
PDB_ID = "7T83"
CELL = {"a": 80.11, "b": 85.40, "c": 84.44, "beta_deg": 108.64}
Z_SU = 2  # P 1 21 1 -> 2 symmetry ops per unit cell
N_CHAINS_ASU = 8
CITRATE_MASS = 189.11  # FLC, C6H8O7

SEQ = ("MAEVQLVESGGGLVQPGGSLRLSCAASGYIYRRYRMGWYRQAPGKGREFVAAIS"
       "GGSSTNYADSVKGRFTISRDNSKNTVYLQMNSLRAEDTAVYYCAAYRIVSDPRVY"
       "WGQGTQVTVSSLEHHHHHH")

RES_MASS = {"G": 57.0519, "A": 71.0788, "S": 87.0782, "P": 97.1167,
            "V": 99.1326, "T": 101.1051, "C": 103.1388, "L": 113.1594,
            "I": 113.1594, "N": 114.1038, "D": 115.0886, "Q": 128.1307,
            "K": 128.1741, "E": 129.1155, "M": 131.1926, "H": 137.1411,
            "F": 147.1766, "R": 156.1875, "Y": 163.176, "W": 186.2132}

KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}

REGIONS = [("FR1", 0, 25), ("CDR1", 25, 33), ("FR2", 33, 52),
           ("CDR2", 52, 60), ("FR3", 60, 98), ("CDR3", 98, 110),
           ("FR4", 110, 121), ("His-tag", 121, 128)]

CDR_MOTIFS = {"CDR1": "GYIYRRYR", "CDR2": "ISGGSSTN", "CDR3": "RIVSDPRVY"}
RESIDUE_CHARGE = {aa: 0 for aa in set(KD)}
RESIDUE_CHARGE.update({a: 1 for a in "KR"})
RESIDUE_CHARGE.update({a: -1 for a in "DE"})


def plot(fig):
    st.pyplot(fig)
    plt.close(fig)


# --------------------------------------------------------------------------
# Small pure-python helpers
# --------------------------------------------------------------------------
def chain_mass(seq):
    return sum(RES_MASS[aa] for aa in seq) + 18.015


def cell_volume(cell):
    return cell["a"] * cell["b"] * cell["c"] * math.sin(
        math.radians(cell["beta_deg"]))


def matthews(n_chains, v, mass):
    vm = v / (n_chains * mass * Z_SU)
    solvent = 1 - 1.23 / vm
    return vm, solvent


def load_xvg_file(uploaded):
    """Parse a GROMACS .xvg upload (first two columns)."""
    data = []
    for line in uploaded.getvalue().decode("utf-8", "replace").splitlines():
        s = line.strip()
        if not s or s.startswith(("#", "@", "&")):
            continue
        tok = s.split()
        try:
            data.append([float(tok[0]), float(tok[1])])
        except (IndexError, ValueError):
            continue
    return (np.asarray(data)[:, 0], np.asarray(data)[:, 1]) if data else None


def read_pdb_heavy(path):
    """[(resnum, resname, (N,3) heavy-atom coords)] for one PDB file."""
    out, cur = [], None
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM") or line[16] != " ":
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


def min_pair_distance(a, b):
    d = a[:, None, :] - b[None, :, :]
    return float(np.sqrt((d * d).sum(-1)).min())


# --------------------------------------------------------------------------
# Pages
# --------------------------------------------------------------------------
def page_overview():
    st.markdown("""
**Nanobody AT118i4h32** is an antagonist nanobody of the angiotensin II
type 1 receptor (AT1R), deposited as **PDB 7T83** at 2.10 A resolution
(X-ray, P 1 21 1). The crystal holds 8 copies of the same 128-residue
construct (chains A-H, monomer assembly).
""")
    cols = st.columns(4)
    cols[0].metric("Resol.", "2.10 A")
    cols[1].metric("Chains", "8 x 128 aa")
    cols[2].metric("Ligand", "Citrate x4")
    cols[3].metric("R-work / R-free", "0.208 / 0.257")

    v = cell_volume(CELL)
    m = chain_mass(SEQ)
    cols = st.columns(3)
    cols[0].metric("Cell volume", f"{v:,.0f} A^3")
    cols[1].metric("Chain mass", f"{m/1000:.2f} kDa")
    cols[2].metric("Dep. ASU mass", f"{8*m + 4*CITRATE_MASS:.2f} kDa")

    st.markdown("""
Reference facts: a=80.11 A, b=85.40 A, c=84.44 A, beta=108.64 deg.

Sections:
- **Crystallography** - Matthews coefficient / solvent-content sweep for 1-12 chains per ASU.
- **Sequence** - 128-res construct, CDR loops, charge and Kyte-Doolittle hydropathy.
- **Box planner** - water/atom budget for a GROMACS run box at a chosen padding.
- **MD observables** - upload `rmsd.xvg` / `gyrate.xvg` / `hbond.xvg`, or view an illustrative demo.
- **8TH3 interface** - static nanobody-receptor contact map and table.
""")


def page_crystallography():
    st.subheader("Cell volume & Matthews coefficient")
    v = cell_volume(CELL)
    m = chain_mass(SEQ)
    st.markdown(
        f"Monoclinic cell: *V = a b c sin(beta)* = "
        f"{CELL['a']} x {CELL['b']} x {CELL['c']} x sin({CELL['beta_deg']}) = "
        f"**{v:,.0f} A^3**.")

    rows = []
    for n in [1, 2, 4, 6, 8, 12]:
        vm, solv = matthews(n, v, m)
        rows.append({"chains/ASU": n, "Vm (A^3/Da)": round(vm, 2),
                     "solvent": f"{solv:.1%}",
                     "verdict": ("in usual 1.7-3.5 range" if 1.7 <= vm <= 3.5
                                 else "outside")})
    st.dataframe(pd.DataFrame(rows), hide_index=True)

    ns = np.array([1, 2, 4, 6, 8, 12])
    vms = np.array([matthews(n, v, m)[0] for n in ns])
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axhspan(1.7, 3.5, alpha=0.2, color="k")
    ax.plot(ns, vms, "o-")
    for x, y in zip(ns, vms):
        ax.text(x, y + 0.3, f"{y:.2f}", ha="center")
    ax.set_xlabel("chains per asymmetric unit")
    ax.set_ylabel("Matthews Vm (A^3/Da)")
    ax.set_title("Only ~8 chains/ASU sits in the usual range")
    ax.grid(alpha=0.3)
    plot(fig)
    vm8, solv8 = matthews(8, v, m)
    st.markdown(f"**Verdict:** 8 chains/ASU gives Vm = {vm8:.2f} A^3/Da and "
                f"about {solv8:.0%} solvent - the deposited packing.")


def page_sequence():
    st.subheader("Sequence anatomy (128 residues)")
    st.markdown("Construct: " + SEQ)

    df = []
    for name, s, e in REGIONS:
        frag = SEQ[s:e]
        q = sum(RESIDUE_CHARGE[x] for x in frag)
        h = sum(KD[x] for x in frag) / len(frag)
        df.append({"region": name, "range": f"{s+1}-{e}", "len": len(frag),
                   "seq": frag, "charge": q, "mean KD": round(h, 2)})
    st.dataframe(pd.DataFrame(df), hide_index=True)

    q_full = sum(RESIDUE_CHARGE[x] for x in SEQ)
    st.markdown(
        f"Full charge (K+R minus D+E, His neutral): **{q_full:+d}**. "
        f"CDR1 `{CDR_MOTIFS['CDR1']}`, CDR2 `{CDR_MOTIFS['CDR2']}`, "
        f"CDR3 `{CDR_MOTIFS['CDR3']}`. Cysteines at "
        f"{[i+1 for i, x in enumerate(SEQ) if x == 'C']} (conserved C24-C97 "
        f"disulfide); histidines only in the C-terminal His-tag.")

    win = st.slider("Hydropathy window (+/-)", 1, 10, 5)
    prof = np.array([np.mean([KD[SEQ[j]] for j in
                              range(max(0, i - win),
                                    min(len(SEQ), i + win + 1))])
                     for i in range(len(SEQ))])
    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(np.arange(1, len(SEQ) + 1), prof)
    for name, s, e in REGIONS:
        if "CDR" in name:
            ax.axvspan(s + 1, e, alpha=0.25, color="#2980b9")
    ax.set_xlabel("residue number")
    ax.set_ylabel(f"mean Kyte-Doolittle (window +/-{win})")
    ax.set_title("Hydropathy profile, shaded = CDR loops")
    ax.grid(alpha=0.3)
    plot(fig)


def page_box():
    st.subheader("GROMACS box-size estimator")
    pad = st.slider("Padding to box edge (nm)", 0.5, 3.0, 1.0, 0.1)
    prot_d = st.slider("Protein diameter (nm)", 2.0, 8.0, 4.0, 0.1)
    prot_atoms = st.slider("Protein atoms", 500, 6000, 1900, 100)

    edge = prot_d + 2 * pad
    box_v = edge ** 3
    waters = int((box_v - 17.0) * 33.4)
    atoms = prot_atoms + waters * 3
    ion_pairs = int(waters * 0.15 / 55.5)

    c = st.columns(4)
    c[0].metric("Edge", f"{edge:.1f} nm")
    c[1].metric("Box volume", f"{box_v:.0f} nm^3")
    c[2].metric("Water", f"{waters:,}")
    c[3].metric("Atoms", f"{atoms:,}")
    st.caption(f"0.15 M NaCl estimate: ~{ion_pairs:,} Na+ / Cl- ion pairs.")

    base = 1.0
    st.markdown(f"Cost scales as volume: a {pad:.1f} nm pad is roughly "
                f"**{(edge / (prot_d + 2*base)) ** 3:.1f}x** the atoms of a "
                f"{base:g} nm pad.")


def page_md():
    st.subheader("MD observables (.xvg)")
    use_demo = st.toggle("Use illustrative synthetic data", value=True)

    labels = {"rmsd": ("RMSD (nm)", "#2e86c1"),
              "gyrate": ("Rg (nm)", "#1e8449"),
              "hbond": ("H-bonds", "#b03a2e")}
    series = {}
    if use_demo:
        rng = np.random.default_rng(7)
        t = np.arange(0, 100, 0.1)
        series["rmsd"] = (t, 0.18 + 0.06 * (1 - np.exp(-t / 12)) +
                          0.015 * np.sin(t * 0.7) + rng.normal(0, 0.008, len(t)))
        series["gyrate"] = (t, 1.52 + 0.015 * (1 - np.exp(-t / 20)) +
                            rng.normal(0, 0.004, len(t)))
        series["hbond"] = (t, 95 + 6 * np.sin(t * 0.15) +
                           rng.normal(0, 2.5, len(t)))
        note = "ILLUSTRATIVE synthetic observables - upload real .xvg files to replace them."
    else:
        note = "Upload GROMACS .xvg files (two columns: time, value)."
        for key in labels:
            up = st.file_uploader(labels[key][0], type=["xvg", "txt"],
                                  key=f"xvg_{key}")
            if up is not None:
                parsed = load_xvg_file(up)
                if parsed is not None:
                    series[key] = parsed

    if note:
        st.caption(note)

    if not series:
        st.info("No data to plot yet - select the demo or upload files.")
        return

    fig, axes = plt.subplots(len(series), 1, figsize=(9, 2.4 * len(series)),
                             sharex=True)
    if len(series) == 1:
        axes = [axes]
    for ax, key in zip(axes, series):
        t, y = series[key]
        ax.plot(t, y, color=labels[key][1], lw=0.8)
        ax.set_ylabel(labels[key][0])
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("time")
    fig.tight_layout()
    plot(fig)

    res = np.arange(1, len(SEQ) + 1)
    rmsf = (0.08 + 0.05 * np.exp(-((res - 104) / 6) ** 2) +
            0.03 * np.exp(-((res - 29) / 4) ** 2) + 0.04 * (res > 121) +
            np.random.default_rng(7).normal(0, 0.008, len(res)))
    rmsf = np.clip(rmsf, 0.03, None)
    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.plot(res, rmsf)
    ax.set_xlabel("residue")
    ax.set_ylabel("RMSF (nm)")
    ax.set_title("ILLUSTRATIVE per-residue RMSF - expect peaks in CDR loops / termini")
    ax.grid(alpha=0.3)
    plot(fig)


def page_interface():
    st.subheader("8TH3 nanobody-on-receptor interface")
    csv_path = BASE / "simulation" / "results" / "interface.csv"
    pdb_path = BASE / "simulation" / "data" / "8TH3_complex.pdb"
    if not (csv_path.is_file() and pdb_path.is_file()):
        st.info("interface.csv / 8TH3_complex.pdb not committed yet.")
        return

    meta = st.expander("About this model")
    meta.markdown("""
The cryo-EM deposit 8TH3 leaves the nanobody ~2.4 nm off AT1R; it was
rigidly translated (no rotation) to re-establish a physical interface
before analysis. Contact table uses a 0.45 nm heavy-atom cutoff.
""")

    df = pd.read_csv(csv_path)
    nb = df[df["chain"] == "B"]
    rec = df[df["chain"] == "A"]
    cols = st.columns(3)
    cols[0].metric("Interface res (nb)", int((nb.in_interface == 1).sum()))
    cols[1].metric("Interface res (receptor)", int((rec.in_interface == 1).sum()))
    cols[2].metric("Total contacts", int(df.n_contacts.sum()))

    with st.expander("Interface-ready residue table"):
        st.dataframe(nb[nb.in_interface == 1][["chain", "resnum", "resname",
                                               "n_contacts"]],
                     hide_index=True)
        st.dataframe(rec[rec.in_interface == 1][["chain", "resnum", "resname",
                                                 "n_contacts"]],
                     hide_index=True)

    nb_atoms = read_pdb_heavy(str(pdb_path))
    chain_b = [(r, n, a) for r, n, a in nb_atoms if n == "B"]
    chain_a = [(r, n, a) for r, n, a in nb_atoms if n == "A"]
    if chain_b and chain_a:
        contacts = np.zeros((len(chain_b), len(chain_a)), dtype=float)
        for i, (r1, n1, a) in enumerate(chain_b):
            for j, (r2, n2, b) in enumerate(chain_a):
                if min_pair_distance(a, b) <= 0.45:
                    contacts[i, j] = 1.0
        nbc = contacts.sum(1)
        idx = np.nonzero(nbc)[0]
        if idx.size:
            fig, ax = plt.subplots(figsize=(9, 4))
            ax.imshow(contacts[idx][:, np.nonzero(contacts.sum(0))[0]].T,
                      aspect="auto", cmap="cividis", interpolation="nearest")
            ax.set_yticks(range(len(np.nonzero(contacts.sum(0))[0])))
            ax.set_yticklabels([chain_a[j][0] for j in
                                np.nonzero(contacts.sum(0))[0]])
            ax.set_xticks(range(len(idx)))
            ax.set_xticklabels([chain_b[i][0] for i in idx], rotation=90)
            ax.set_xlabel("nanobody residue (chain B)")
            ax.set_ylabel("receptor residue (chain A)")
            ax.set_title(f"Static contact map, {int((contacts > 0).sum())} "
                         f"residue pairs < 0.45 nm")
            plot(fig)
        else:
            st.info("No contacts within 0.45 nm in the reference frame.")


# --------------------------------------------------------------------------
# App shell
# --------------------------------------------------------------------------
def main():
    st.set_page_config(page_title="Nanobody AT118i4h32 (7T83)",
                       page_icon=None,
                       layout="wide")
    st.title("Nanobody AT118i4h32 (PDB 7T83) - Structure & Dynamics")
    st.caption("AT1R antagonist nanobody | X-ray 2.10 A | 8 x 128-residue "
               "chains | GROMACS MD")

    tabs = st.tabs(["Overview", "Crystallography", "Sequence", "Box planner",
                    "MD observables", "8TH3 interface"])
    pages = [page_overview, page_crystallography, page_sequence, page_box,
             page_md, page_interface]
    for tab, page in zip(tabs, pages):
        with tab:
            page()

    st.divider()
    st.markdown("Web: https://Harshini-Suresha.github.io/at1r-nanobody-7T83/ "
                "| Repo: github.com/Harshini-Suresha/at1r-nanobody-7T83 | "
                "RCSB: rcsb.org/structure/7T83")


if __name__ == "__main__":
    main()