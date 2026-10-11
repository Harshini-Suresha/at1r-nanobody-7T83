# Nanobody AT118i4h32 (PDB 7T83) — Structure & Dynamics

An interactive structural and molecular-dynamics study of **AT118i4h32**, an antagonist nanobody of the angiotensin II type 1 receptor (AT1R), built from **PDB entry [7T83](https://www.rcsb.org/structure/7T83)** at 2.10 Å.

> One nanobody. Eight copies. A simulation — crystal, sequence, binding loops, and stability in water.

## 🌐 Website (GitHub Pages)

**Live site:** https://Harshini-Suresha.github.io/at1r-nanobody-7T83/

The site is a single self-contained `index.html` (React + in-browser calculations: cell volume, Matthews coefficient, sequence charge/hydropathy, box-size estimator, run planner).

## 🚀 Streamlit app

An interactive companion app reproduces the site's in-browser calculations
plus a 7T83-vs-7T84 variant comparison and the 8TH3 nanobody-receptor contact
map, with GROMACS `.xvg` upload for your own trajectories.

```bash
pip install -r requirements.txt
streamlit run app.py
```

Or deploy on [Streamlit Community Cloud](https://share.streamlit.io) with
this repo, main file `app.py`.

## 🚀 Open in Colab

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Harshini-Suresha/at1r-nanobody-7T83/blob/main/notebooks/AT118i4h32_analysis.ipynb)

- **Colab link:** https://colab.research.google.com/github/Harshini-Suresha/at1r-nanobody-7T83/blob/main/notebooks/AT118i4h32_analysis.ipynb
- **Raw ipynb:** https://github.com/Harshini-Suresha/at1r-nanobody-7T83/blob/main/notebooks/AT118i4h32_analysis.ipynb
- **Raw (download):** https://raw.githubusercontent.com/Harshini-Suresha/at1r-nanobody-7T83/main/notebooks/AT118i4h32_analysis.ipynb
- **nbviewer:** https://nbviewer.org/github/Harshini-Suresha/at1r-nanobody-7T83/blob/main/notebooks/AT118i4h32_analysis.ipynb

## 📁 Contents

```
at1r-nanobody-7T83/
├── index.html                        # project website (deployed to Pages)
├── app.py                            # Streamlit companion app
├── requirements.txt
├── notebooks/
│   ├── AT118i4h32_analysis.ipynb     # Colab-ready analysis (open badge above)
│   └── MMP7_docking.ipynb            # MMP-7 docking lab
├── simulation/                       # GROMACS MD suite (scripts + data/results)
├── README.md
└── .nojekyll
```

## 🔬 What the notebook does

1. Fetches 7T83 metadata live from RCSB (method, resolution, chains, ligands, cell).
2. Recomputes crystallographic checks: cell volume, chain mass, Matthews coefficient / solvent content.
3. Analyses the 128-residue sequence (`MAEVQLVES...LEHHHHHH`): CDR/FR regions, charge, Kyte–Doolittle hydropathy, disulfide, His-tag.
4. Plots hydropathy, box-size estimator, and example MD observables (RMSD / RMSF / Rg / H-bonds) with hooks to load your own GROMACS `.xvg` files.
5. Links back to the website, RCSB, and PDB download.

Run it locally:

```bash
pip install biopython matplotlib numpy requests
jupyter notebook notebooks/AT118i4h32_analysis.ipynb
```

## 🧪 Key facts (from RCSB 7T83)

| Item | Value |
|------|-------|
| PDB | 7T83 — X-ray, 2.10 Å |
| Construct | 128 residues × 8 chains (A–H), monomer assembly |
| Ligand | Citrate (FLC) × 4 |
| Cell | a=80.11 Å, b=85.40 Å, c=84.44 Å, β=108.64°, P 1 21 1 |
| R-work / R-free | 0.208 / 0.257 |
| Total weight | ~114.56 kDa |
| Tools | BIOVIA Discovery Studio, PyMOL, GROMACS |

Related receptor-containing entries are listed on the website under *Related structures*.

## 🔗 Links

- Website: https://Harshini-Suresha.github.io/at1r-nanobody-7T83/
- Repo: https://github.com/Harshini-Suresha/at1r-nanobody-7T83
- Colab: https://colab.research.google.com/github/Harshini-Suresha/at1r-nanobody-7T83/blob/main/notebooks/AT118i4h32_analysis.ipynb
- RCSB 7T83: https://www.rcsb.org/structure/7T83
- PDB file: https://files.rcsb.org/download/7T83.pdb
- GROMACS: https://manual.gromacs.org/

## 📄 Citation

If you use this work, cite PDB 7T83 and the associated polyreactivity study (see website *Methodology* / *How the polyreactivity work relates*), plus GROMACS.
