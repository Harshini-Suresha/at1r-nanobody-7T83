#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Shared configuration for the AT118i4h32 (7T83 / 7T84) + 8TH3 simulation suite.
# Source this from every other simulation script:  source config.sh
# ---------------------------------------------------------------------------
set -euo pipefail

# --- workspace -------------------------------------------------------------
export SIM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATA="${SIM_ROOT}/data"          # downloaded PDBs + split chains
export WORK="${SIM_ROOT}/work"          # per-system GROMACS working dirs
export RESULTS="${SIM_ROOT}/results"    # xvg / csv / figures
mkdir -p "${DATA}" "${WORK}" "${RESULTS}"

# --- structures ------------------------------------------------------------
export PDB_WT="7T83"          # AT118i4h32, X-ray 2.10 A, chains A-H (nanobody only)
export PDB_VAR="7T84"         # AT118i4h32 G26D T57I, X-ray 1.60 A, chains A-H
export PDB_COMPLEX="8TH3"     # AT118-H nanobody + AT1R (b562/BRIL fusion) + anti-BRIL Fab, cryo-EM 3.0 A

export NB_CHAIN="A"           # one nanobody copy is enough for a free-in-water study

# --- force field / solvent -------------------------------------------------
export FF="${FF:-amber99sb-ildn}"
export WATER="tip3p"
export BOX_PAD_NM="${BOX_PAD_NM:-0.9}"   # solute-to-box padding (0.9 = POC-friendly)
export ION_CONC_M="0.15"      # NaCl
export PBC="xyz"

# --- MD protocol -----------------------------------------------------------
export TEMP_K="300"
export NVT_PS="200"           # position-restrained heat
export NPT_PS="500"           # position-restrained density equilibration
export PROD_NS="${PROD_NS:-50}"   # production per replicate
export DT_PS="0.002"

# --- replicates: identical protocol, different starting velocities ---------
# Each replicate uses a distinct GROMACS gen_seed, i.e. independent Maxwell-
# Boltzmann initial velocities from the same equilibrated structure.
export N_REPLICATES="${N_REPLICATES:-3}"
export SEEDS=(101 202 303 404 505)   # one per replicate (override with -s if needed)

# --- gmx_MMPBSA ------------------------------------------------------------
export GB_MODEL="5"           # igb=5 (OBC-II), standard for protein-protein
export SALTCON="0.150"
export MMPBSA_START_FRAME="1" # after equil, frames are already production-only
export MMPBSA_END_FRAME="999999"
export MMPBSA_INTERVAL="10"   # use every 10th frame (200 ps @ 2 fs, 1 ps save)

# --- membrane (8TH3) -------------------------------------------------------
export MEMB_LIPID="POPC"
export MEMB_FF="${MEMB_FF:-charmm36-jul2022}"

# --- helpers ---------------------------------------------------------------
log()  { printf '\033[1;34m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }
die()  { printf '\033[1;31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }

require() { for c in "$@"; do command -v "$c" >/dev/null 2>&1 || die "missing executable: $c"; done; }
