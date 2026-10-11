#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 04_run_mmgbsa.sh -- per-residue MM-GBSA decomposition of the nanobody/AT1R
# binding free energy (single-trajectory) from the 8TH3 production run.
#
# Inputs (from 06_run_membrane.sh / CHARMM-GUI):
#   $WORK/8TH3/prod/md.tpr, md.xtc, topol.top
#
# Notes
#  * gmx_MMPBSA -cg notation is <Receptor group> <Ligand group>; the complex
#    is their union, so no separate "complex" group is needed.
#  * The trajectory must be PBC-removed for solvation energies -- we do that
#    here with trjconv (center + -pbc mol) before the MM-PBSA step.
#  * OBC (igb=5) has no explicit bilayer term for a membrane-embedded GPCR:
#    take the ABSOLUTE dG as qualitative and read the per-residue
#    decomposition for "which nanobody residues drive binding".
#
# Run:  bash simulation/04_run_mmgbsa.sh
# ---------------------------------------------------------------------------
source "$(dirname "$0")/config.sh"
require gmx_MMPBSA gmx

F="$WORK/8TH3/prod"
[ -f "$F/md.tpr" ]   || die "no $F/md.tpr -- run 06_run_membrane.sh first"
[ -f "$F/md.xtc" ]   || die "no $F/md.xtc"
[ -f "$F/topol.top" ] || die "no $F/topol.top"
cd "$F"

# ---- PBC-correct the trajectory (center + mol-wise minimum image) ----------
if [ ! -f md_center.xtc ]; then
  log "centering / removing pbc (md.xtc -> md_center.xtc)"
  gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc \
        -center -pbc mol <<< "Protein" > trjconv.log 2>&1 || die "trjconv pbc failed"
fi

# ---- index groups: Receptor (chain A) and Nanobody (chain B) ---------------
if [ ! -f index_mmgbsa.ndx ]; then
  log "building Receptor / Nanobody groups with gmx select"
  gmx select -s md.tpr -on index_mmgbsa.ndx \
        <<< $'chain A and not resname POPC and not resname SOL and not ion\nname 1 Receptor\nchain B and not resname POPC and not resname SOL and not ion\nname 2 Nanobody\nq' \
        > sel_mmgbsa.log 2>&1 \
      || die "gmx select failed; check which chains the CHARMM-GUI system uses"
fi

# ---- run gmx_MMPBSA --------------------------------------------------------
log "running gmx_MMPBSA (igb=${GB_MODEL}, salt=0.15 M)"
gmx_MMPBSA -O -i "$SIM_ROOT/04_mmgbsa/mmpbsa.in" \
    -cs md.tpr -ct md_center.xtc -cp topol.top \
    -ci index_mmgbsa.ndx -cg 1 2 \
    -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv \
    > gmx_MMPBSA.log 2>&1 || die "gmx_MMPBSA failed; see gmx_MMPBSA.log"

log "binding energetics:"
awk '/DELTA TOTAL|DELTA G gas|DELTA G solv/{print}' FINAL_RESULTS_MMPBSA.dat | head -6
log "per-residue decomposition in FINAL_RESULTS_MMPBSA.csv"