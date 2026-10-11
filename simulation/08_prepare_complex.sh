#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 08_prepare_complex.sh -- build the solvated, ionized 8TH3 nanobody-AT1R
# complex (chains A = receptor, B = nanobody) ready for short production MD.
#
# The refit pose (data/8TH3_complex.pdb) is used as-is.  pdb2gmx writes
# posre.itp for chain A only, and 01's topology patch is NOT applied here:
# instead the run script adds a receptor-only position-restraint block, so the
# nanobody stays free while the receptor is held at the crystallographic pose.
#
# Outputs in $WORK/8TH3/step1_padding/:  top.top, step1.gro, index.ndx, posre.itp
# Run:  bash simulation/08_prepare_complex.sh
# ---------------------------------------------------------------------------
source "$(dirname "$0")/config.sh"
require gmx

SYSTEM="8TH3"
IN="$DATA/${SYSTEM}_complex.pdb"
W="$WORK/${SYSTEM}"; P="$W/step1_padding"
mkdir -p "$P"
log "preparing ${SYSTEM} complex from ${IN}"

gmx pdb2gmx -f "$IN" -o "$P/1_prot.gro" -p "$P/top.top" \
      -i "$P/posre.itp" -ff "${FF}" -water "${WATER}" -ignh \
      <<< "0" > "$W/pdb2gmx.log" 2>&1 \
      || die "pdb2gmx failed for ${SYSTEM}; see $W/pdb2gmx.log"

gmx editconf -f "$P/1_prot.gro" -o "$P/2_box.gro" -c -d "${BOX_PAD_NM}" -bt cubic \
    > "$W/editconf.log" 2>&1 || die "editconf failed for ${SYSTEM}"

gmx solvate -cp "$P/2_box.gro" -cs spc216.gro -o "$P/3_solv.gro" -p "$P/top.top" \
    > "$W/solvate.log" 2>&1 || die "solvate failed for ${SYSTEM}"

gmx grompp -f "$SIM_ROOT/mdp/em.mdp" -c "$P/3_solv.gro" -p "$P/top.top" \
      -o "$P/ions.tpr" -pp "$P/ions_only.top" -maxwarn 1 \
      > "$W/grompp_ions.log" 2>&1 || die "grompp (ions) failed for ${SYSTEM}"

gmx genion -s "$P/ions.tpr" -o "$P/4_ion.gro" -p "$P/ions_only.top" \
      -pname NA -nname CL -neutral -conc "${ION_CONC_M}" \
      <<< "SOL" > "$W/genion.log" 2>&1 || die "genion failed for ${SYSTEM}"
mv "$P/ions_only.top" "$P/top.top"

# index groups: 0 System, 1 Protein, ... 11 non-Protein (defaults from gmx).
# Chain-based Receptor/Nanobody groups are built with gmx select at run time.
gmx make_ndx -f "$P/4_ion.gro" -o "$P/index.ndx" \
      <<< $'name 11 Non-Protein\nq' \
      > "$W/make_ndx.log" 2>&1 || die "make_ndx failed for ${SYSTEM}"

cp "$P/4_ion.gro" "$P/step1.gro"
log "  -> ${P}/top.top, step1.gro, index.ndx, posre.itp"
log "done."