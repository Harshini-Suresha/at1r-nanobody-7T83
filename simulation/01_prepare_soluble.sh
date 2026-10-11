#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_prepare_soluble.sh -- build solvated, ionized GROMACS systems for the
# free nanobody variants 7T83 (WT) and 7T84 (G26D T57I).
#
# Pipeline per system:  pdb2gmx -> editconf (cubic box) -> solvate -> genion
#                          -> make_ndx (Protein / Non-Protein groups)
#
# Outputs in $WORK/<SYSTEM>/step1_padding/:  top.top, step1.gro, index.ndx ...
# Run:  bash simulation/01_prepare_soluble.sh
# ---------------------------------------------------------------------------
source "$(dirname "$0")/config.sh"
require gmx

prepare_system() {
  local SYSTEM="$1"
  local IN="$DATA/${SYSTEM}_nb.pdb"
  local W="$WORK/${SYSTEM}"; local P="$W/step1_padding"
  mkdir -p "$P"
  log "preparing ${SYSTEM} from ${IN}"

  # topology / structure with H atoms (non-interactive, accept defaults).
  # Here-string (not a pipe) so set -o pipefail never sees a SIGPIPE'ed yes.
  gmx pdb2gmx -f "$IN" -o "$P/1_prot.gro" -p "$P/top.top" \
        -i "$P/posre.itp" -ff "${FF}" -water "${WATER}" -ignh \
        <<< "0" > "$W/pdb2gmx.log" 2>&1 \
        || die "pdb2gmx failed for ${SYSTEM}; see $W/pdb2gmx.log"

  # cubic box with padding on every side
  gmx editconf -f "$P/1_prot.gro" -o "$P/2_box.gro" -c -d "${BOX_PAD_NM}" -bt cubic \
      > "$W/editconf.log" 2>&1 || die "editconf failed for ${SYSTEM}"

  # solvate
  gmx solvate -cp "$P/2_box.gro" -cs spc216.gro -o "$P/3_solv.gro" -p "$P/top.top" \
      > "$W/solvate.log" 2>&1 || die "solvate failed for ${SYSTEM}"

  # ions: grompp needs a tpr; em.mdp is fine.
  # (-maxwarn 1: Ewald net-charge warning is expected here -- genion -neutral
  #  neutralises the box immediately afterwards.)
  gmx grompp -f "$SIM_ROOT/mdp/em.mdp" -c "$P/3_solv.gro" -p "$P/top.top" \
        -o "$P/ions.tpr" -pp "$P/ions_only.top" -maxwarn 1 \
        > "$W/grompp_ions.log" 2>&1 \
        || die "grompp (ions) failed for ${SYSTEM}"

  # replace some SOL with Na+/Cl- to 0.15 M (and neutralise)
  gmx genion -s "$P/ions.tpr" -o "$P/4_ion.gro" -p "$P/ions_only.top" \
        -pname NA -nname CL -neutral -conc "${ION_CONC_M}" \
        <<< "SOL" > "$W/genion.log" 2>&1 || die "genion failed for ${SYSTEM}"
  mv "$P/ions_only.top" "$P/top.top"

  # grompp -pp (above) expands/skips the #ifdef POSRES block away, so re-add
  # the guarded include of posre.itp inside the protein moleculetype (before
  # its first [ bonds ]) on the FINAL topology.  An absolute path keeps grompp
  # able to resolve it no matter the run directory.
  /usr/bin/env python3 - "$P/top.top" "$P/posre.itp" <<'PY'
import sys
p, itp = sys.argv[1], sys.argv[2]
s = open(p).read()
if "#ifdef POSRES" in s:
    sys.exit(0)                      # already patched
marker = "\n[ bonds ]"
i = s.find(marker)
if i < 0:
    sys.exit("could not find protein [ bonds ] in " + p)
blk = ("\n; Include Position restraint file"
       "\n#ifdef POSRES\n#include \"%s\"\n#endif" % itp)
s = s[:i] + blk + s[i:]
open(p, "w").write(s)
PY
  [ -s "$P/top.top" ] || die "topology patch failed for ${SYSTEM}"

  # index groups: 0 System, 1 Protein, ... 11 non-Protein (defaults).  Rename
  # group 11 -> "Non-Protein" so the mdps' tc-grps = Protein Non-Protein works.
  gmx make_ndx -f "$P/4_ion.gro" -o "$P/index.ndx" \
        <<< $'name 11 Non-Protein\nq' \
        > "$W/make_ndx.log" 2>&1 || die "make_ndx failed for ${SYSTEM}"

  # final step-1 coordinate renamed for the run scripts
  cp "$P/4_ion.gro" "$P/step1.gro"
  log "  -> ${P}/top.top, step1.gro, index.ndx, posre.itp"
}

for S in "${PDB_WT}" "${PDB_VAR}"; do
  prepare_system "$S"
done
log "done."