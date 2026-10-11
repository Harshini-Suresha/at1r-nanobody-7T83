#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 06_run_membrane.sh -- production MD + analysis for the 8TH3 nanobody/AT1R
# membrane system.
#
# Expects $WORK/8TH3/prod/ to contain the CHARMM-GUI GROMACS system staged by
# 05_prepare_membrane.sh (step5_assembly.gro/.top, index.ndx, toppar/, *.itp).
#
#  * If the CHARMM-GUI step6 equilibration already ran, production continues
#    from step6.6_equilibration; otherwise we run a self-contained
#    em -> nvt -> npt (restrained) before production.
#  * Parameters are read from config.sh (PROD_NS, TEMP_K=310 uses membrane
#    mdps; T in the CHARMM-GUI system may differ -- set TEMP_K to match).
#
# Run:  bash simulation/06_run_membrane.sh
# ---------------------------------------------------------------------------
source "$(dirname "$0")/config.sh"
require gmx

F="$WORK/8TH3/prod"
[ -f "$F/step5_assembly.gro" ] || die "no staged system -- run 05_prepare_membrane.sh first"

cd "$F"
# membrane simulations run warmer; CHARMM-GUI systems are pre-equilibrated at
# the temperature set during the build (default 310 K)
TMP_REF="${TEMP_K:-310}"

# ---- continue from CHARMM-GUI equilibration if it is present --------------
EQ=""
for cand in step6.6_equilibration/step6.6_equilibration.gro \
            step6.7_equilibration/step6.7_equilibration.gro \
            step6_equilibration.gro; do
  [ -f "$cand" ] && EQ="$cand" && break
done

if [ -n "$EQ" ]; then
  log "continuing production from CHARMM-GUI equilibration: $EQ"
  START_FILE="$EQ"
else
  log "no CHARMM-GUI step6 output found -- running built-in em/nvt/npt"
  gmx grompp -f "$SIM_ROOT/mdp/em_mem.mdp" -c step5_assembly.gro -p topol.top \
        -o em_mem.tpr > grompp_em.log 2>&1 || die "grompp em failed"
  gmx mdrun -deffnm em_mem ${GMX_ORG} > mdrun_em.log 2>&1 || die "mdrun em failed"

  sed -e "s/^gen_temp         = .*/gen_temp         = ${TMP_REF}/" \
      -e "s/^ref_t            = 310 310 310/ref_t            = ${TMP_REF} ${TMP_REF} ${TMP_REF}/" \
      "$SIM_ROOT/mdp/nvt_mem.mdp" > nvt_mem.mdp
  gmx grompp -f nvt_mem.mdp -c em_mem.gro -r em_mem.gro -p topol.top \
        -n index.ndx -o nvt_mem.tpr > grompp_nvt.log 2>&1 || die "grompp nvt failed"
  gmx mdrun -deffnm nvt_mem ${GMX_ORG} > mdrun_nvt.log 2>&1 || die "mdrun nvt failed"

  sed -e "s/^ref_t            = 310 310 310/ref_t            = ${TMP_REF} ${TMP_REF} ${TMP_REF}/" \
      "$SIM_ROOT/mdp/npt_mem.mdp" > npt_mem.mdp
  gmx grompp -f npt_mem.mdp -c nvt_mem.gro -r nvt_mem.gro -p topol.top \
        -n index.ndx -o npt_mem.tpr > grompp_npt.log 2>&1 || die "grompp npt failed"
  gmx mdrun -deffnm npt_mem ${GMX_ORG} > mdrun_npt.log 2>&1 || die "mdrun npt failed"
  START_FILE=npt_mem.gro
fi

# ---- production ------------------------------------------------------------
NSTEPS=$(( PROD_NS * 1000000 ))
sed "s/^nsteps           = .*/nsteps           = ${NSTEPS}/" \
    "$SIM_ROOT/mdp/md_mem.mdp" > md_prod.mdp
gmx grompp -f md_prod.mdp -c "${START_FILE}" -p topol.top -n index.ndx \
      -o md.tpr > grompp_md.log 2>&1 || die "grompp md failed"
gmx mdrun -deffnm md ${GMX_ORG} > mdrun_md.log 2>&1 || die "mdrun md failed"

# ---- analysis (receptor RMSD, per-chain residue RMSF, Rg, H-bonds) ---------
gmx rms -s md.tpr -f md.xtc -o rmsd_prot.xvg -tu ns <<< "Protein" > anal.log 2>&1
gmx rmsf -s md.tpr -f md.xtc -o rmsf_bb.xvg -res -fit <<< "Backbone" >> anal.log 2>&1
gmx gyrate -s md.tpr -f md.xtc -o gyrate_prot.xvg -tu ns <<< "Protein" >> anal.log 2>&1
gmx hbond -s md.tpr -f md.xtc -num hbond_popc.xvg <<< "POPC" >> anal.log 2>&1

log "membrane MD done: md.xtc, rmsd_prot.xvg, rmsf_bb.xvg, gyrate_prot.xvg"
log "next: 07_analyze_interface.py (nanobody-on-receptor contacts) and"
log "       04_run_mmgbsa.sh (binding energetics)."