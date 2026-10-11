#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02_run_replicates.sh -- restrained equilibration + production MD for every
# replicate of the free-nanobody systems (7T83, 7T84).
#
# "Replicates" here = identical protocol but a distinct GROMACS gen_seed in
# the NVT step, i.e. independent Maxwell-Boltzmann starting velocities drawn
# from the SAME equilibrated structure.  The spread across replicates is the
# sampling uncertainty reported by 03_analyze.py.
#
# Equilibration:  em -> nvt(restrained, gen_vel) -> npt(restrained)
# Production:     md (unrestrained, PROD_NS nanoseconds)
# Analysis:       rmsd / residue-rmsf / gyrate / hbond .xvg per replicate.
#
# Run:
#   bash simulation/02_run_replicates.sh                    # 3 replicates x 50 ns
#   PROD_NS=20 N_REPLICATES=2 bash simulation/02_run_replicates.sh
# ---------------------------------------------------------------------------
source "$(dirname "$0")/config.sh"
require gmx

NSTEPS_MD="${NSTEPS_MD:-$(( PROD_NS * 1000000 ))}"   # 2 fs steps -> (NSTEPS_MD*2) ps

make_nvt_mdp() { # $1=seed -> writes nvt_${1}.mdp
  sed "s/^gen_seed         = .*/gen_seed         = ${1}/" "$SIM_ROOT/mdp/nvt.mdp" > "${1}_nvt.mdp"
}

run_replicate() {
  local SYSTEM="$1" N="$2" SEED="$3"
  local P="$WORK/${SYSTEM}/step1_padding"
  local R="$WORK/${SYSTEM}/rep${N}_seed${SEED}"
  mkdir -p "$R"
  cd "$R"
  log "[${SYSTEM} rep${N}] seed=${SEED}  starting in $(pwd)"

  # ---- energy minimization -------------------------------------------------
  gmx grompp -f "$SIM_ROOT/mdp/em.mdp" -c "$P/step1.gro" -p "$P/top.top" \
        -o em.tpr > grompp_em.log 2>&1 \
        || die "grompp em failed (${SYSTEM} rep${N})"
  gmx mdrun -deffnm em ${GMX_ORG} > mdrun_em.log 2>&1 \
        || die "mdrun em failed (${SYSTEM} rep${N})"

  # ---- NVT: distinct starting velocities via gen_seed ----------------------
  make_nvt_mdp "$SEED"
  NVT_STEPS="${NVT_STEPS:-100000}"
  sed "s/^nsteps           = .*/nsteps           = ${NVT_STEPS}/" "./${SEED}_nvt.mdp" > nvt_used.mdp
  gmx grompp -f nvt_used.mdp -c em.gro -r em.gro -p "$P/top.top" \
        -n "$P/index.ndx" -o nvt.tpr > grompp_nvt.log 2>&1 \
        || die "grompp nvt failed (${SYSTEM} rep${N})"
  gmx mdrun -deffnm nvt ${GMX_ORG} > mdrun_nvt.log 2>&1 \
        || die "mdrun nvt failed (${SYSTEM} rep${N})"

  # ---- NPT: density equilibration (restrained, velocities carried over) ----
  NPT_STEPS="${NPT_STEPS:-250000}"
  sed "s/^nsteps           = .*/nsteps           = ${NPT_STEPS}/" \
      "$SIM_ROOT/mdp/npt.mdp" > npt_used.mdp
  gmx grompp -f npt_used.mdp -c nvt.gro -r nvt.gro -p "$P/top.top" \
        -n "$P/index.ndx" -o npt.tpr > grompp_npt.log 2>&1 \
        || die "grompp npt failed (${SYSTEM} rep${N})"
  gmx mdrun -deffnm npt ${GMX_ORG} > mdrun_npt.log 2>&1 \
        || die "mdrun npt failed (${SYSTEM} rep${N})"

  # ---- production ----------------------------------------------------------
  sed "s/^nsteps           = .*/nsteps           = ${NSTEPS_MD}/" \
      "$SIM_ROOT/mdp/md.mdp" > md_prod.mdp
  gmx grompp -f md_prod.mdp -c npt.gro -p "$P/top.top" -n "$P/index.ndx" \
        -o md.tpr > grompp_md.log 2>&1 \
        || die "grompp md failed (${SYSTEM} rep${N})"
  gmx mdrun -deffnm md ${GMX_ORG} > mdrun_md.log 2>&1 \
        || die "mdrun md failed (${SYSTEM} rep${N})"

  # ---- per-replicate observables (xvg for 03_analyze.py) -------------------
  gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns \
          <<< "Backbone" > rms.log 2>&1
  gmx rmsf -s md.tpr -f md.xtc -o rmsf_bb.xvg -res -fit \
          <<< "Backbone" >> rms.log 2>&1
  gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg -tu ns \
          <<< "Protein" >> rms.log 2>&1
  gmx hbond -s md.tpr -f md.xtc -num hbond.xvg \
          <<< "Protein" >> rms.log 2>&1
  log "[${SYSTEM} rep${N}] MD done; xvg observables written"
  cd - > /dev/null
}

N=0
for S in "${PDB_WT}" "${PDB_VAR}"; do
  for SEED in "${SEEDS[@]:0:${N_REPLICATES}}"; do
    N=$((N+1))
    run_replicate "$S" "$N" "$SEED"
  done
done
log "replicates done: ${N} production runs."