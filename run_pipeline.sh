#!/bin/bash
# ============================================================================
#  run_pipeline.sh -- Orchestration complete de l'etude
# ----------------------------------------------------------------------------
#  1) compile le solveur
#  2) lance le PRECURSEUR (turbulence seule, INJECT=0) -> dump "end"
#  3) pour lvl = 7,8,9 : stage "end" -> "restart" puis lance la BULLE (INJECT=1)
#
#  Deux modes :
#    - SLURM (defaut)  : sbatch avec dependance afterok sur le precurseur.
#    - LOCAL=1         : execution directe via mpirun (poste/noeud interactif).
#
#  Variables surchargeable en ligne de commande, ex :
#    NTASKS=128 PARTITION=skylake MAXTIME_BUB=80 ./run_pipeline.sh
# ============================================================================
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)

# ----------------------------- Configuration --------------------------------
NTASKS=${NTASKS:-64}              # nb de coeurs MPI par defaut (fallback)
# Coeurs PAR NIVEAU : chaque niveau octree = ~x8 cellules en 3D, donc lvl9 >> lvl8 >> lvl7.
# Vise ~100k cellules/coeur (voir cells.dat). Si non defini -> retombe sur NTASKS.
NTASKS_PRE=${NTASKS_PRE:-$NTASKS} # precurseur (lvl 9 conseille -> beaucoup de coeurs)
NTASKS_L7=${NTASKS_L7:-$NTASKS}   # bulle lvl7
NTASKS_L8=${NTASKS_L8:-$NTASKS}   # bulle lvl8
NTASKS_L9=${NTASKS_L9:-$NTASKS}   # bulle lvl9
PARTITION=${PARTITION:-compute}   # partition SLURM
ACCOUNT=${ACCOUNT:-}              # compte projet (optionnel)
TIME_PRE=${TIME_PRE:-12:00:00}    # walltime precurseur
TIME_BUB=${TIME_BUB:-24:00:00}    # walltime runs bulle
MAXTIME_PRE=${MAXTIME_PRE:-40}    # temps physique precurseur
MAXTIME_BUB=${MAXTIME_BUB:-60}    # temps physique runs bulle
R0=${R0:-8.0}                     # rayon de la bulle
FORCED=${FORCED:-1}               # 1 = turbulence forcee, 0 = libre
KE_TARGET=${KE_TARGET:-24}        # energie cinetique cible (forcage controle) ; calibrer pour Re_lambda~60-80
PRE_LEVEL=${PRE_LEVEL:-8}         # niveau du precurseur UNIFORME (8=256^3 resout Re_l~70 ; 9=512^3=1.3e8 mailles, lourd)
LEVELS=${LEVELS:-"7 8 9"}         # niveaux pour l'etude de convergence
LOCAL=${LOCAL:-0}                 # 1 = mpirun direct au lieu de SLURM
EXE="$ROOT/bubble"
PRE="$ROOT/simulations/precursor"

acct_line=""
[ -n "$ACCOUNT" ] && acct_line="#SBATCH --account=$ACCOUNT"

# Coeurs alloues pour un niveau de bulle donne.
ntasks_for () {
  case "$1" in
    7) echo "$NTASKS_L7" ;;
    8) echo "$NTASKS_L8" ;;
    9) echo "$NTASKS_L9" ;;
    *) echo "$NTASKS"    ;;
  esac
}

# ------------------------------- Build --------------------------------------
echo "[build] make"
make -C "$ROOT"

mkdir -p "$PRE"
cp "$EXE" "$PRE/"

# =================================================================
#  Mode LOCAL : tout en sequentiel (precurseur puis les 3 niveaux)
# =================================================================
if [ "$LOCAL" = "1" ]; then
  echo "[local] precurseur (lvl $PRE_LEVEL, $NTASKS_PRE coeurs, t<=$MAXTIME_PRE)"
  ( cd "$PRE" && mpirun -np "$NTASKS_PRE" ./bubble \
        "$PRE_LEVEL" "$MAXTIME_PRE" "$R0" "$FORCED" "$KE_TARGET" 0 \
        >log.out 2>log.err )

  for lvl in $LEVELS; do
    W="$ROOT/simulations/lvl$lvl"
    NT=$(ntasks_for "$lvl")
    mkdir -p "$W"
    cp "$EXE" "$W/"
    cp "$PRE/end" "$W/restart"
    echo "[local] bulle lvl$lvl ($NT coeurs, t<=$MAXTIME_BUB)"
    ( cd "$W" && mpirun -np "$NT" ./bubble \
          "$lvl" "$MAXTIME_BUB" "$R0" "$FORCED" "$KE_TARGET" 1 \
          >log.out 2>log.err )
  done
  echo "[local] termine."
  exit 0
fi

# =================================================================
#  Mode SLURM : soumission avec dependances
# =================================================================
cat > "$PRE/job.slurm" <<EOF
#!/bin/bash
#SBATCH --job-name=precursor
#SBATCH --ntasks=$NTASKS_PRE
#SBATCH --partition=$PARTITION
#SBATCH --time=$TIME_PRE
$acct_line
#SBATCH --output=$PRE/slurm-%j.out
cd "$PRE"
srun ./bubble $PRE_LEVEL $MAXTIME_PRE $R0 $FORCED $KE_TARGET 0
EOF

PRE_ID=$(sbatch --parsable "$PRE/job.slurm")
echo "[slurm] precurseur soumis (lvl $PRE_LEVEL, $NTASKS_PRE coeurs) : job $PRE_ID"

for lvl in $LEVELS; do
  W="$ROOT/simulations/lvl$lvl"
  NT=$(ntasks_for "$lvl")
  mkdir -p "$W"
  cp "$EXE" "$W/"
  cat > "$W/job.slurm" <<EOF
#!/bin/bash
#SBATCH --job-name=lvl$lvl
#SBATCH --ntasks=$NT
#SBATCH --partition=$PARTITION
#SBATCH --time=$TIME_BUB
$acct_line
#SBATCH --output=$W/slurm-%j.out
cd "$W"
cp "$PRE/end" restart
srun ./bubble $lvl $MAXTIME_BUB $R0 $FORCED $KE_TARGET 1
EOF
  ID=$(sbatch --parsable --dependency=afterok:$PRE_ID "$W/job.slurm")
  echo "[slurm] bulle lvl$lvl soumise ($NT coeurs, apres $PRE_ID) : job $ID"
done

echo "[slurm] tout est soumis. Suivi : squeue -u \$USER"
