#!/bin/bash
# ============================================================================
#  run_pipeline.sh -- Lancement et orchestration de la chaîne de calcul
# ============================================================================
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)

NTASKS=${NTASKS:-64}
NTASKS_PRE=${NTASKS_PRE:-$NTASKS}
NTASKS_L7=${NTASKS_L7:-$NTASKS}
NTASKS_L8=${NTASKS_L8:-$NTASKS}
NTASKS_L9=${NTASKS_L9:-$NTASKS}
PARTITION=${PARTITION:-compute}
ACCOUNT=${ACCOUNT:-}
TIME_PRE=${TIME_PRE:-12:00:00}
TIME_BUB=${TIME_BUB:-24:00:00}
MAXTIME_PRE=${MAXTIME_PRE:-40}
MAXTIME_BUB=${MAXTIME_BUB:-60}
R0=${R0:-8.0}
FORCED=${FORCED:-1}
KE_TARGET=${KE_TARGET:-24}
PRE_LEVEL=${PRE_LEVEL:-8}
LEVELS=${LEVELS:-"7 8 9"}
LOCAL=${LOCAL:-0}
EXE="$ROOT/bubble"
PRE="$ROOT/simulations/precursor"

acct_line=""
[ -n "$ACCOUNT" ] && acct_line="#SBATCH --account=$ACCOUNT"

ntasks_for () {
  case "$1" in
    7) echo "$NTASKS_L7" ;;
    8) echo "$NTASKS_L8" ;;
    9) echo "$NTASKS_L9" ;;
    *) echo "$NTASKS"    ;;
  esac
}

echo "[build] make"
make -C "$ROOT"

mkdir -p "$PRE"
cp "$EXE" "$PRE/"

if [ "$LOCAL" = "1" ]; then
  echo "[local] précurseur (lvl $PRE_LEVEL, $NTASKS_PRE cœurs)"
  ( cd "$PRE" && mpirun -np "$NTASKS_PRE" ./bubble \
        "$PRE_LEVEL" "$MAXTIME_PRE" "$R0" "$FORCED" "$KE_TARGET" 0 \
        >log.out 2>log.err )

  for lvl in $LEVELS; do
    W="$ROOT/simulations/convturb/lvl$lvl"
    NT=$(ntasks_for "$lvl")
    mkdir -p "$W"
    cp "$EXE" "$W/"
    cp "$PRE/end" "$W/restart"
    echo "[local] bulle lvl$lvl ($NT cœurs)"
    ( cd "$W" && mpirun -np "$NT" ./bubble \
          "$lvl" "$MAXTIME_BUB" "$R0" "$FORCED" "$KE_TARGET" 1 \
          >log.out 2>log.err )
  done
  echo "[local] terminé."
  exit 0
fi

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
echo "[slurm] précurseur soumis : job $PRE_ID"

for lvl in $LEVELS; do
  W="$ROOT/simulations/convturb/lvl$lvl"
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
  echo "[slurm] bulle lvl$lvl soumise : job $ID"
done

echo "[slurm] tout est soumis."
