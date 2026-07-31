#!/bin/bash
# ============================================================================
#  deploy_compare.sh -- cas de comparaison : bulle en domaine BORNE (fond +
#  surface libre) avec reduced.h, memes parametres, laminaire, lvl7, 4 coeurs.
#  A confronter au cas periodique + force volumique (lam7_frame).
#
#  Usage :  bash deploy_compare.sh            # DRY-RUN
#           LAUNCH=1 bash deploy_compare.sh   # soumet (1 run, 4c, shared-cpu)
# ============================================================================
set -euo pipefail
ROOT=/work/p0910/bonjeanf/bulles
D=$ROOT/compare_wall
HERE=$(cd "$(dirname "$0")" && pwd)
export BASILISK=/users/p0910/bonjeanf/basilisk2/basilisk/src
export PATH="$BASILISK:$PATH"

[ -f "$HERE/compare_wall.c" ] || { echo "ERREUR: compare_wall.c absent (push)"; exit 1; }
grep -q '#include "reduced.h"' "$HERE/compare_wall.c" || { echo "ERREUR: reduced.h attendu ici !"; exit 1; }

echo "== compilation (reduced.h, domaine borne) =="
cd "$HERE"
CC99='mpicc -std=c99' qcc -D_MPI=1 -O3 -Wall -D_GNU_SOURCE -disable-dimensions \
    compare_wall.c -o compare_wall -lm
echo "[ok] compare_wall"

rm -rf "$D"; mkdir -p "$D"
cp -f "$HERE/compare_wall" "$D/"
{
  echo "#!/bin/bash"
  echo "#SBATCH --job-name=cmpwall"
  echo "#SBATCH --ntasks=4"
  echo "#SBATCH --partition=shared-cpu"
  echo "#SBATCH --time=12:00:00"
  echo "#SBATCH --output=slurm-%j.out"
  echo "cd \"$D\""
  echo "srun ./compare_wall 7 18"     # maxlevel 7, MAXTIME=18, domaine L0=180
} > "$D/job.slurm"
echo "[prep] bulle bornee (L0=180) + reduced.h : lvl7, laminaire, g=4, Bo=1, 4c -> $D"

if [ "${LAUNCH:-0}" = 1 ]; then
  ( cd "$D" && sbatch job.slurm )
  echo; squeue -u "$USER" -o "%.8i %.10P %.10j %.2t %.10M %R" | grep -E "cmpwall|JOBID" || true
else
  echo "DRY-RUN. LAUNCH=1 bash deploy_compare.sh pour soumettre."
fi
