#!/bin/bash
# ============================================================================
#  deploy_ensemble.sh -- MOYENNE D'ENSEMBLE pour des points FIABLES.
#
#  Constat (campagne weber_fix, 1 realisation/point) : a beta>=0.3 le scatter
#  realisation-unique est enorme (SEM +/-11-15% ; wt21 -39% vs wt25 -24% a beta
#  quasi egal). Cause : u'/u_inf ~ 0.5 -> la bulle est ballottee, une seule
#  trajectoire ne fixe pas la moyenne. Remede : plusieurs injections a des
#  positions DECALEES dans le MEME precurseur (turbulence homogene -> realisations
#  independantes), puis moyenne. L'erreur sur la moyenne ~ SEM/sqrt(N_membres).
#
#  Ici : MEMBER=1 et MEMBER=2 (2 nouvelles positions) pour CHACUN des 5 Weber.
#  Avec le run existant (MEMBER=0 dans weber_fix/) -> 3 membres/point.
#  10 runs x 4 coeurs = 40 <= garde-fou 64. Duree t=160->280 (120 u.t. : plateau
#  sur [220,280] ~ 4 tau, suffisant par membre car c'est l'ENSEMBLE qui resserre).
#
#  Usage :  bash deploy_ensemble.sh            # DRY-RUN
#           LAUNCH=1 bash deploy_ensemble.sh   # soumet (10 runs, shared-cpu)
# ============================================================================
set -euo pipefail
ROOT=/work/p0910/bonjeanf/bulles
BASE=$ROOT/ensemble
PREBASE=$ROOT/designA
HERE=$(cd "$(dirname "$0")" && pwd)
export BASILISK=/users/p0910/bonjeanf/basilisk2/basilisk/src
export PATH="$BASILISK:$PATH"
GRAV0=${GRAV0:-4}
TMAX=${TMAX:-280}
NT=${NT:-4}
PART=${PART:-shared-cpu}
MEMBERS=${MEMBERS:-"1 2"}     # membres supplementaires (0 = deja fait dans weber_fix)

[ -f "$HERE/main.c" ] || { echo "ERREUR: main.c absent"; exit 1; }
grep -q 'MEMBER' "$HERE/main.c" || { echo "ERREUR: support MEMBER absent de main.c"; exit 1; }
grep -q 'event move_frame' "$HERE/main.c" || { echo "ERREUR: repere mobile absent"; exit 1; }

echo "== compilation =="
cd "$HERE"
CC99='mpicc -std=c99' qcc -D_MPI=1 -O3 -Wall -D_GNU_SOURCE -disable-dimensions \
    main.c -o bubble -lm
echo "[ok] bubble"

LABS=(wt05 wt11 wt21 wt25 wt32)
PRES=(we03 we07 we13 we16 we20)
KES=( 5    11   21   25   32)

mkdir -p "$BASE"
JOBS=""
for i in "${!LABS[@]}"; do
  pre=${PRES[$i]}; ke=${KES[$i]}
  [ -f "$PREBASE/${pre}_pre/end" ] || { echo "ERREUR: ${pre}_pre/end absent"; exit 1; }
  for m in $MEMBERS; do
    lab=${LABS[$i]}_g${GRAV0}_m${m}
    D=$BASE/${lab}_bub
    rm -rf "$D"; mkdir -p "$D"
    cp -f "$HERE/bubble" "$D/"
    cp    "$PREBASE/${pre}_pre/end" "$D/restart"
    { echo "#!/bin/bash"
      echo "#SBATCH --job-name=$lab"
      echo "#SBATCH --ntasks=$NT"
      echo "#SBATCH --partition=$PART"
      echo "#SBATCH --time=48:00:00"
      echo "#SBATCH --output=slurm-%j.out"
      echo "cd \"$D\""
      #                    maxlvl MAXTIME R0 FORCED KE INJ BOND OMECO GRAV NU MEMBER
      echo "srun ./bubble 7 $TMAX 8.0 1 $ke 1 1.0 0.20 $GRAV0 -1 $m"
    } > "$D/job.slurm"
    echo "[prep] $lab : KE=$ke (${pre}) GRAV=$GRAV0 MEMBER=$m -> $D"
    JOBS="$JOBS $D"
  done
done

NR=$(echo $JOBS | wc -w)
echo
echo "Recap : $NR runs ($NT coeurs, $PART, t=160->$TMAX), membres [$MEMBERS] x 5 Weber."
echo "  -> avec le membre 0 (weber_fix), 3 mesures/point ; erreur ~ SEM/sqrt(3)=0.58x."

if [ "${LAUNCH:-0}" = 1 ]; then
  for d in $JOBS; do ( cd "$d" && sbatch job.slurm ); done
  echo; squeue -u "$USER" -o "%.8i %.10P %.12j %.2t %.10M %R" | head -20
else
  echo "DRY-RUN : rien soumis. LAUNCH=1 bash deploy_ensemble.sh pour soumettre."
fi
