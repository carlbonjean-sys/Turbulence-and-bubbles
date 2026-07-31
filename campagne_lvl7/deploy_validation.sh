#!/bin/bash
# ============================================================================
#  CAMPAGNE DE VALIDATION -- 2026-07-22
#  (feu vert user : "convergence en maillage : alors fais la" + "courants
#   parasites : on peut faire ca aussi")
#
#  AUCUNE PRODUCTION ICI. Uniquement des cas dont la reponse est connue a
#  l'avance ou dont on veut la convergence.
#
#  V3  CONVERGENCE EN MAILLAGE du cas LAMINAIRE (le trou principal du dossier) :
#      meme run que lam7_frame (g=4, Bo=1, liquide au repos, repere mobile),
#      a maxlevel 6 et 8. Avec le lvl7 deja fait -> 3 points.
#      Enjeu : a lvl7 la couche limite fait 1.64 maille (delta ~ D/sqrt(Re_b)).
#      On veut savoir de combien u_inf = 12.27 bouge.
#        lvl6 : D/dx =  8.5   lvl7 : D/dx = 17.1   lvl8 : D/dx = 34.1
#
#  V1  COURANTS PARASITES (src/static_bubble.c) : bulle spherique au repos,
#      sans gravite -> la vitesse exacte est ZERO. Tout ce qu'on mesure est de
#      l'erreur du schema de tension de surface. A lvl 6, 7 et 8.
#      Jamais mesure dans ce projet.
#
#  Le binaire de production est recompile depuis main.c CORRIGE (event bubble
#  reecrit le 2026-07-22 : sommes completes ponderees (1-f), colonnes 14-16).
#
#  Usage :  bash deploy_validation.sh            -> DRY-RUN (montre, ne soumet pas)
#           LAUNCH=1 bash deploy_validation.sh   -> soumet
# ============================================================================
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT=/work/p0910/bonjeanf/bulles
VAL="$ROOT/validation"
PART=${PART:-shared-cpu}

# --- garde-fous ---------------------------------------------------------
[ -f "$HERE/main.c" ]          || { echo "ERREUR: main.c absent";          exit 1; }
[ -f "$HERE/static_bubble.c" ] || { echo "ERREUR: static_bubble.c absent"; exit 1; }
[ -f "$ROOT/designA/we07_pre/end" ] || { echo "ERREUR: precurseur we07_pre/end absent"; exit 1; }
grep -q '#include "reduced.h"' "$HERE/main.c" && { echo "ERREUR: reduced.h de retour !"; exit 1; }
grep -q 'event move_frame' "$HERE/main.c" || { echo "ERREUR: repere mobile absent"; exit 1; }
grep -q 'n_regions V_tag straddle' "$HERE/main.c" || { echo "ERREUR: main.c n'est PAS la version corrigee du 2026-07-22"; exit 1; }

echo "== compilation =="
cd "$HERE"
CC99='mpicc -std=c99' qcc -D_MPI=1 -O3 -Wall -D_GNU_SOURCE -disable-dimensions \
    main.c -o bubble -lm
echo "[ok] bubble (main.c corrige)"
CC99='mpicc -std=c99' qcc -D_MPI=1 -O3 -Wall -D_GNU_SOURCE -disable-dimensions \
    static_bubble.c -o static_bubble -lm
echo "[ok] static_bubble"

mkdir -p "$VAL"

# ---------------------------------------------------------------- V3 : lvl6 / lvl8
#   lam7_frame a tourne t=160->260 sur 4 coeurs en ~2h20.
#   lvl8 : ~8x plus de mailles et dt/2.8 (limite capillaire) => ~23x le cout.
#          On raccourcit a t=160->250 (90 u.t.) : transitoire ~10 + plateau ~30.
prep_lam () {           # $1=level  $2=maxtime  $3=ntasks  $4=walltime
  local D="$VAL/lam$1"
  rm -rf "$D"; mkdir -p "$D"
  cp -f "$HERE/bubble" "$D/"
  cp "$ROOT/designA/we07_pre/end" "$D/restart"
  {
    echo "#!/bin/bash"
    echo "#SBATCH --job-name=lam$1"
    echo "#SBATCH --ntasks=$3"
    echo "#SBATCH --partition=$PART"
    echo "#SBATCH --time=$4"
    echo "#SBATCH --output=slurm-%j.out"
    echo "cd \"$D\""
    #                        maxlvl MAXTIME R0 FORCED KE INJ BOND OMECO GRAV
    echo "srun ./bubble $1 $2 8.0 0 0 1 1.0 0.20 4"
  } > "$D/job.slurm"
  echo "[prep] V3 laminaire lvl$1 : t=160->$2, $3 coeurs, $4  -> $D"
}
prep_lam 6 260 2  12:00:00
prep_lam 8 250 16 48:00:00

# ---------------------------------------------------------------- V1 : parasites
prep_stat () {          # $1=level  $2=ntasks  $3=walltime
  local D="$VAL/spur$1"
  rm -rf "$D"; mkdir -p "$D"
  cp -f "$HERE/static_bubble" "$D/"
  {
    echo "#!/bin/bash"
    echo "#SBATCH --job-name=spur$1"
    echo "#SBATCH --ntasks=$2"
    echo "#SBATCH --partition=$PART"
    echo "#SBATCH --time=$3"
    echo "#SBATCH --output=slurm-%j.out"
    echo "cd \"$D\""
    echo "srun ./static_bubble $1 5.0"
  } > "$D/job.slurm"
  echo "[prep] V1 courants parasites lvl$1 : t=0->5, $2 coeurs, $3  -> $D"
}
prep_stat 6 1 02:00:00
prep_stat 7 2 04:00:00
prep_stat 8 4 08:00:00

TOT=$((2 + 16 + 1 + 2 + 4))
echo
echo "TOTAL : 5 jobs, $TOT coeurs (garde-fou 64), walltime max 48h (garde-fou 72h)"
echo

if [ "${LAUNCH:-0}" = 1 ]; then
  for D in "$VAL/lam6" "$VAL/lam8" "$VAL/spur6" "$VAL/spur7" "$VAL/spur8"; do
    ( cd "$D" && sbatch job.slurm )
  done
  echo; squeue -u "$USER" -o "%.9i %.12j %.10P %.9T %.11M %.5D %R"
else
  echo "DRY-RUN. Relancer avec LAUNCH=1 pour soumettre."
fi
