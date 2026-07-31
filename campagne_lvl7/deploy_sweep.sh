#!/bin/bash
# ============================================================================
#  deploy_sweep.sh -- BALAYAGE β DENSE : 9 nouveaux points d'intensite turbulente
#  pour completer la courbe u_inf_turb/u_inf(β), avec 5 membres d'ensemble chacun.
#
#  β = u'/u_inf visés : 0.05 0.10 0.18 0.20 0.26 0.29 0.35 0.40 0.45
#  KE_TARGET = 1.5*(β*u_inf)^2 avec u_inf=12.24 (laminaire mesure).
#
#  Pipeline 2 etages par point (comme Design A) :
#   1) PRECURSEUR (INJECT=0, monophasique) : reprend end_we1, re-equilibre au KE
#      cible -> end (t=160). Rapide (dt monophasique). 4 coeurs.
#   2) 5 BULLES (INJECT=1, MEMBER=0..4) : injectees dans CE precurseur a positions
#      decalees (realisations independantes), g=4 Bo=1 lvl7, t=160->280. afterok
#      precurseur. 2 coeurs (econome), shared-cpu.
#
#  cout : 9 pre + 45 bulle = 54 jobs. Kairos ne lance que ~4-6 jobs a la fois
#  (AssocGrpJobsLimit) -> la file met du temps (jours-semaines). Les dependances
#  afterok gerent l'ordre : rien a surveiller.
#
#  Usage :  bash deploy_sweep.sh            # DRY-RUN
#           LAUNCH=1 bash deploy_sweep.sh   # soumet (precurseurs + bulles dependantes)
# ============================================================================
set -euo pipefail
ROOT=${ROOT:-/work/p0910/bonjeanf/bulles}
BASE=$ROOT/sweep
HERE=$(cd "$(dirname "$0")" && pwd)
export BASILISK=${BASILISK:-/users/p0910/bonjeanf/basilisk2/basilisk/src}
export PATH="$BASILISK:$PATH"
PART=${PART:-shared-cpu}
NT_PRE=${NT_PRE:-4}
NT_BUB=${NT_BUB:-2}
TPRE=${TPRE:-160}
TBUB=${TBUB:-280}
GRAV0=${GRAV0:-4}
NMEMBERS=${NMEMBERS:-5}     # membres par point (m0..m{N-1})

[ -f "$ROOT/end_we1" ] || { echo "ERREUR: $ROOT/end_we1 absent"; exit 1; }
[ -f "$HERE/main.c" ]  || { echo "ERREUR: main.c absent"; exit 1; }
grep -q 'MEMBER' "$HERE/main.c" || { echo "ERREUR: support MEMBER absent"; exit 1; }

# points : label / KE (= 1.5*(β*12.24)^2)
LABS=(sw05 sw10 sw18 sw20 sw26 sw29 sw35 sw40 sw45)
BETA=(0.05 0.10 0.18 0.20 0.26 0.29 0.35 0.40 0.45)
KES=( 0.56 2.25 7.28 8.99 15.19 18.9 27.5 36.0 45.5)

echo "== compilation =="
cd "$HERE"
CC99='mpicc -std=c99' qcc -D_MPI=1 -O3 -Wall -D_GNU_SOURCE -disable-dimensions \
    main.c -o bubble -lm
echo "[ok] bubble"

mkjob () {  # $1 dir $2 nom $3 nt $4 wall $5 dep(vide|__DEP__) $6...=lignes
  local dir=$1 name=$2 nt=$3 wall=$4 dep=$5; shift 5
  { echo "#!/bin/bash"; echo "#SBATCH --job-name=$name"; echo "#SBATCH --ntasks=$nt"
    echo "#SBATCH --partition=$PART"; echo "#SBATCH --time=$wall"
    echo "#SBATCH --output=slurm-%j.out"
    [ -n "$dep" ] && echo "#SBATCH --dependency=afterok:$dep"
    echo "cd \"$dir\""; printf '%s\n' "$@"; } > "$dir/job.slurm"
}

mkdir -p "$BASE"
echo "== preparation (9 points x [1 pre + $NMEMBERS bulles]) =="
for i in "${!LABS[@]}"; do
  lab=${LABS[$i]}; ke=${KES[$i]}; b=${BETA[$i]}
  PRE=$BASE/${lab}_pre
  mkdir -p "$PRE"; cp -f "$HERE/bubble" "$PRE/"; cp -n "$ROOT/end_we1" "$PRE/dump"
  #                                       maxlvl MAXTIME R0 FORCED KE INJ BOND OMECO GRAV
  mkjob "$PRE" "${lab}pre" "$NT_PRE" "06:00:00" "" \
    "srun ./bubble 7 $TPRE 8.0 1 $ke 0 1.0 0.20 1.0"
  for m in $(seq 0 $((NMEMBERS-1))); do
    BUB=$BASE/${lab}_m${m}_bub
    mkdir -p "$BUB"; cp -f "$HERE/bubble" "$BUB/"
    #                                       maxlvl MAXTIME R0 FORCED KE INJ BOND OMECO GRAV NU MEMBER
    mkjob "$BUB" "${lab}m${m}" "$NT_BUB" "48:00:00" "__DEP__" \
      "cp ../${lab}_pre/end restart" \
      "srun ./bubble 7 $TBUB 8.0 1 $ke 1 1.0 0.20 $GRAV0 -1 $m"
  done
  echo "[prep] $lab  β=$b  KE=$ke  (1 pre + $NMEMBERS bulles)"
done

echo
NPRE=${#LABS[@]}; NBUB=$((NPRE*NMEMBERS))
echo "TOTAL : $NPRE precurseurs (${NT_PRE}c) + $NBUB bulles (${NT_BUB}c) = $((NPRE+NBUB)) jobs, $PART."
if [ "${LAUNCH:-0}" = 1 ]; then
  for i in "${!LABS[@]}"; do
    lab=${LABS[$i]}
    pid=$(cd "$BASE/${lab}_pre" && sbatch --parsable job.slurm)
    for m in $(seq 0 $((NMEMBERS-1))); do
      sed -i "s/__DEP__/$pid/" "$BASE/${lab}_m${m}_bub/job.slurm"
      ( cd "$BASE/${lab}_m${m}_bub" && sbatch --parsable job.slurm >/dev/null )
    done
    echo "[submit] $lab : pre=$pid + $NMEMBERS bulles (afterok $pid)"
  done
  echo; squeue -u "$USER" -h | wc -l | awk '{print $1" jobs dans la file"}'
else
  echo "DRY-RUN : rien soumis. LAUNCH=1 bash deploy_sweep.sh pour soumettre."
fi
