#!/bin/bash
# ============================================================================
#  deploy_beta_full.sh -- 2026-07-24  (feu vert user explicite)
#
#  Complete la courbe u_inf_turb/u_inf (beta) dans les trois directions :
#
#   A) n=5 PARTOUT   : membres m3,m4 pour wt11/wt21/wt25/wt32 (wt05 les a deja).
#                      Avec n=3, l'incertitude SUR l'ecart-type est de +/-50 %
#                      (1/sqrt(2(n-1))) : on ne connait pas nos barres d'erreur.
#                      A n=5 elle tombe a ~35 %. Precurseurs designA existants,
#                      pas de dependance.
#
#   B) POINTS BAS    : beta = 0.05 / 0.075 / 0.10. Aujourd'hui la parabole
#                      1-3.6 beta^2 est ajustee sur des points tous >= 0.15,
#                      ancres seulement par le point laminaire : on SUPPOSE le
#                      comportement quadratique a l'origine sans le verifier.
#
#   C) POINTS HAUTS  : beta = 0.50 / 0.65 / 0.85 / 1.00, soit ~3x la plage
#                      actuelle (0.38). Le raisonnement "KE=200 fragmente donc
#                      on est bloque" ne s'applique PLUS : ce run etait a GRAV=1
#                      (sigma=255.7), on est maintenant a GRAV=4 (sigma=1022.8),
#                      donc We_t = rho u'^2 D / sigma est 4x plus petit a
#                      turbulence egale. Aux cibles ci-dessous We_t <= 2.36,
#                      sous le We_c=3 de Riviere et al. (2021). Le point haut
#                      approche donc le seuil sans le franchir : s'il fragmente
#                      quand meme, ce n'est pas un echec, c'est une MESURE de
#                      We_c -- resultat dans les deux cas.
#
#  KE_TARGET = 1.5*(beta*u_inf)^2 avec u_inf = 12.27 (laminaire, valide en
#  maillage lvl7 vs lvl8 a 0.32 %).
#
#  ---------------------------------------------------------------------------
#  GARDE-FOU D'EQUILIBRAGE (nouveau -- corrige un piege documente de
#  deploy_sweep.sh : "afterok ne verifie que l'exit code")
#  ---------------------------------------------------------------------------
#  Chaque precurseur se CONTROLE lui-meme avant de rendre la main :
#    - ke atteint dans [0.70, 1.40] x cible   (sinon le point n'est pas ou on croit)
#    - k_max.eta >= 1.5                        (sinon la turbulence est sous-resolue)
#  S'il echoue, il sort en erreur => afterok BLOQUE ses 5 bulles. Sans ce garde-fou
#  un precurseur mal equilibre lance 5 runs de 48 h sur un point fantome.
#
#  Risques connus, annonces :
#   * les precurseurs designA sous-atteignent leur cible de 8 a 15 % (mesure :
#     we03 4.62/5, we20 27.2/32) -- systematique, sans gravite car beta est
#     MESURE a posteriori sur u' reel, jamais suppose. Attendre le meme biais ici.
#   * points BAS : il faut DECROITRE de ke=15 (end_we1) jusqu'a 0.56, soit ~27x,
#     en 84 u.t. Le temps de retournement grandit quand ke tombe => c'est le point
#     le plus a risque de ne pas converger. Le garde-fou le dira (les precurseurs
#     sont rapides). Si echec : decroissance en deux etages depuis we03_pre/end.
#   * point beta=1.00 : k_max.eta extrapole a ~1.9 (seuil 1.5), marge x1.25 --
#     le plus serre de la campagne. Le garde-fou le verifie sur la VRAIE valeur.
#
#  Usage :  bash deploy_beta_full.sh            # DRY-RUN (compile + prepare)
#           LAUNCH=1 bash deploy_beta_full.sh   # soumet
#           PHASES="A" LAUNCH=1 bash ...        # une phase seulement (A, B, C)
# ============================================================================
set -euo pipefail

ROOT=${ROOT:-/work/p0910/bonjeanf/bulles}
BASE=$ROOT/betasweep
PREBASE=$ROOT/designA
ENSBASE=$ROOT/ensemble
HERE=$(cd "$(dirname "$0")" && pwd)
export BASILISK=${BASILISK:-/users/p0910/bonjeanf/basilisk2/basilisk/src}
export PATH="$BASILISK:$PATH"

PART=${PART:-shared-cpu}
GRAV0=${GRAV0:-4}
TPRE=${TPRE:-160}          # les bulles demarrent a t=160 -> plateau t>=220 (uniforme)
TBUB=${TBUB:-280}
NMEM=${NMEM:-5}
PHASES=${PHASES:-"A B C"}

# --- garde-fous d'entree -------------------------------------------------
[ -f "$HERE/main.c" ] || { echo "ERREUR: main.c absent"; exit 1; }
grep -q 'n_regions V_tag straddle' "$HERE/main.c" || {
  echo "ERREUR: main.c n'est PAS la version corrigee (event bubble du 2026-07-22)."
  echo "        Sans elle les nouveaux runs re-ecriraient des volumes faux."; exit 1; }
grep -q 'event move_frame' "$HERE/main.c" || { echo "ERREUR: repere mobile absent"; exit 1; }
grep -q '#include "reduced.h"' "$HERE/main.c" && { echo "ERREUR: reduced.h de retour !"; exit 1; }
[ -f "$ROOT/end_we1" ] || { echo "ERREUR: $ROOT/end_we1 absent"; exit 1; }

echo "== compilation =="
cd "$HERE"
CC99='mpicc -std=c99' qcc -D_MPI=1 -O3 -Wall -D_GNU_SOURCE -disable-dimensions \
    main.c -o bubble -lm
echo "[ok] bubble (main.c corrige : bubble.dat a 16 colonnes)"

# --- le garde-fou d'equilibrage, injecte dans chaque job precurseur ------
GUARD='tail -1 stats.dat | awk -v tgt=__KE__ '"'"'{
  kmax = 3.14159265*128/120; ratio = $3/tgt; keta = kmax*$5;
  printf "[GUARD] ke=%.4f cible=%.4f (%.0f%%)  eta=%.3f  k_max.eta=%.2f\n", $3, tgt, 100*ratio, $5, keta;
  if (ratio < 0.70) { print "[GUARD] ECHEC : forcage non equilibre (<70% de la cible)"; exit 1 }
  if (ratio > 1.40) { print "[GUARD] ECHEC : energie non redescendue (>140% de la cible)"; exit 1 }
  if (keta  < 1.50) { print "[GUARD] ECHEC : turbulence sous-resolue (k_max.eta < 1.5)"; exit 1 }
  print "[GUARD] OK -> les bulles peuvent partir"
}'"'"''

mkjob () {   # $1 dir  $2 nom  $3 ntasks  $4 wall  $5 dep("" ou __DEP__)  $6... lignes
  local dir=$1 name=$2 nt=$3 wall=$4 dep=$5; shift 5
  { echo "#!/bin/bash"
    echo "#SBATCH --job-name=$name"
    echo "#SBATCH --ntasks=$nt"
    echo "#SBATCH --partition=$PART"
    echo "#SBATCH --time=$wall"
    echo "#SBATCH --output=slurm-%j.out"
    [ -n "$dep" ] && echo "#SBATCH --dependency=afterok:$dep"
    echo "cd \"$dir\""
    printf '%s\n' "$@"
  } > "$dir/job.slurm"
}

SUB_NODEP=""      # dossiers a soumettre sans dependance
SUB_PRE=""        # dossiers precurseurs (label:dir)
NA=0; NB=0; NC=0

# ======================================================================== A
if [[ " $PHASES " == *" A "* ]]; then
  echo
  echo "== PHASE A : n=5 partout (membres m3,m4) =="
  LABS=(wt11 wt21 wt25 wt32)
  PRES=(we07 we13 we16 we20)
  KES=(  11   21   25   32)
  for i in "${!LABS[@]}"; do
    pre=${PRES[$i]}; ke=${KES[$i]}; lab=${LABS[$i]}
    [ -f "$PREBASE/${pre}_pre/end" ] || { echo "ERREUR: ${pre}_pre/end absent"; exit 1; }
    for m in 3 4; do
      D=$ENSBASE/${lab}_g${GRAV0}_m${m}_bub
      # ces dossiers existent (prepares en juillet) mais sont VIDES : les jobs
      # 39165-39172 ont ete annules avant de demarrer. On refuse d'ecraser un
      # dossier qui contiendrait des donnees.
      if [ -f "$D/frame.dat" ]; then
        echo "  [skip] $lab m$m : frame.dat deja present, on ne touche pas"
        continue
      fi
      mkdir -p "$D"
      cp -f "$HERE/bubble" "$D/"                    # binaire CORRIGE
      cp -f "$PREBASE/${pre}_pre/end" "$D/restart"
      rm -f "$D/bubble.dat" "$D/frame.dat" "$D/stats.dat" "$D/momentum.dat" "$D/dump"
      mkjob "$D" "${lab}m${m}" 4 "48:00:00" "" \
        "srun ./bubble 7 $TBUB 8.0 1 $ke 1 1.0 0.20 $GRAV0 -1 $m"
      SUB_NODEP="$SUB_NODEP $D"; NA=$((NA+1))
      echo "  [prep] $lab m$m  (KE=$ke, precurseur $pre)"
    done
  done
fi

# ====================================================================== B+C
prep_point () {    # $1 label  $2 beta  $3 KE  $4 nt_pre  $5 wall_pre  $6 nt_bub
  local lab=$1 b=$2 ke=$3 ntp=$4 wp=$5 ntb=$6
  local PRE=$BASE/${lab}_pre
  mkdir -p "$PRE"; cp -f "$HERE/bubble" "$PRE/"
  cp -f "$ROOT/end_we1" "$PRE/dump"          # re-equilibrage monophasique
  rm -f "$PRE/stats.dat" "$PRE/cells.dat" "$PRE/levels.dat"
  mkjob "$PRE" "${lab}pre" "$ntp" "$wp" "" \
    "srun ./bubble 7 $TPRE 8.0 1 $ke 0 1.0 0.20 1.0" \
    "${GUARD//__KE__/$ke}"
  for m in $(seq 0 $((NMEM-1))); do
    local B=$BASE/${lab}_m${m}_bub
    mkdir -p "$B"; cp -f "$HERE/bubble" "$B/"
    rm -f "$B/bubble.dat" "$B/frame.dat" "$B/stats.dat" "$B/momentum.dat" "$B/dump"
    mkjob "$B" "${lab}m${m}" "$ntb" "48:00:00" "__DEP__" \
      "cp ../${lab}_pre/end restart" \
      "srun ./bubble 7 $TBUB 8.0 1 $ke 1 1.0 0.20 $GRAV0 -1 $m"
  done
  SUB_PRE="$SUB_PRE ${lab}"
  # We_t = rho u'^2 D / sigma  avec u'^2 = 2ke/3, D=16, sigma=1022.8
  awk -v b=$b -v ke=$ke 'BEGIN{ printf "  [prep] %-6s beta=%-5s KE=%-7s We_t=%.2f  (1 pre + '"$NMEM"' bulles)\n", "'"$lab"'", b, ke, (2*ke/3)*16/1022.8 }'
}

if [[ " $PHASES " == *" B "* ]]; then
  echo
  echo "== PHASE B : points BAS (teste la loi quadratique a l'origine) =="
  prep_point lo050 0.050 0.56  4 "06:00:00" 4; NB=$((NB+6))
  prep_point lo075 0.075 1.27  4 "06:00:00" 4; NB=$((NB+6))
  prep_point lo100 0.100 2.26  4 "06:00:00" 4; NB=$((NB+6))
fi

if [[ " $PHASES " == *" C "* ]]; then
  echo
  echo "== PHASE C : points HAUTS (triple la plage ; approche We_c sans le franchir) =="
  prep_point hi050 0.50  56.5   6 "12:00:00" 6; NC=$((NC+6))
  prep_point hi065 0.65  95.4   6 "12:00:00" 6; NC=$((NC+6))
  prep_point hi085 0.85 163.2   6 "12:00:00" 6; NC=$((NC+6))
  prep_point hi100 1.00 225.8   6 "12:00:00" 6; NC=$((NC+6))
fi

# ===================================================================== recap
TOT=$((NA+NB+NC))
echo
echo "======================================================================"
echo "  RECAP :  phase A = $NA jobs   phase B = $NB jobs   phase C = $NC jobs"
echo "           TOTAL = $TOT jobs sur $PART"
echo "  Coeurs simultanes max ~6 jobs x 6c = 36  (garde-fou 64) ; walltime 48 h max (garde-fou 72 h)."
echo "  Les bulles de B et C sont en afterok du precurseur de leur point"
echo "  (+ garde-fou d'equilibrage : un precurseur rate BLOQUE ses 5 bulles)."
echo "======================================================================"

if [ "${LAUNCH:-0}" = 1 ]; then
  echo
  for d in $SUB_NODEP; do ( cd "$d" && sbatch --parsable job.slurm | xargs -I{} echo "[submit] $(basename $d) = {}" ); done
  for lab in $SUB_PRE; do
    pid=$(cd "$BASE/${lab}_pre" && sbatch --parsable job.slurm)
    for m in $(seq 0 $((NMEM-1))); do
      sed -i "s/__DEP__/$pid/" "$BASE/${lab}_m${m}_bub/job.slurm"
      ( cd "$BASE/${lab}_m${m}_bub" && sbatch --parsable job.slurm >/dev/null )
    done
    echo "[submit] $lab : pre=$pid + $NMEM bulles (afterok $pid)"
  done
  echo
  squeue -u "$USER" -h | wc -l | awk '{print "  -> "$1" jobs dans la file"}'
else
  echo
  echo "DRY-RUN : RIEN SOUMIS.  LAUNCH=1 bash deploy_beta_full.sh pour soumettre."
fi
