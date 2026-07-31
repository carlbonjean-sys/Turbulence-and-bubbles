#!/bin/bash
# ============================================================================
#  deploy_mechanism.sh -- 2026-07-26
#  Mesure le MECANISME du ralentissement (src/mechanism.c) sur les snapshots
#  DEJA presents : echantillonnage preferentiel vs trainee non lineaire.
#
#  Pur post-traitement en LECTURE (aucune simu) : 1 coeur, ecrit mechanism.dat
#  dans chaque dossier de run. Serie, comme diag_volume.
#
#  Runs traites (trend en beta + un controle) :
#    lam7_frame   beta=0    CONTROLE : pas de turbulence -> signal DOIT etre ~0
#    wt05         beta=0.15 (-6 %)
#    wt21         beta=0.31 (-37 %)
#    wt32         beta=0.38 (-39 %)
#  Si le signal preferentiel croit avec beta ET est nul en laminaire, il est reel.
#
#  Usage :  bash deploy_mechanism.sh            # DRY-RUN
#           LAUNCH=1 bash deploy_mechanism.sh   # soumet (1 coeur)
#           STEP=2 ... pour 1 snapshot sur 2 (moitie du cout)
# ============================================================================
set -euo pipefail
ROOT=${ROOT:-/work/p0910/bonjeanf/bulles}
D=$ROOT/mechanism
HERE=$(cd "$(dirname "$0")" && pwd)
export BASILISK=${BASILISK:-/users/p0910/bonjeanf/basilisk2/basilisk/src}
export PATH="$BASILISK:$PATH"
STEP=${STEP:-1}
PART=${PART:-shared-cpu}

[ -f "$HERE/mechanism.c" ] || { echo "ERREUR: mechanism.c absent"; exit 1; }
grep -q 'periodic (right)' "$HERE/mechanism.c" || {
  echo "ERREUR: mechanism.c ne declare pas periodic() -- piege connu"; exit 1; }

echo "== compilation serie =="
cd "$HERE"
qcc -O3 -Wall -D_GNU_SOURCE -disable-dimensions mechanism.c -o mechanism -lm
echo "[ok] mechanism"

mkdir -p "$D"; cp -f "$HERE/mechanism" "$D/"

# label -> chemin des runs (ROOT vaut deja .../bulles)
declare -A RUNPATH=(
  [lam_b000]="$ROOT/lam7_frame"
  [wt05_b015]="$ROOT/weber_fix/wt05_g4_bub"
  [wt21_b031]="$ROOT/weber_fix/wt21_g4_bub"
  [wt32_b038]="$ROOT/weber_fix/wt32_g4_bub"
)

FOUND=""
for lab in "${!RUNPATH[@]}"; do
  p="${RUNPATH[$lab]}"
  n=$(ls "$p"/snapshot-* 2>/dev/null | wc -l)
  if [ "$n" -gt 0 ]; then echo "  [ok]   $lab : $n snapshots  ($p)"; FOUND="$FOUND $lab"
  else echo "  [SKIP] $lab : aucun snapshot ($p)"; fi
done
[ -n "$FOUND" ] || { echo "ERREUR: aucun run avec snapshots"; exit 1; }

# script du job : boucle sur runs puis snapshots
JOB="$D/job.slurm"
{
  echo "#!/bin/bash"
  echo "#SBATCH --job-name=mecanism"
  echo "#SBATCH --ntasks=1"
  echo "#SBATCH --partition=$PART"
  echo "#SBATCH --time=06:00:00"
  echo "#SBATCH --output=slurm-%j.out"
  echo "cd \"$D\""
  echo "HDR='# label z_c u_bulle uz_domaine c1 c2 c3 c4 h1 h2 h3 h4 f1 f2 f3 f4 Vb'"
  echo "HDR2='# c=cone_sup h=hemisphere f=coquille ; bandes r/D=[1,1.5][1.5,2][2,2.5][2.5,3] ; STOCKE (repere mobile)'"
  for lab in $FOUND; do
    p="${RUNPATH[$lab]}"
    echo "out=\"$p/mechanism.dat\""
    echo "echo \"\$HDR\"  > \"\$out\"; echo \"\$HDR2\" >> \"\$out\""
    echo "k=0"
    echo "for s in \"$p\"/snapshot-*; do"
    echo "  k=\$((k+1)); [ \$(( (k-1) % $STEP )) -eq 0 ] || continue"
    echo "  t=\${s##*snapshot-}"
    echo "  ./mechanism \"\$s\" \"\$t\" >> \"\$out\" 2>/dev/null || true"
    echo "done"
    echo "echo \"[fait] $lab : \$(grep -vc '^#' \"\$out\") lignes\""
  done
} > "$JOB"

echo
echo "1 job, 1 coeur, 06:00:00. Ecrit mechanism.dat dans chaque dossier de run."
if [ "${LAUNCH:-0}" = 1 ]; then
  ( cd "$D" && sbatch --parsable job.slurm | xargs -I{} echo "[submit] mechanism = {}" )
  squeue -u "$USER" -h | wc -l | awk '{print "  -> "$1" jobs dans la file"}'
else
  echo "DRY-RUN. LAUNCH=1 bash deploy_mechanism.sh pour soumettre."
fi
