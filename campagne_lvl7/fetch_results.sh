#!/bin/bash
# ============================================================================
#  fetch_results.sh -- Rapatrie les .dat de la campagne (a executer DANS WSL)
#  Leger : fichiers de diagnostic seulement, jamais les snapshots/dumps.
# ============================================================================
set -euo pipefail
LOCAL="/mnt/c/Users/carlb/OneDrive/Desktop/code stage IMFT"
HOST=${HOST:-kairoslogin}
ROOT=/work/p0910/bonjeanf/bulles

for run in lvl7 lvl7_we2 lvl7_we3 sens_ome008 sens_ome01 sens_ome02 sens_ome04; do
  mkdir -p "$LOCAL/simulations/$run"
  for f in bubble.dat stats.dat cells.dat levels.dat perfs; do
    if scp -q "$HOST:$ROOT/$run/$f" "$LOCAL/simulations/$run/" 2>/dev/null; then
      echo "  [ok] $run/$f"
    else
      echo "  [--] $run/$f absent"
    fi
  done
done

# --- Campagne Design A (balayage intensite/We_t) : bulles + precurseurs ---
for lab in we03 we07 we13 we16 we20; do
  # run bulle : weXX_bub -> simulations/designA/weXX
  mkdir -p "$LOCAL/simulations/designA/$lab"
  for f in bubble.dat stats.dat cells.dat levels.dat; do
    scp -q "$HOST:$ROOT/designA/${lab}_bub/$f" "$LOCAL/simulations/designA/$lab/" 2>/dev/null \
      && echo "  [ok] designA/$lab/$f" || echo "  [--] designA/$lab/$f absent"
  done
  # precurseur : weXX_pre/stats.dat -> simulations/designA/weXX_pre
  mkdir -p "$LOCAL/simulations/designA/${lab}_pre"
  scp -q "$HOST:$ROOT/designA/${lab}_pre/stats.dat" "$LOCAL/simulations/designA/${lab}_pre/" 2>/dev/null \
    && echo "  [ok] designA/${lab}_pre/stats.dat" || echo "  [--] designA/${lab}_pre/stats.dat absent"
done

echo
echo "Post-traitement (Windows ou WSL, depuis la racine du repo) :"
echo "  python scripts/analyse_designA.py --root .   # <-- Design A (fig + tableaux sec9)"
echo "  python scripts/postprocess.py --root . --levels 7 8"
echo "  python scripts/plot_levels.py simulations/sens_ome008 simulations/sens_ome01 \\"
echo "      simulations/sens_ome02 simulations/sens_ome04 --out scripts/figures/amr_sensibilite.png"
echo "  python scripts/plot_levels.py simulations/lvl7 --out scripts/figures/levels_lvl7.png"
