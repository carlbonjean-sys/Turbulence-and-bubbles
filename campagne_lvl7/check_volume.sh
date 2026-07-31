#!/bin/bash
# ============================================================================
#  check_volume.sh -- surveille la NON-FRAGMENTATION des bulles Design A
#  (a executer DANS WSL : wsl bash campagne_lvl7/check_volume.sh)
#  V_theo = 4/3 pi R0^3 = 2144.66 (R0=8). Alerte si volume < 85 % (fragmentation)
#  ou si chi (col 13) explose. Lit la derniere ligne de chaque bubble.dat.
# ============================================================================
set -uo pipefail
HOST=${HOST:-kairoslogin}
BASE=/work/p0910/bonjeanf/bulles/designA
ssh "$HOST" 'bash -s' <<'EOF'
BASE=/work/p0910/bonjeanf/bulles/designA
VTHEO=2144.66
printf "%-8s %8s %10s %7s %7s  %s\n" run t_dernier vol %V_theo chi etat
for d in "$BASE"/*_bub; do
  [ -f "$d/bubble.dat" ] || { printf "%-8s %8s\n" "$(basename "$d")" "(pas encore)"; continue; }
  last=$(grep -v '^#' "$d/bubble.dat" | tail -1)
  [ -z "$last" ] && continue
  t=$(echo "$last"  | awk '{print $1}')
  v=$(echo "$last"  | awk '{print $2}')
  chi=$(echo "$last"| awk '{print $13}')
  pct=$(awk -v v="$v" -v vt="$VTHEO" 'BEGIN{printf "%.1f", 100*v/vt}')
  etat=$(awk -v p="$pct" -v c="$chi" 'BEGIN{
           if (p<85) print "!! FRAGMENTE ?"; else if (c>3) print "! tres deforme";
           else print "ok intacte"}')
  printf "%-8s %8.1f %10.1f %6s%%  %6.2f  %s\n" "$(basename "$d" _bub)" "$t" "$v" "$pct" "$chi" "$etat"
done
EOF
