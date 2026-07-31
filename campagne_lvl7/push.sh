#!/bin/bash
# ============================================================================
#  push.sh -- Pousse la campagne lvl7 sur Kairos (a executer DANS WSL, au labo)
#  (chemins litteraux volontairement : les variables inline se font manger en
#   traversant Git Bash -> wsl -> ssh)
# ============================================================================
set -euo pipefail
LOCAL="/mnt/c/Users/carlb/OneDrive/Desktop/code stage IMFT"
HOST=${HOST:-kairoslogin}
DEST=/work/p0910/bonjeanf/bulles/campagne_lvl7

ssh "$HOST" "mkdir -p $DEST"
scp "$LOCAL/src/main.c" "$LOCAL/campagne_lvl7/deploy.sh" "$HOST:$DEST/"

echo
echo "[ok] pousse sur $HOST:$DEST"
echo "Suite, sur Kairos :"
echo "  ssh $HOST"
echo "  cd $DEST && bash deploy.sh          # dry-run (compile + prepare)"
echo "  LAUNCH=1 bash deploy.sh             # soumet les 6 jobs"
