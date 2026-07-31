#!/usr/bin/env python3
"""
dataio.py -- Chargement et lecture des fichiers .dat Basilisk.
"""
import os
import sys
from collections import Counter
import numpy as np

def load_table(path, verbose=False):
    """Charge un fichier .dat et renvoie un tableau numpy."""
    if not os.path.isfile(path):
        return None

    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        try:
            rows.append([float(x) for x in s.split()])
        except ValueError:
            continue
    if not rows:
        return None

    widths = Counter(len(r) for r in rows)
    formats = sorted(w for w, c in widths.items() if c >= 2)
    if not formats:
        formats = sorted(widths)

    ncommon = min(formats)
    kept = [r[:ncommon] for r in rows if len(r) >= ncommon]

    if verbose and len(formats) > 1:
        print(f"[dataio] {os.path.basename(path)}: multiple formats detected, truncating to {ncommon} cols", file=sys.stderr)

    return np.array(kept) if kept else None

if __name__ == "__main__":
    for p in sys.argv[1:]:
        a = load_table(p)
        print(f"{p} -> {None if a is None else a.shape}")
