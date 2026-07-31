#!/usr/bin/env python3
"""
plot_levels.py -- Visualisation de la répartition et de l'évolution des niveaux de raffinement AMR.
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load_levels(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                vals = [float(p) for p in line.split()]
            except ValueError:
                continue
            rows.append(vals)
    if not rows:
        return None
    return np.array(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+", help="dossiers de run contenant levels.dat")
    ap.add_argument("--out", default="scripts/figures/amr_levels.png")
    args = ap.parse_args()

    runs = {}
    for d in args.dirs:
        p = os.path.join(d, "levels.dat")
        arr = load_levels(p)
        if arr is None:
            continue
        runs[os.path.basename(os.path.normpath(d))] = arr
    if not runs:
        raise SystemExit("aucun levels.dat lisible")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))

    for name, arr in runs.items():
        ax[0].plot(arr[:, 0], arr[:, 1:].sum(axis=1), label=name)
    ax[0].set_xlabel("temps t")
    ax[0].set_ylabel("cellules (total)")
    ax[0].grid(True, alpha=0.3)
    ax[0].legend(fontsize=8)

    nrun = len(runs)
    width = 0.8 / nrun
    for k, (name, arr) in enumerate(runs.items()):
        counts = arr[-1, 1:]
        lvls = np.nonzero(counts)[0]
        lo, hi = (lvls.min(), lvls.max()) if len(lvls) else (0, 1)
        levels = np.arange(lo, hi + 1)
        ax[1].bar(levels + (k - (nrun - 1) / 2) * width, counts[lo:hi + 1],
                  width=width, label=f"{name} (t={arr[-1, 0]:.0f})")
    ax[1].set_xlabel("niveau de raffinement")
    ax[1].set_ylabel("cellules")
    ax[1].set_yscale("log")
    ax[1].grid(True, alpha=0.3, axis="y")
    ax[1].legend(fontsize=8)

    plt.tight_layout()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    plt.savefig(args.out, dpi=150)
    print(f"Figure AMR générée : {args.out}")

if __name__ == "__main__":
    main()
