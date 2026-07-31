#!/usr/bin/env python3
"""
plot_trajectoires_cote_a_cote.py -- Trajectoires de la bulle aux trois maillages.

Comparaison des trajectoires de la bulle pour les niveaux de raffinement 6, 7 et 8.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
from pathlib import Path

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'figure.titlesize': 12
})

sys.path.insert(0, str(Path(__file__).parent))
from dataio import load_table

OUT = Path("scripts/figures/trajectoires_cote_a_cote.png")
L0 = 120.0
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a84"
COL = {6: "#2a78d6", 7: "#eb6834", 8: "#008300"}
TMAX = 209.0


def unwrap(c):
    return np.unwrap(c * 2 * np.pi / L0) * L0 / (2 * np.pi)


d = {}
for lv in (6, 7, 8):
    bu = load_table(f"simulations/convturb/lvl{lv}/bubble.dat", verbose=False)
    fr = load_table(f"simulations/convturb/lvl{lv}/frame.dat", verbose=False)
    m = bu[:, 0] <= TMAX
    t = bu[m, 0]
    x = unwrap(bu[m, 2]); x -= x[0]          # positions LAB horizontales
    y = unwrap(bu[m, 3]); y -= y[0]
    zc = unwrap(bu[m, 4])
    z = np.interp(t, fr[:, 0], fr[:, 2]) + (zc - zc[0]); z -= z[0]
    d[lv] = dict(t=t, x=x, y=y, z=z)

fig = plt.figure(figsize=(13.5, 8.5))
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.36, wspace=0.32)

allx = np.concatenate([d[l]["x"] for l in (6, 7, 8)])
ally = np.concatenate([d[l]["y"] for l in (6, 7, 8)])
lim = 1.1 * max(np.abs(allx).max(), np.abs(ally).max())


def style(a):
    a.grid(alpha=.25, lw=.6, color=INK3); a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_color(INK3)
    a.tick_params(colors=INK2, labelsize=11)


# HAUT : errance horizontale par niveau
for i, lv in enumerate((6, 7, 8)):
    a = fig.add_subplot(gs[0, i]); style(a)
    a.plot(d[lv]["x"], d[lv]["y"], lw=1.4, color=COL[lv], alpha=.85)
    a.plot(0, 0, "o", ms=8, color="k", mec="white", mew=1.5, zorder=5)
    a.plot(d[lv]["x"][-1], d[lv]["y"][-1], "s", ms=8, color=COL[lv],
           mec="white", mew=1.5, zorder=5)
    a.set_aspect("equal")
    a.set_xlim(-lim, lim); a.set_ylim(-lim, lim)
    a.axhline(0, color=INK3, lw=.6, ls=":"); a.axvline(0, color=INK3, lw=.6, ls=":")
    a.set_title(f"Niveau {lv}", color=COL[lv], pad=6, fontweight="bold")
    a.set_xlabel(r"$x - x_0$  [u.l.]", color=INK)
    if i == 0:
        a.set_ylabel(r"$y - y_0$  [u.l.]", color=INK)

# BAS GAUCHE : les trois superposées
a = fig.add_subplot(gs[1, 0]); style(a)
for lv in (6, 7, 8):
    a.plot(d[lv]["x"], d[lv]["y"], lw=1.6, color=COL[lv], label=f"Niveau {lv}", alpha=.9)
    a.plot(d[lv]["x"][-1], d[lv]["y"][-1], "s", ms=7, color=COL[lv],
           mec="white", mew=1.2, zorder=5)
a.plot(0, 0, "o", ms=9, color="k", mec="white", mew=1.5, zorder=6)
a.annotate("Injection $(x_0, y_0)$", xy=(0, 0), xytext=(lim * .1, -lim * .55),
           color=INK2, arrowprops=dict(arrowstyle="->", color=INK3, lw=1))
a.set_aspect("equal"); a.set_xlim(-lim, lim); a.set_ylim(-lim, lim)
a.set_xlabel(r"$x - x_0$  [u.l.]", color=INK)
a.set_ylabel(r"$y - y_0$  [u.l.]", color=INK)
a.set_title("Trajectoires horizontales superposées", color=INK, pad=6)
a.legend(frameon=False, labelcolor=INK2, loc="upper left")

# BAS MILIEU : la montée z(t)
a = fig.add_subplot(gs[1, 1]); style(a)
for lv in (6, 7, 8):
    a.plot(d[lv]["t"] - d[lv]["t"][0], d[lv]["z"], lw=1.8, color=COL[lv],
           label=f"Niveau {lv}")
a.set_xlabel("t - t$_{inj}$  [u.t.]", color=INK)
a.set_ylabel(r"$z - z_0$  [u.l.]", color=INK)
a.set_title(r"Position verticale $z - z_0$", color=INK, pad=6)
a.legend(frameon=False, labelcolor=INK2, loc="upper left")

# BAS DROITE : écart horizontal entre paires
a = fig.add_subplot(gs[1, 2]); style(a)
tc = np.linspace(0, TMAX - 161, 400)


def sep(la, lb):
    xa = np.interp(tc, d[la]["t"] - d[la]["t"][0], d[la]["x"])
    ya = np.interp(tc, d[la]["t"] - d[la]["t"][0], d[la]["y"])
    xb = np.interp(tc, d[lb]["t"] - d[lb]["t"][0], d[lb]["x"])
    yb = np.interp(tc, d[lb]["t"] - d[lb]["t"][0], d[lb]["y"])
    return np.hypot(xa - xb, ya - yb)


a.plot(tc, sep(6, 7), lw=1.8, color="#4a3aa7", label="Niveau 6 vs 7")
a.plot(tc, sep(7, 8), lw=1.8, color="#e34948", label="Niveau 7 vs 8")
a.axhline(16, color=INK3, ls=":", lw=1.2)
a.text(1, 17, "Diamètre $D = 16$", color=INK2)
a.set_xlabel("t - t$_{inj}$  [u.t.]", color=INK)
a.set_ylabel("Écart horizontal  [u.l.]", color=INK)
a.set_title("Écart inter-trajectoires avec saturation", color=INK, pad=6)
a.legend(frameon=False, labelcolor=INK2, loc="upper left")

fig.suptitle("Trajectoires de la bulle selon le niveau de maillage ($We_t \\approx 1$)",
             color=INK, y=0.99)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=150, facecolor="#fcfcfb", bbox_inches="tight")
print(f"Position initiale commune (0,0) ; carre = position finale ({TMAX:.0f})")
for lv in (6, 7, 8):
    print(f"  lvl{lv} : position finale (x,y)=({d[lv]['x'][-1]:+.1f},{d[lv]['y'][-1]:+.1f})"
          f"  z={d[lv]['z'][-1]:.0f}")
print(f"\n-> {OUT}")
