#!/usr/bin/env python3
"""
plot_trajectoires_3d.py -- Trajectoires 3D de la bulle aux 3 maillages.
Produit deux figures :
  - trajectoires_xz.png : le plan x-z (vue de cote), cote a cote + superpose
  - trajectoires_3d.png : le cube en perspective, trajectoire + ombres portees

Runs convturb (we07, KE=11, beta=0.22, meme injection). Fenetre commune [161,209]
(lvl8 partiel). Coordonnees deroulees (domaine periodique), origine ramenee a 0.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dataio import load_table

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
    # x,y : positions LAB horizontales (le repere mobile ne bouge que selon z)
    x = unwrap(bu[m, 2]); x -= x[0]
    y = unwrap(bu[m, 3]); y -= y[0]
    # z_lab : VRAIE montee = deplacement du repere (frame_z) + residuel dans la grille.
    # zc seul ~ constant (le repere suit la bulle), il ne donne PAS la montee.
    zc = unwrap(bu[m, 4])
    frame_z = np.interp(t, fr[:, 0], fr[:, 2])
    z = frame_z + (zc - zc[0]); z -= z[0]
    d[lv] = dict(t=t, x=x, y=y, z=z)


def style2d(a):
    a.grid(alpha=.25, lw=.6, color=INK3); a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_color(INK3)
    a.tick_params(colors=INK2, labelsize=9)


# ============================================================ (1) plan x-z
zmax = max(d[l]["z"].max() for l in (6, 7, 8))
xlo = min(d[l]["x"].min() for l in (6, 7, 8)) - 5
xhi = max(d[l]["x"].max() for l in (6, 7, 8)) + 5

fig, ax = plt.subplots(1, 4, figsize=(14, 5.2),
                       gridspec_kw=dict(width_ratios=[1, 1, 1, 1.15]))
for i, lv in enumerate((6, 7, 8)):
    a = ax[i]; style2d(a)
    a.plot(d[lv]["x"], d[lv]["z"], lw=1.5, color=COL[lv])
    a.plot(0, 0, "o", ms=9, color="k", mec="white", mew=1.5, zorder=5)
    a.plot(d[lv]["x"][-1], d[lv]["z"][-1], "s", ms=9, color=COL[lv],
           mec="white", mew=1.5, zorder=5)
    a.set_xlim(xlo, xhi); a.set_ylim(-15, zmax * 1.05)
    a.set_title(f"niveau {lv}" + ("" if lv < 8 else " (partiel)"),
                fontsize=11, color=COL[lv], pad=6, fontweight="bold")
    a.set_xlabel(r"$x - x_0$", fontsize=9.5, color=INK)
    if i == 0:
        a.set_ylabel(r"montée  $z - z_0$", fontsize=10, color=INK)

a = ax[3]; style2d(a)
for lv in (6, 7, 8):
    a.plot(d[lv]["x"], d[lv]["z"], lw=1.6, color=COL[lv], label=f"niveau {lv}")
    a.plot(d[lv]["x"][-1], d[lv]["z"][-1], "s", ms=7, color=COL[lv],
           mec="white", mew=1.2, zorder=5)
a.plot(0, 0, "o", ms=10, color="k", mec="white", mew=1.5, zorder=6)
a.set_xlim(xlo, xhi); a.set_ylim(-15, zmax * 1.05)
a.set_title("superposées", fontsize=11, color=INK, pad=6)
a.set_xlabel(r"$x - x_0$", fontsize=9.5, color=INK)
a.legend(fontsize=8.5, frameon=False, labelcolor=INK2, loc="lower left")

fig.suptitle("Plan $x$--$z$ (vue de côté) : la bulle monte en serpentant, "
             "chaque maillage suit un chemin différent", fontsize=12, color=INK, y=1.0)
fig.tight_layout()
fig.savefig("scripts/figures/trajectoires_xz.png", dpi=150, facecolor="#fcfcfb")
print("-> scripts/figures/trajectoires_xz.png")

# ============================================================ (2) cube 3D
fig = plt.figure(figsize=(13, 6.5))

xr = (min(d[l]["x"].min() for l in (6, 7, 8)) - 5,
      max(d[l]["x"].max() for l in (6, 7, 8)) + 5)
yr = (min(d[l]["y"].min() for l in (6, 7, 8)) - 5,
      max(d[l]["y"].max() for l in (6, 7, 8)) + 5)
zr = (-10, zmax * 1.05)


def draw_cube(ax):
    """trajectoires + ombres portees sur les 3 parois, cube en perspective"""
    for lv in (6, 7, 8):
        x, y, z = d[lv]["x"], d[lv]["y"], d[lv]["z"]
        # ombres portees (projections) sur les parois -- lisent le chemin en 2D
        ax.plot(x, y, zr[0], color=COL[lv], lw=.8, alpha=.25)          # sol x-y
        ax.plot(x, yr[1], z, color=COL[lv], lw=.8, alpha=.25)          # paroi x-z
        ax.plot(xr[0], y, z, color=COL[lv], lw=.8, alpha=.25)          # paroi y-z
        # la trajectoire 3D
        ax.plot(x, y, z, color=COL[lv], lw=2.0, label=f"niveau {lv}")
        ax.scatter([x[-1]], [y[-1]], [z[-1]], color=COL[lv], s=45,
                   edgecolors="white", linewidths=1.2, depthshade=False, zorder=6)
    ax.scatter([0], [0], [0], color="k", s=70, edgecolors="white",
               linewidths=1.5, depthshade=False, zorder=7)
    ax.set_xlim(xr); ax.set_ylim(yr); ax.set_zlim(zr)
    ax.set_xlabel(r"$x-x_0$", fontsize=10, color=INK, labelpad=6)
    ax.set_ylabel(r"$y-y_0$", fontsize=10, color=INK, labelpad=6)
    ax.set_zlabel(r"montée $z-z_0$", fontsize=10, color=INK, labelpad=8)
    ax.tick_params(colors=INK2, labelsize=8)
    try:
        ax.set_box_aspect((1, 1, 2.6))
    except Exception:
        pass
    ax.xaxis.pane.set_alpha(.03); ax.yaxis.pane.set_alpha(.03)
    ax.zaxis.pane.set_alpha(.03)


ax1 = fig.add_subplot(1, 2, 1, projection="3d")
draw_cube(ax1)
ax1.view_init(elev=18, azim=-58)
ax1.set_title("Vue en perspective", fontsize=11, color=INK, pad=2)
ax1.legend(fontsize=9, frameon=False, labelcolor=INK2, loc="upper left")
ax1.text2D(0.02, 0.02, "point noir = injection commune\n"
           "traits pâles = ombres portées sur les parois",
           transform=ax1.transAxes, fontsize=7.5, color=INK2)

ax2 = fig.add_subplot(1, 2, 2, projection="3d")
draw_cube(ax2)
ax2.view_init(elev=6, azim=-90)     # presque de face (plan x-z)
ax2.set_title("Vue quasi de côté (plan $x$--$z$)", fontsize=11, color=INK, pad=2)

fig.suptitle("Trajectoires 3D de la bulle aux trois maillages — "
             "elle monte ($z$) en serpentant ($x,y$), chemins distincts",
             fontsize=12.5, color=INK, y=0.99)
fig.savefig("scripts/figures/trajectoires_3d.png", dpi=150, facecolor="#fcfcfb",
            bbox_inches="tight")
print("-> scripts/figures/trajectoires_3d.png")
