#!/usr/bin/env python3
"""
plot_weber_turbulent.py -- Sécurité vis-à-vis de la fragmentation.

Calcul et tracé du Weber turbulent We_t = rho * u'^2 * D / sigma en fonction de beta
pour valider l'absence de fragmentation (We_t <= 2.36 < We_c = 3).
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

OUT = Path("scripts/figures/weber_turbulent_points.png")
U_LAM, G, D = 12.27, 4.0, 16.0
SIGMA = (1 - 1 / 850) * G * (2 * 8.0) ** 2 / 1.0     # 1022.8

C_OLD, C_NEW = "#2a78d6", "#eb6834"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a84"

EXIST_BETA = np.array([0.05, 0.075, 0.10, 0.15, 0.22, 0.31, 0.33, 0.38, 0.50, 0.65])
EXT_BETA = np.array([0.85, 1.00])

def calc_we(b):
    u_prime = b * U_LAM
    ke = 1.5 * u_prime ** 2
    return (2.0 * ke / 3.0) * D / SIGMA

we_exist = calc_we(EXIST_BETA)
we_ext = calc_we(EXT_BETA)

fig, ax = plt.subplots(figsize=(8.5, 5.2))

ax.grid(alpha=.25, lw=.6, color=INK3)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(INK3)
ax.tick_params(colors=INK2, labelsize=11)

# Tracé continu We_t(beta)
b_continuous = np.linspace(0, 1.05, 300)
we_continuous = calc_we(b_continuous)
ax.plot(b_continuous, we_continuous, color=INK3, lw=1.5, ls=":", label=r"$We_t(\beta) = \frac{\rho (u'_\infty \beta)^2 D}{\sigma}$")

# Points acquis et extension
ax.plot(EXIST_BETA, we_exist, "o", ms=9, color=C_OLD, mec="white", mew=1.8, label=r"Points acquis ($\beta \leq 0{,}65$)", zorder=4)
ax.plot(EXT_BETA, we_ext, "s", ms=8, color=C_NEW, mec="white", mew=1.8, label=r"Extension ciblée ($\beta \leq 1{,}00$)", zorder=4)

# Seuil critique We_c = 3
ax.axhline(3.0, color="#e34948", ls="--", lw=1.8, zorder=3, label=r"Seuil de rupture critique $We_c = 3{,}0$")
ax.axhline(2.36, color=INK2, ls="-.", lw=1.2, zorder=2)

ax.text(0.02, 3.08, r"Zone de fragmentation ($We_t > We_c$)", color="#e34948", fontweight="bold")
ax.text(0.02, 2.44, r"Valeur maximale de la campagne : $We_{t,max} = 2{,}36$", color=INK)

ax.set_xlabel(r"Intensité turbulente $\beta = u'/u_\infty$", color=INK)
ax.set_ylabel(r"Nombre de Weber turbulent $We_t = \rho u'^2 D / \sigma$", color=INK)
ax.set_title(r"Évolution du nombre de Weber turbulent avec l'intensité turbulente", color=INK, pad=10)
ax.set_xlim(-0.02, 1.05)
ax.set_ylim(0, 3.5)
ax.legend(frameon=False, loc="center left", labelcolor=INK2)

fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=150, facecolor="#fcfcfb")
print(f"-> {OUT}")
