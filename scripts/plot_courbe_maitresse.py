#!/usr/bin/env python3
"""
plot_courbe_maitresse.py -- Vitesse d'ascension relative en fonction de beta et Fr'.
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

RUNS = Path("simulations/betasweep")
OUT1 = Path("scripts/figures/ralentissement_beta_froude.png")
OUT2 = Path("scripts/figures/courbe_maitresse.png")
U_LAM, G, D, T0 = 12.27, 4.0, 16.0, 220.0
SQ_GD = np.sqrt(G * D)

C_DATA, C_LIT1, C_LIT2 = "#2a78d6", "#eb6834", "#9b59b6"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a84"

POINTS = [
    ("beta0050", 0.050, Path("simulations/betasweep/beta0050_m{m}_bub/frame.dat")),
    ("beta0075", 0.075, Path("simulations/betasweep/beta0075_m{m}_bub/frame.dat")),
    ("beta0100", 0.100, Path("simulations/betasweep/beta0100_m{m}_bub/frame.dat")),
    ("beta0150", 0.150, Path("simulations/betasweep/beta0150_m{m}_bub/frame.dat")),
    ("beta0220", 0.220, Path("simulations/betasweep/beta0220_m{m}_bub/frame.dat")),
    ("beta0310", 0.310, Path("simulations/betasweep/beta0310_m{m}_bub/frame.dat")),
    ("beta0330", 0.330, Path("simulations/betasweep/beta0330_m{m}_bub/frame.dat")),
    ("beta0380", 0.380, Path("simulations/betasweep/beta0380_m{m}_bub/frame.dat")),
    ("beta0500", 0.500, Path("simulations/betasweep/beta0500_m{m}_bub/frame.dat")),
    ("beta0650", 0.650, Path("simulations/betasweep/beta0650_m{m}_bub/frame.dat")),
]

rows = []
for tag, beta, ppat in POINTS:
    vals = []
    for m in range(5):
        f = Path(str(ppat).format(m=m))
        a = load_table(f, verbose=False)
        if a is None:
            continue
        sel = a[:, 0] >= T0
        if sel.sum() > 2:
            vals.append(a[sel, 1].mean())
    v = np.array(vals)
    n = len(v)
    sem = v.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    rows.append(dict(tag=tag, beta=beta, n=n, vals=v, mean=v.mean(), sem=sem))

beta = np.array([r["beta"] for r in rows])
ratio = np.array([r["mean"] for r in rows]) / U_LAM
err = np.array([r["sem"] for r in rows]) / U_LAM
fr = beta * U_LAM / SQ_GD

w = 1.0 / err ** 2
c = np.sum(w * (1 - ratio) * beta ** 2) / np.sum(w * beta ** 4)

fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2))
for a in ax:
    a.grid(alpha=.25, lw=.6, color=INK3)
    a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_color(INK3)

a = ax[0]
bb = np.linspace(0, 0.70, 200)
a.plot(bb, 1 - c * bb ** 2, color=INK3, lw=1.6, ls="--", zorder=1,
       label=rf"Ajustement $1 - {c:.2f}\,\beta^2$")
for r in rows:
    a.plot(np.full(r["n"], r["beta"]), r["vals"] / U_LAM, "o", ms=4.5,
           color=C_DATA, alpha=.30, mec="none", zorder=2)
a.errorbar(beta, ratio, yerr=err, fmt="o", ms=7, color=C_DATA, ecolor=C_DATA,
           elinewidth=1.6, capsize=3.5, label="DNS (moyenne d'ensemble)", zorder=3)
a.set_xlabel(r"Intensité turbulente $\beta = u'/u_\infty$")
a.set_ylabel(r"Vitesse d'ascension relative $\bar{u}_{\infty,\text{turb}} / u_\infty$")
a.set_ylim(0.55, 1.05)
a.legend(frameon=True, facecolor="white", edgecolor="none")

a = ax[1]
ff = bb * U_LAM / SQ_GD
a.plot(ff, 1 - c * bb ** 2, color=INK3, lw=1.6, ls="--", zorder=1)
a.errorbar(fr, ratio, yerr=err, fmt="o", ms=7, color=C_DATA, ecolor=C_DATA,
           elinewidth=1.6, capsize=3.5, zorder=3)
a.set_xlabel(r"Nombre de Froude turbulent $Fr' = u'/\sqrt{gD}$")
a.set_ylabel(r"$\bar{u}_{\infty,\text{turb}} / u_\infty$")
a.set_ylim(0.55, 1.05)

plt.tight_layout()
OUT1.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT1, dpi=300)
plt.close()

fig, ax = plt.subplots(figsize=(7, 5))
ax.grid(alpha=.25, lw=.6, color=INK3)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(INK3)

ax.errorbar(beta, ratio, yerr=err, fmt="o", ms=7.5, color=C_DATA, ecolor=C_DATA,
            elinewidth=1.6, capsize=3.5, label="DNS présent travail (Bo=1, Ga=70)", zorder=3)
ax.plot(bb, 1 - c * bb ** 2, color=INK, lw=1.8, ls="--",
        label=rf"Loi quadratique $1 - {c:.2f}\,\beta^2$")

ax.set_xlabel(r"Intensité turbulente relative $\beta = u'/u_\infty$")
ax.set_ylabel(r"Vitesse d'ascension relative $\bar{u}_{\infty,\text{turb}} / u_\infty$")
ax.set_ylim(0.55, 1.05)
ax.legend(frameon=True, facecolor="white", edgecolor="none")

plt.tight_layout()
plt.savefig(OUT2, dpi=300)
plt.close()
print("Figure maître générée dans scripts/figures/")
