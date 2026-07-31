#!/usr/bin/env python3
"""
plot_courbe_maitresse.py -- Effet de la turbulence sur la vitesse d'ascension.

Vitesse d'ascension relative u_turb / u_inf en fonction de l'intensité
turbulente beta = u'/u_inf et du nombre de Froude turbulent Fr' = u'/sqrt(gD).
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

RUNS = Path("simulations/weber_ensemble")
OUT1 = Path("scripts/figures/ralentissement_beta_froude.png")
OUT2 = Path("scripts/figures/courbe_maitresse.png")
U_LAM, G, D, T0 = 12.27, 4.0, 16.0, 220.0
SQ_GD = np.sqrt(G * D)

C_DATA, C_LIT1, C_LIT2 = "#2a78d6", "#eb6834", "#9b59b6"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a84"

POINTS = [
    ("lo050", 0.050, Path("simulations/betasweep/lo050_m{m}_bub/frame.dat")),
    ("lo075", 0.075, Path("simulations/betasweep/lo075_m{m}_bub/frame.dat")),
    ("lo100", 0.100, Path("simulations/betasweep/lo100_m{m}_bub/frame.dat")),
    ("wt05",  0.150, Path("simulations/weber_ensemble/wt05_m{m}/frame.dat")),
    ("wt11",  0.220, Path("simulations/weber_ensemble/wt11_m{m}/frame.dat")),
    ("wt21",  0.310, Path("simulations/weber_ensemble/wt21_m{m}/frame.dat")),
    ("wt25",  0.330, Path("simulations/weber_ensemble/wt25_m{m}/frame.dat")),
    ("wt32",  0.380, Path("simulations/weber_ensemble/wt32_m{m}/frame.dat")),
    ("hi050", 0.500, Path("simulations/betasweep/hi050_m{m}_bub/frame.dat")),
    ("hi065", 0.650, Path("simulations/betasweep/hi065_m{m}_bub/frame.dat")),
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

# Ajustement quadratique pondéré
w = 1.0 / err ** 2
c = np.sum(w * (1 - ratio) * beta ** 2) / np.sum(w * beta ** 4)

print(f"{'pt':6} {'beta':>5} {'Fr':>6} {'n':>2} {'u_moy':>7} {'ratio':>6} "
      f"{'reduc':>7} {'SEM':>6}   membres")
for r, rr, ee, ff in zip(rows, ratio, err, fr):
    print(f"{r['tag']:6} {r['beta']:5.2f} {ff:6.3f} {r['n']:2d} {r['mean']:7.3f} "
          f"{rr:6.3f} {100*(rr-1):+6.1f}% {100*ee:5.1f}%   "
          + " ".join(f"{x:.2f}" for x in r["vals"]))
print(f"\nLoi ajustée : ratio = 1 - {c:.2f} beta^2")

fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2))
for a in ax:
    a.grid(alpha=.25, lw=.6, color=INK3)
    a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_color(INK3)
    a.tick_params(colors=INK2, labelsize=11)

# Panneau A : vs beta
a = ax[0]
bb = np.linspace(0, 0.70, 200)
a.plot(bb, 1 - c * bb ** 2, color=INK3, lw=1.6, ls="--", zorder=1,
       label=rf"Ajustement $1 - {c:.2f}\,\beta^2$")
for r in rows:
    a.plot(np.full(r["n"], r["beta"]), r["vals"] / U_LAM, "o", ms=4.5,
           color=C_DATA, alpha=.30, mec="none", zorder=2)
a.errorbar(beta, ratio, yerr=err, fmt="o", ms=8, color=C_DATA, lw=0,
           elinewidth=1.6, capsize=4, ecolor=C_DATA,
           mec="white", mew=2, zorder=3, label="Moyenne d'ensemble ($n=5$)")
a.plot([0], [1], "s", ms=8, color=C_DATA, mec="white", mew=2, zorder=3)
a.annotate("Référence laminaire\n$u_\\infty = 12{,}27$", xy=(0, 1), xytext=(0.03, 1.03),
           color=INK2, arrowprops=dict(arrowstyle="-", color=INK3, lw=.9))

a.annotate(f"{100*(ratio[0]-1):+.1f} %", xy=(beta[0], ratio[0]),
           xytext=(beta[0] + 0.015, ratio[0] + 0.012),
           color=INK, fontweight="bold")
a.annotate(f"{100*(ratio[-1]-1):+.1f} %", xy=(beta[-1], ratio[-1]),
           xytext=(beta[-1] - 0.09, ratio[-1] + 0.04),
           color=INK, fontweight="bold")
a.axhline(1, color=INK3, lw=.8, ls=":", zorder=0)
a.set_xlabel(r"Intensité turbulente $\beta = u'/u_\infty$", color=INK)
a.set_ylabel(r"Vitesse relative $\bar{u}_\infty^{\,turb}\,/\,u_\infty$", color=INK)
a.set_title(r"Ralentissement en fonction de $\beta$", color=INK, pad=9)
a.set_xlim(-0.02, 0.70); a.set_ylim(0.20, 1.10)
a.legend(frameon=False, loc="lower left", labelcolor=INK2)

# Panneau B : vs Fr'
a = ax[1]
ff = np.linspace(0.02, 1.05, 300)

# Branche Spelt & Biesheuvel (1997) : 1 - 1.55 * Fr'^2
a.plot(ff, 1 - 1.55 * ff ** 2, color=C_LIT2, lw=1.8, ls="-.", zorder=2,
       label=r"Spelt & Biesheuvel (1997) : $1 - 1{,}55\,Fr'^2$")

# Branche Liu & Deike (2024) : 0.37 / Fr'
m = ff > 0.30
a.plot(ff[m], 0.37 / ff[m], color=C_LIT1, lw=2.0, zorder=2,
       label=r"Liu & Deike (2024) : $0{,}37/Fr'$")

# Ajustement local DNS
cf = c * (SQ_GD / U_LAM) ** 2
a.plot(ff, 1 - cf * ff ** 2, color=INK3, lw=1.6, ls="--", zorder=1,
       label=rf"Ajustement DNS : $1 - {cf:.2f}\,Fr'^2$")

for r, ffr in zip(rows, fr):
    a.plot(np.full(r["n"], ffr), r["vals"] / U_LAM, "o", ms=4.5,
           color=C_DATA, alpha=.30, mec="none", zorder=3)
a.errorbar(fr, ratio, yerr=err, fmt="o", ms=8, color=C_DATA, lw=0,
           elinewidth=1.6, capsize=4, ecolor=C_DATA,
           mec="white", mew=2, zorder=4, label="Données DNS ($n=5$)")
a.axhline(1, color=INK3, lw=.8, ls=":", zorder=0)
a.set_xlabel(r"Froude turbulent $Fr' = u'/\sqrt{gD}$", color=INK)
a.set_ylabel(r"Vitesse relative $\bar{u}_\infty^{\,turb}\,/\,u_\infty$", color=INK)
a.set_title(r"Comparaison aux régimes de la littérature", color=INK, pad=9)
a.set_xlim(0, 1.05); a.set_ylim(0.20, 1.10)
a.legend(frameon=False, loc="lower left", labelcolor=INK2)

fig.suptitle(r"Vitesse d'ascension relative d'une bulle en turbulence isotrope ($Bo = 1$, $Ga = 70$)",
             color=INK, y=0.99)
fig.tight_layout()
OUT1.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT1, dpi=150, facecolor="#fcfcfb")
fig.savefig(OUT2, dpi=150, facecolor="#fcfcfb")
print(f"\n-> {OUT1}")
print(f"-> {OUT2}")
