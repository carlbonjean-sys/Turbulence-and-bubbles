#!/usr/bin/env python3
"""
plot_campagne_beta.py -- Plan de couverture de la campagne d'intensité turbulente.

Présentation des points d'intensité turbulente acquis et des extensions ciblées.
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
OUT1 = Path("scripts/figures/couverture_campagne.png")
OUT2 = Path("scripts/figures/campagne_beta.png")
U_LAM, G, D, T0 = 12.27, 4.0, 16.0, 220.0
SQ_GD = np.sqrt(G * D)
SIGMA = (1 - 1 / 850) * G * (2 * 8.0) ** 2 / 1.0     # 1022.8

C_OLD, C_NEW, C_LIT = "#2a78d6", "#eb6834", "#008300"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a84"

EXIST = [
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
NOUV = [("hi085", 0.85), ("hi100", 1.00)]

rows = []
for tag, b, ppat in EXIST:
    v = []
    for m in range(5):
        f = Path(str(ppat).format(m=m))
        a = load_table(f, verbose=False)
        if a is None:
            continue
        s = a[:, 0] >= T0
        if s.sum() > 2:
            v.append(a[s, 1].mean())
    v = np.array(v)
    rows.append(dict(tag=tag, beta=b, n=len(v), mean=v.mean(),
                     sem=v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan))

beta = np.array([r["beta"] for r in rows])
ratio = np.array([r["mean"] for r in rows]) / U_LAM
err = np.array([r["sem"] for r in rows]) / U_LAM
w = 1 / err ** 2
c = np.sum(w * (1 - ratio) * beta ** 2) / np.sum(w * beta ** 4)

print(f"{'point':7} {'beta':>5} {'Fr':>6} {'We_t':>6} {'n':>2}  {'ratio':>6}")
for r in rows:
    ke = 1.5 * (r["beta"] * U_LAM) ** 2
    print(f"{r['tag']:7} {r['beta']:5.2f} {r['beta']*U_LAM/SQ_GD:6.3f} "
          f"{(2*ke/3)*D/SIGMA:6.2f} {r['n']:2d}  {r['mean']/U_LAM:6.3f}")
print("  --- à venir ---")
for tag, b in NOUV:
    ke = 1.5 * (b * U_LAM) ** 2
    print(f"{tag:7} {b:5.3f} {b*U_LAM/SQ_GD:6.3f} {(2*ke/3)*D/SIGMA:6.2f}  5")

fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2),
                       gridspec_kw=dict(width_ratios=[1.45, 1]))
for a in ax:
    a.grid(alpha=.25, lw=.6, color=INK3); a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_color(INK3)
    a.tick_params(colors=INK2, labelsize=11)

# Panneau A : Couverture en beta
a = ax[0]
bb = np.linspace(0.005, 1.05, 400)
a.plot(bb, 1 - c * bb ** 2, color=INK3, lw=1.5, ls="--", zorder=1,
       label=rf"Loi quadratique $1 - {c:.2f}\beta^2$")
frb = bb * U_LAM / SQ_GD
m = frb > 0.42
a.plot(bb[m], 0.37 / frb[m], color=C_LIT, lw=2.0, zorder=2,
       label=r"Régime asymptotique $0{,}37/Fr'$ (Liu & Deike 2024)")
a.errorbar(beta, ratio, yerr=err, fmt="o", ms=8, color=C_OLD, lw=0,
           elinewidth=1.6, capsize=4, ecolor=C_OLD, mec="white", mew=2, zorder=4,
           label=r"Points acquis ($\beta \in [0{,}15 \,;\, 0{,}38]$)")
a.plot([0], [1], "s", ms=8, color=C_OLD, mec="white", mew=2, zorder=4)

for tag, b in NOUV:
    fr = b * U_LAM / SQ_GD
    yy = 0.37 / fr if fr > 0.42 else 1 - c * b ** 2
    a.plot(b, yy, "o", ms=8, mfc="none", mec=C_NEW, mew=2.0, zorder=3)
a.plot([], [], "o", ms=8, mfc="none", mec=C_NEW, mew=2.0,
       label=r"Extension ciblée ($\beta \in [0{,}05 \,;\, 1{,}00]$)")

a.axvspan(0.04, 0.11, color=C_NEW, alpha=.07, lw=0)
a.axvspan(0.46, 1.04, color=C_NEW, alpha=.07, lw=0)
a.text(0.075, 1.04, "Faible $\\beta$\n(validation quadratique)", color=C_NEW,
       ha="center")
a.text(0.75, 1.04, "Fort $\\beta$\n(régime $0{,}37/Fr'$)",
       color=C_NEW, ha="center")
a.axhline(1, color=INK3, lw=.8, ls=":", zorder=0)
a.set_xlabel(r"Intensité turbulente $\beta = u'/u_\infty$", color=INK)
a.set_ylabel(r"Vitesse relative $\bar{u}_\infty^{\,turb}/u_\infty$", color=INK)
a.set_title(r"Plage de l'intensité turbulente $\beta$", color=INK, pad=9)
a.legend(frameon=False, loc="lower left", labelcolor=INK2)
a.set_xlim(-0.03, 1.06); a.set_ylim(0.15, 1.12)

# Panneau B : We_t par point
a = ax[1]
allb = np.array(sorted([r["beta"] for r in rows] + [b for _, b in NOUV]))
ke_all = 1.5 * (allb * U_LAM) ** 2
we_all = (2 * ke_all / 3) * D / SIGMA
isnew = np.array([b not in [r["beta"] for r in rows] for b in allb])
a.bar(np.arange(len(allb))[~isnew], we_all[~isnew], color=C_OLD, width=.62, zorder=2,
      label="Acquis")
a.bar(np.arange(len(allb))[isnew], we_all[isnew], color=C_NEW, width=.62, zorder=2,
      alpha=.85, label="Extension")
a.axhline(3.0, color="#e34948", ls="--", lw=1.6, zorder=3,
          label=r"Seuil de rupture $We_c = 3$")
a.set_xticks(np.arange(len(allb)))
a.set_xticklabels([f"{b:.2f}" for b in allb], rotation=45)
a.set_xlabel(r"Intensité turbulente $\beta$", color=INK)
a.set_ylabel(r"Weber turbulent $We_t = \rho u'^2 D / \sigma$", color=INK)
a.set_title(r"Nombre de Weber turbulent $We_t$", color=INK, pad=9)
a.legend(frameon=False, loc="upper left", labelcolor=INK2)
a.set_ylim(0, 3.6)
a.annotate(f"$We_{{t,max}} = {we_all.max():.2f}$", xy=(len(allb) - 1, we_all.max()),
           xytext=(len(allb) - 4.5, 2.75), color=INK,
           arrowprops=dict(arrowstyle="->", color=INK3, lw=1.0))

fig.suptitle(r"Couverture de la campagne numérique en fonction de $\beta$",
             color=INK, y=0.99)
fig.tight_layout()
OUT1.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT1, dpi=150, facecolor="#fcfcfb")
fig.savefig(OUT2, dpi=150, facecolor="#fcfcfb")
print(f"\n-> {OUT1}")
print(f"-> {OUT2}")
