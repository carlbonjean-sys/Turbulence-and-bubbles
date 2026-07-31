#!/usr/bin/env python3
"""
plot_campagne_beta.py -- Cartographie des points d'ensemble acquis dans l'espace (beta, Fr').
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

OUT = Path("scripts/figures/campagne_beta_froude.png")
U_LAM, G, D, T0 = 12.27, 4.0, 16.0, 220.0
SQ_GD = np.sqrt(G * D)

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
fr = beta * U_LAM / SQ_GD
ratio = np.array([r["mean"] for r in rows]) / U_LAM

fig, ax = plt.subplots(figsize=(7.5, 5))
ax.grid(alpha=.25, lw=.6, color="#8a8a84")
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

sc = ax.scatter(beta, fr, c=ratio, s=80, cmap="viridis", edgecolors="black", linewidths=0.8, zorder=3)
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label(r"Vitesse d'ascension relative $\bar{u}_{\infty,\text{turb}} / u_\infty$")

for tag, b, f_val in zip([r["tag"] for r in rows], beta, fr):
    ax.annotate(tag, (b, f_val), textcoords="offset points", xytext=(5, 5), fontsize=9)

ax.set_xlabel(r"Intensité turbulente relative $\beta = u'/u_\infty$")
ax.set_ylabel(r"Nombre de Froude turbulent $Fr' = u'/\sqrt{gD}$")

plt.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT, dpi=300)
plt.close()
print("Figure de campagne générée dans scripts/figures/")
