#!/usr/bin/env python3
"""
plot_convergence_maillage.py -- Convergence en maillage du cas laminaire.

Étude de sensibilité à la résolution spatiale pour la montée d'une bulle en liquide au repos.
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

OUT = Path("scripts/figures/convergence_maillage.png")
L0, R0, G = 120.0, 8.0, 4.0
D = 2 * R0
MU = 0.01 * (L0 / (2 * np.pi)) ** 2 / 2
V_THEO = 4 / 3 * np.pi * R0 ** 3
T0_PLATEAU = 220.0

RUNS = [
    (6, "simulations/validation/lam6", "#2a78d6"),
    (7, "simulations/lam7_frame",      "#eb6834"),
    (8, "simulations/validation/lam8", "#008300"),
]
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a84"


def smooth(y, n=15):
    return np.convolve(np.pad(y, (n // 2, n - 1 - n // 2), mode="edge"),
                       np.ones(n) / n, mode="valid")


res = []
for lv, d, col in RUNS:
    fr = load_table(f"{d}/frame.dat", verbose=False)
    bu = load_table(f"{d}/bubble.dat", verbose=False)
    t = fr[:, 0]
    zc_u = np.unwrap(bu[:, 4] * 2 * np.pi / L0) * L0 / (2 * np.pi)
    z_lab = fr[:, 2] + (zc_u - zc_u[0])
    w = smooth(np.gradient(z_lab, t))
    m = t >= T0_PLATEAU
    u_inf = w[m].mean()
    dx = L0 / 2 ** lv
    Re_b = u_inf * D / MU
    res.append(dict(lv=lv, col=col, t=t, z=z_lab, w=w, u=u_inf, dx=dx,
                    Ddx=D / dx, ddx=D / np.sqrt(Re_b) / dx, Re=Re_b,
                    CD=4 * G * D * (1 - 1 / 850) / (3 * u_inf ** 2),
                    chi=bu[m, 12].mean(), vol=100 * bu[m, 1].mean() / V_THEO,
                    ncol=bu.shape[1]))

u6, u7, u8 = (r["u"] for r in res)
u_rich = u8 + (u8 - u7) / 3

fig = plt.figure(figsize=(13.5, 8.0))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.36, wspace=0.32)


def style(a):
    a.grid(alpha=.25, lw=.6, color=INK3); a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_color(INK3)
    a.tick_params(colors=INK2, labelsize=11)


# (a) Trajectoires superposées
a = fig.add_subplot(gs[0, :2]); style(a)
for r in res:
    a.plot(r["t"] - r["t"][0], r["z"], lw=2.0, color=r["col"],
           label=rf"Niveau {r['lv']} ($D/\Delta x = {r['Ddx']:.0f}$)")
a.set_xlabel("t - t$_{inj}$  [u.t.]", color=INK)
a.set_ylabel(r"$z_{lab}$  [u.l.]", color=INK)
a.set_title("Position verticale $z_{lab}(t)$ en régime laminaire", color=INK, pad=8)
a.legend(frameon=False, loc="upper left", labelcolor=INK2)

# (b) u_inf vs résolution
a = fig.add_subplot(gs[0, 2]); style(a)
dd = [r["ddx"] for r in res]
uu = [r["u"] for r in res]
a.axhline(u_rich, color=INK3, ls="--", lw=1.2,
          label=rf"Extrapolation $u_\infty = {u_rich:.2f}$")
for r in res:
    a.plot(r["ddx"], r["u"], "o", ms=10, color=r["col"], mec="white", mew=2, zorder=3)
    a.annotate(f"Niveau {r['lv']}", xy=(r["ddx"], r["u"]),
               xytext=(r["ddx"] + 0.14, r["u"] + 0.10), color=r["col"],
               fontweight="bold")
a.plot(dd, uu, "-", lw=1.2, color=INK3, alpha=.5, zorder=1)
a.axvspan(0, 1.0, color="#e34948", alpha=.09, lw=0)
a.text(0.5, 11.55, "Zone sous-résolue\n($\\delta/\\Delta x < 1$)", color="#e34948",
       ha="center")
a.set_xlabel(r"Résolution de couche limite $\delta/\Delta x$", color=INK)
a.set_ylabel(r"$u_\infty$", color=INK)
a.set_title("Vitesse finale vs résolution de couche limite", color=INK, pad=8)
a.legend(frameon=False, loc="lower right", labelcolor=INK2)
a.set_xlim(0, 3.9)

# (c) Vitesse instantanée
a = fig.add_subplot(gs[1, :2]); style(a)
for r in res:
    a.plot(r["t"] - r["t"][0], r["w"], lw=1.6, color=r["col"], label=f"Niveau {r['lv']}")
a.axhline(u_rich, color=INK3, ls="--", lw=1.0)
a.set_xlabel("t - t$_{inj}$  [u.t.]", color=INK)
a.set_ylabel(r"$w = \mathrm{d}z_{lab}/\mathrm{d}t$", color=INK)
a.set_title("Profil de vitesse d'ascension instantanée", color=INK, pad=8)
a.legend(frameon=False, loc="lower right", labelcolor=INK2)
a.set_ylim(0, 15.5)

# (d) Écarts relatifs
a = fig.add_subplot(gs[1, 2]); style(a)
labels = ["Niveau 6\n vers 7", "Niveau 7\n vers 8"]
vals = [100 * (u7 / u6 - 1), 100 * (u8 / u7 - 1)]
cols = ["#2a78d6", "#008300"]
bars = a.bar(labels, vals, color=cols, width=.55, zorder=2)
for b, v in zip(bars, vals):
    a.text(b.get_x() + b.get_width() / 2, v + 0.45, f"{v:+.2f} %",
           ha="center", color=INK, fontweight="bold")
a.axhline(0, color=INK3, lw=.8)
a.set_ylabel(r"Variation de $u_\infty$  [%]", color=INK)
a.set_title("Écart relatif entre maillages", color=INK, pad=8)
a.set_ylim(-1, 15)

fig.suptitle("Étude de convergence en maillage en régime laminaire ($Bo = 1$, $Ga = 70$)",
             color=INK, y=0.99)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=150, facecolor="#fcfcfb", bbox_inches="tight")
print(f"\n-> {OUT}")
