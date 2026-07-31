#!/usr/bin/env python3
"""
plot_conv_turb.py -- Convergence en maillage en turbulence (lvl 6/7/8).

Analyse de la sensibilité au maillage de la vitesse d'ascension et des
statistiques de turbulence.
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

OUT = Path("scripts/figures/convergence_turbulente.png")
L0, R0 = 120.0, 8.0
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a84"
COL = {6: "#2a78d6", 7: "#eb6834", 8: "#008300"}
WIN = (180, 209)             # Fenêtre commune d'analyse

d = {}
for lv in (6, 7, 8):
    fr = load_table(f"simulations/convturb/lvl{lv}/frame.dat", verbose=False)
    bu = load_table(f"simulations/convturb/lvl{lv}/bubble.dat", verbose=False)
    st = load_table(f"simulations/convturb/lvl{lv}/stats.dat", verbose=False)
    d[lv] = dict(fr=fr, bu=bu, st=st, tmax=fr[-1, 0])


def uinterp(t, c, td):
    cu = np.unwrap(c * 2 * np.pi / L0) * L0 / (2 * np.pi)
    return np.interp(td, t, cu)


tc = np.linspace(161, WIN[1], 500)


def hdiv(a, b):
    dx = uinterp(a[:, 0], a[:, 2], tc) - uinterp(b[:, 0], b[:, 2], tc)
    dy = uinterp(a[:, 0], a[:, 3], tc) - uinterp(b[:, 0], b[:, 3], tc)
    return np.hypot(dx - dx[0], dy - dy[0])


d67 = hdiv(d[6]["bu"], d[7]["bu"])
d78 = hdiv(d[7]["bu"], d[8]["bu"])

for lv in (6, 7, 8):
    fr = d[lv]["fr"]
    m = (fr[:, 0] >= WIN[0]) & (fr[:, 0] <= WIN[1])
    d[lv]["u"] = fr[m, 1].mean()
    d[lv]["ustd"] = fr[m, 1].std()
    st = d[lv]["st"]
    ms = (st[:, 0] >= WIN[0]) & (st[:, 0] <= WIN[1])
    d[lv]["rel"] = st[ms, 3].mean()
    d[lv]["eta"] = st[ms, 4].mean()
    d[lv]["keta"] = np.pi * 2 ** lv / L0 * d[lv]["eta"]

print(f"{'lvl':>3} {'tmax':>5} {'u_inf':>7} {'std':>6} {'Re_lam':>7} {'k_max.eta':>10}")
for lv in (6, 7, 8):
    print(f"{lv:>3} {d[lv]['tmax']:>5.0f} {d[lv]['u']:>7.2f} {d[lv]['ustd']:>6.2f} "
          f"{d[lv]['rel']:>7.1f} {d[lv]['keta']:>10.2f}")

fig = plt.figure(figsize=(13.5, 8.0))
gs = fig.add_gridspec(2, 3, hspace=0.36, wspace=0.32)


def style(a):
    a.grid(alpha=.25, lw=.6, color=INK3); a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_color(INK3)
    a.tick_params(colors=INK2, labelsize=11)


# (a) Trajectoires verticales (montée)
ax = fig.add_subplot(gs[0, :2]); style(ax)
for lv in (6, 7, 8):
    fr = d[lv]["fr"]
    lab = f"Niveau {lv}"
    ax.plot(fr[:, 0], fr[:, 2], lw=1.8, color=COL[lv], label=lab)
ax.axvspan(WIN[0], WIN[1], color=INK3, alpha=.08, lw=0)
ax.set_xlabel("t  [u.t.]", color=INK)
ax.set_ylabel(r"$z_{lab}$  [u.l.]", color=INK)
ax.set_title("Position verticale de la bulle $z_{lab}(t)$", color=INK, pad=8)
ax.legend(frameon=False, labelcolor=INK2, loc="upper left")

# (b) u_inf par niveau
ax = fig.add_subplot(gs[0, 2]); style(ax)
for lv in (6, 7, 8):
    ax.errorbar([lv], [d[lv]["u"]], yerr=[d[lv]["ustd"] / np.sqrt(30 / 14)],
                fmt="o", ms=10, color=COL[lv], mec="white", mew=2, capsize=5,
                elinewidth=1.6, zorder=3)
ax.axhspan(d[7]["u"] - .3, d[7]["u"] + .3, color=COL[7], alpha=.08, lw=0)
ax.set_xticks([6, 7, 8]); ax.set_xticklabels(["Niveau 6", "Niveau 7", "Niveau 8"])
ax.set_xlabel(r"Niveau de raffinement", color=INK)
ax.set_ylabel(r"$u_\infty$ sur [180, 209]", color=INK)
ax.set_title(r"Vitesse moyenne d'ascension $u_\infty$", color=INK, pad=8)
ax.set_xlim(5.5, 8.5)

# (c) Divergence des trajectoires
ax = fig.add_subplot(gs[1, 0]); style(ax)
ax.semilogy(tc - tc[0], np.maximum(d67, 1e-2), lw=1.6, color="#4a3aa7",
            label="Niveau 6 vs 7")
ax.semilogy(tc - tc[0], np.maximum(d78, 1e-2), lw=1.6, color="#e34948",
            label="Niveau 7 vs 8")
ax.axhline(2 * R0, color=INK3, ls=":", lw=1.2)
ax.text(2, 2 * R0 * 1.2, "Diamètre $D$", color=INK2)
ax.set_xlabel("t - t$_{inj}$  [u.t.]", color=INK)
ax.set_ylabel(r"$|\Delta_h|$  [u.l.]", color=INK)
ax.set_title(r"Écart horizontal $|\Delta_h(t)|$", color=INK, pad=8)
ax.legend(frameon=False, labelcolor=INK2, loc="lower right")

# (d) Statistiques de turbulence
ax = fig.add_subplot(gs[1, 1]); style(ax)
for lv in (6, 7, 8):
    st = d[lv]["st"]
    ax.plot(st[:, 0], st[:, 3], lw=1.5, color=COL[lv], label=f"Niveau {lv}")
ax.set_xlabel("t  [u.t.]", color=INK)
ax.set_ylabel(r"$Re_\lambda$", color=INK)
ax.set_title(r"Nombre de Reynolds de Taylor $Re_\lambda$", color=INK, pad=8)
ax.legend(frameon=False, labelcolor=INK2)

# (e) Résolution turbulente
ax = fig.add_subplot(gs[1, 2]); style(ax)
kk = [d[lv]["keta"] for lv in (6, 7, 8)]
bars = ax.bar(["Niveau 6", "Niveau 7", "Niveau 8"], kk, color=[COL[6], COL[7], COL[8]],
              width=.55, zorder=2)
for b, v in zip(bars, kk):
    ax.text(b.get_x() + b.get_width() / 2, v + .25, f"{v:.2f}", ha="center",
            color=INK, fontweight="bold")
ax.axhline(1.5, color="#e34948", ls="--", lw=1.4, label="Seuil DNS (1.5)")
ax.set_ylabel(r"$k_{max}\,\eta$", color=INK)
ax.set_title(r"Critère de résolution spatiale $k_{max}\eta$", color=INK, pad=8)
ax.legend(frameon=False, labelcolor=INK2, loc="upper left")
ax.set_ylim(0, max(kk) * 1.25)

fig.suptitle("Convergence en maillage en écoulement turbulent ($We_t \\approx 1$, $\\beta = 0{,}22$)",
             color=INK, y=0.99)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=150, facecolor="#fcfcfb", bbox_inches="tight")
print(f"\n-> {OUT}")
