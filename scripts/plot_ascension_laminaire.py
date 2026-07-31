#!/usr/bin/env python3
"""
plot_ascension_laminaire.py -- Validation en régime laminaire.

Vitesse d'ascension et trajectoire de la bulle en liquide au repos.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

RUN   = Path(sys.argv[1] if len(sys.argv) > 1 else "simulations/lam7_frame")
OUT   = Path("scripts/figures/ascension_laminaire.png")

G, R0, RHO1, RHOR, WIDTH, MAXLEVEL = 4.0, 8.0, 1.0, 850.0, 120.0, 7
RHO2 = RHO1 / RHOR
MU1  = 0.01 * (WIDTH / (2 * np.pi)) ** 2 / 2.0
D    = 2 * R0

fr = np.loadtxt(RUN / "frame.dat")
bu = np.loadtxt(RUN / "bubble.dat")
t        = fr[:, 0]
frame_uz = fr[:, 1]
frame_z  = fr[:, 2]
u_h      = np.hypot(bu[:, 5], bu[:, 6])
t0       = t[0]

def unwrap(c):
    return np.unwrap(c * 2 * np.pi / WIDTH) * WIDTH / (2 * np.pi)

zc_u  = unwrap(bu[:, 4])
z_lab = frame_z + (zc_u - zc_u[0])

def smooth(y, n=15):
    k = np.ones(n) / n
    return np.convolve(np.pad(y, (n // 2, n - 1 - n // 2), mode="edge"), k, mode="valid")

w_lab = smooth(np.gradient(z_lab, t))

plateau = t > t0 + 60
u_inf   = w_lab[plateau].mean()

trans = t < t0 + 25
def fit_tau(tt, ww, uinf):
    ok = (ww > 0) & (ww < 0.95 * uinf)
    y  = np.log(1.0 - ww[ok] / uinf)
    return -1.0 / np.polyfit(tt[ok] - t0, y, 1)[0]

tau_fit  = fit_tau(t[trans], w_lab[trans], u_inf)
tau_theo = u_inf / (2 * G)
a_init   = (RHO1 - RHO2) * G / (RHO2 + RHO1 / 2)

Re_b = RHO1 * u_inf * D / MU1
C_D  = 4 * G * D * (1 - RHO2 / RHO1) / (3 * u_inf ** 2)
dx   = WIDTH / 2 ** MAXLEVEL

fig, ax = plt.subplots(2, 2, figsize=(13.5, 8.5))
tt = t - t0

# (1) Vue d'ensemble
a = ax[0, 0]
a.plot(tt, w_lab, lw=1.2, color="C0", label=r"$w_{lab} = \mathrm{d}z_{lab}/\mathrm{d}t$")
a.plot(tt, frame_uz, lw=0.9, color="C7", ls="--", label=r"Repère mobile $u_z$")
a.axhline(u_inf, color="k", ls=":", lw=1.2, label=rf"$u_\infty = {u_inf:.2f}$")
a.plot(tt, u_inf * (1 - np.exp(-tt / tau_fit)), color="C3", lw=1.6, ls="-.",
       label=rf"Modèle $u_\infty(1-e^{{-t/\tau}})$, $\tau = {tau_fit:.2f}$")
a.set_xlabel("t - t$_{inj}$  [u.t.]"); a.set_ylabel(r"$w$  [u.l.]")
a.set_title("Évolution temporelle de la vitesse d'ascension")
a.legend(frameon=False); a.grid(alpha=.3)

# (2) Zoom transitoire
a = ax[0, 1]
m = tt < 20
a.plot(tt[m], w_lab[m], "o-", ms=3, lw=1.0, color="C0", label="DNS")
a.plot(tt[m], u_inf * (1 - np.exp(-tt[m] / tau_fit)), color="C3", lw=1.8, ls="-.",
       label=rf"Ajustement $\tau = {tau_fit:.2f}$")
a.plot(tt[m], u_inf * (1 - np.exp(-tt[m] / tau_theo)), color="C2", lw=1.4, ls=":",
       label=rf"Masse ajoutée $\tau = {tau_theo:.2f}$")
a.axhline(u_inf, color="k", ls=":", lw=1.0)
a.set_xlabel("t - t$_{inj}$  [u.t.]"); a.set_ylabel(r"$w$  [u.l.]")
a.set_title("Régime transitoire initial")
a.legend(frameon=False); a.grid(alpha=.3)

# (3) Test log
a = ax[1, 0]
ok = (tt < 25) & (w_lab > 0) & (w_lab < 0.95 * u_inf)
a.plot(tt[ok], np.log(1 - w_lab[ok] / u_inf), "o", ms=3, color="C0", label="DNS")
a.plot(tt[ok], -tt[ok] / tau_fit, color="C3", lw=1.8, ls="-.", label=rf"Pente $-1/\tau$, $\tau = {tau_fit:.2f}$")
a.set_xlabel("t - t$_{inj}$  [u.t.]"); a.set_ylabel(r"$\ln(1 - w/u_\infty)$")
a.set_title(r"Vérification du comportement exponentiel $\ln(1 - w/u_\infty)$")
a.legend(frameon=False); a.grid(alpha=.3)

# (4) Trajectoire labo
a = ax[1, 1]
a.plot(tt, z_lab, color="C0", lw=1.4, label=r"$z_{lab} = z_{frame} + z_c$")
a.plot(tt, u_inf * tt, color="k", ls=":", lw=1.0, label=rf"$u_\infty t$")
a.set_xlabel("t - t$_{inj}$  [u.t.]"); a.set_ylabel(r"$z_{lab}$  [u.l.]", color="C0")
a.set_title("Trajectoire verticale dans le repère du laboratoire")
a.grid(alpha=.3); a.legend(frameon=False, loc="upper left")
a2 = a.twinx()
a2.plot(tt, u_h, color="C1", lw=0.8, alpha=.7)
a2.set_ylabel(r"Dérive horizontale $|u_h|$", color="C1")

fig.suptitle(f"Ascension laminaire d'une bulle ($Ga = 70$, $Bo = 1$)", fontsize=12)
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=140)
print(f"\n-> {OUT}")
