#!/usr/bin/env python3
"""
analyse_biais_sillage.py -- Analyse du biais de sillage propre en domaine périodique.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dataio import load_table

RUN = Path("simulations/lam7_frame")
OUT = Path("scripts/figures/biais_sillage.png")
U_INF, L0 = 12.27, 120.0

BETA_OK = np.array([0.15, 0.22, 0.31])
RATIO_OK = np.array([0.97, 0.82, 0.62])

st = load_table(RUN / "stats.dat", verbose=False)
fr = load_table(RUN / "frame.dat", verbose=False)
bu = load_table(RUN / "bubble.dat", verbose=False)

if st is None or fr is None or bu is None:
    print("[sillage] Fichiers manquants dans lam7_frame.")
    sys.exit(0)

t_st, ke = st[:, 0], st[:, 2]
u_wake = np.sqrt(2 * ke / 3)
beta_wake = u_wake / U_INF

zc_u = np.unwrap(bu[:, 4] * 2 * np.pi / L0) * L0 / (2 * np.pi)
z_lab = fr[:, 2] + (zc_u - zc_u[0])
w = np.gradient(z_lab, fr[:, 0])
w = np.convolve(np.pad(w, (7, 7), mode="edge"), np.ones(15) / 15, mode="valid")

c = np.sum((1 - RATIO_OK) * BETA_OK ** 2) / np.sum(BETA_OK ** 4)
b_end = beta_wake[-1]
pred = 100 * c * b_end ** 2

i0, i1 = np.argmin(abs(fr[:, 0] - 220)), np.argmin(abs(fr[:, 0] - 260))
mes = 100 * (w[i1] - w[i0]) / U_INF

fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))

a = ax[0]
a.plot(t_st, beta_wake, lw=1.8, color="C0")
a.axhline(b_end, color="C3", ls="--", lw=1.0, label=rf"$\beta_{{sillage}}={b_end:.3f}$")
a.set_xlabel("t [u.t.]")
a.set_ylabel(r"$\beta_{sillage}$")
a.legend(fontsize=9)
a.grid(alpha=.3)

a = ax[1]
bb = np.linspace(0, 0.35, 200)
a.plot(bb, 1 - c * bb ** 2, color="0.5", lw=1.4, label=rf"$1-{c:.1f}\beta^2$")
a.plot(BETA_OK, RATIO_OK, "o", ms=7, color="C0", label="DNS")
a.axvline(b_end, color="C3", ls="--", lw=1.0)
a.plot([b_end], [1 - c * b_end ** 2], "*", ms=14, color="C3", label=rf"sillage : $-{pred:.2f}\,\%$")
a.set_xlabel(r"$\beta$")
a.set_ylabel(r"$\bar{u}_\infty/u_\infty$")
a.legend(fontsize=9)
a.grid(alpha=.3)

plt.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=140)
plt.close()
print("Figure de biais de sillage générée avec succès.")
