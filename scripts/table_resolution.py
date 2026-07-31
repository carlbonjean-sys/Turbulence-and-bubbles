#!/usr/bin/env python3
"""
table_resolution.py -- Tableau de résolution de la campagne turbulente (DNS).
"""
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dataio import load_table

RUNS = Path("simulations/betasweep")
L0, D, NU, T0 = 120.0, 16.0, 0.01 * (120.0 / (2 * np.pi)) ** 2 / 2, 220.0
U_LAM = 12.27
POINTS = [("wt05", 0.15, 5), ("wt11", 0.22, 11), ("wt21", 0.31, 21),
          ("wt25", 0.33, 25), ("wt32", 0.38, 32)]

rows = []
for tag, beta, ke in POINTS:
    st_path = RUNS / f"{tag}_m0_bub" / "stats.dat"
    fr_path = RUNS / f"{tag}_m0_bub" / "frame.dat"
    st = load_table(st_path, verbose=False)
    fr = load_table(fr_path, verbose=False)
    if st is None or fr is None:
        continue
    m = st[:, 0] >= T0
    eps, kel, relam, eta = (st[m, 1].mean(), st[m, 2].mean(),
                            st[m, 3].mean(), st[m, 4].mean())
    u_inf = fr[fr[:, 0] >= T0, 1].mean()
    rows.append(dict(tag=tag, beta=beta, ke=ke, eps=eps, k=kel, relam=relam,
                     eta=eta, u_inf=u_inf, up=np.sqrt(2 * kel / 3)))

print(f"{'pt':6} {'beta':>5} {'Re_l':>6} {'eps':>7} {'eta':>6} "
      f"|{'  k_max.eta lvl7':>16} {'lvl8':>7} |{'  D/dx7':>8} {'d/dx7':>7} {'d/dx8':>7}")
for r in rows:
    line = f"{r['tag']:6} {r['beta']:5.2f} {r['relam']:6.1f} {r['eps']:7.3f} {r['eta']:6.3f} |"
    for lv in (7, 8):
        kmax = np.pi * 2 ** lv / L0
        line += f"{kmax*r['eta']:16.2f}" if lv == 7 else f"{kmax*r['eta']:7.2f}"
    dx7, dx8 = L0 / 2 ** 7, L0 / 2 ** 8
    Re_b = r["u_inf"] * D / NU
    delta = D / np.sqrt(Re_b)
    line += f" |{D/dx7:8.1f} {delta/dx7:7.2f} {delta/dx8:7.2f}"
    print(line)

if rows:
    kmin = min(np.pi * 2 ** 7 / L0 * r["eta"] for r in rows)
    print(f"\nMinimum k_max.eta (lvl7) = {kmin:.2f} (seuil 1.5 -> {'RESOLU' if kmin >= 1.5 else 'NON RESOLU'})")
