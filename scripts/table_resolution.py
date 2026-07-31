#!/usr/bin/env python3
"""
table_resolution.py -- TABLEAU DE RESOLUTION de la campagne turbulente.

Repond a l'objection (legitime) : "vous avez valide la resolution sur le cas
LAMINAIRE ; la turbulence ajoute des petites echelles, est-elle resolue ?"
Demandee par le tuteur depuis le 2026-07-07. Aucun run necessaire : tout est
dans les stats.dat deja produits.

Deux criteres, l'un pour la turbulence, l'autre pour la bulle :

  TURBULENCE   k_max . eta >= 1.5      (critere usuel de DNS resolue)
      k_max = pi / dx = pi 2^level / L0        (nombre d'onde max du maillage)
      eta   = (nu^3/eps)^(1/4)                 (echelle de Kolmogorov, col 5)

  BULLE        D/dx  >= 20 souhaite,  et surtout la couche limite
      delta ~ D/sqrt(Re_b)  doit faire plusieurs mailles.
      La convergence mesuree (lam6/lam7/lam8) montre que delta/dx = 1.64
      suffit deja : u_inf change de 0.32 % entre lvl7 et lvl8, alors qu'il
      change de 12.3 % entre lvl6 (delta/dx = 0.87) et lvl7.

Sort un tableau texte + un tableau LaTeX pret a coller dans le rapport.
"""
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dataio import load_table

RUNS = Path("simulations/weber_ensemble")
OUT_TEX = Path("scripts/figures/table_resolution.tex")
L0, D, NU, T0 = 120.0, 16.0, 0.01 * (120.0 / (2 * np.pi)) ** 2 / 2, 220.0
U_LAM = 12.27
POINTS = [("wt05", 0.15, 5), ("wt11", 0.22, 11), ("wt21", 0.31, 21),
          ("wt25", 0.33, 25), ("wt32", 0.38, 32)]

rows = []
for tag, beta, ke in POINTS:
    st = load_table(RUNS / f"{tag}_m0" / "stats.dat", verbose=False)
    fr = load_table(RUNS / f"{tag}_m0" / "frame.dat", verbose=False)
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

kmin = min(np.pi * 2 ** 7 / L0 * r["eta"] for r in rows)
print(f"\n  minimum de k_max.eta a lvl7 sur toute la campagne : {kmin:.2f}  "
      f"(seuil 1.5)  ->  {'RESOLU' if kmin >= 1.5 else 'INSUFFISANT'}")
print(f"  marge : x{kmin/1.5:.1f}")

with open(OUT_TEX, "w", encoding="utf-8") as f:
    f.write("% genere par scripts/table_resolution.py\n")
    f.write("\\begin{tabular}{lccccccc}\n\\hline\n")
    f.write("Point & $\\beta$ & $Re_\\lambda$ & $\\varepsilon$ & $\\eta$ & "
            "$k_{max}\\eta$ & $D/\\Delta x$ & $\\delta/\\Delta x$ \\\\\n")
    f.write("\\multicolumn{5}{l}{} & \\multicolumn{3}{c}{(niveau 7)} \\\\\n\\hline\n")
    for r in rows:
        kmax = np.pi * 2 ** 7 / L0
        Re_b = r["u_inf"] * D / NU
        f.write(f"{r['tag']} & {r['beta']:.2f} & {r['relam']:.1f} & {r['eps']:.2f} & "
                f"{r['eta']:.2f} & {kmax*r['eta']:.2f} & {D/(L0/2**7):.1f} & "
                f"{D/np.sqrt(Re_b)/(L0/2**7):.2f} \\\\\n")
    f.write("\\hline\n\\end{tabular}\n")
print(f"\n-> {OUT_TEX}")
