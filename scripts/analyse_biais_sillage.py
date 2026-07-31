#!/usr/bin/env python3
"""
analyse_biais_sillage.py -- CHIFFRER la limite "la bulle remonte dans son
propre sillage" (domaine periodique).

LE PROBLEME
-----------
La reference laminaire est censee etre "une bulle en liquide au repos". Mais le
domaine est triplement periodique et la bulle parcourt 10,3 traversees : elle
finit par remonter dans le sillage qu'elle a elle-meme laisse, recycle par la
periodicite. Le liquide n'est donc plus au repos, et `u_inf = 12.27` est biaise.

Combien ? On ne peut pas supprimer l'effet sans changer de domaine, mais on peut
le MESURER, de deux facons independantes qui doivent tomber d'accord.

METHODE 1 -- prediction par notre propre courbe maitre
   `ke` de stats.dat est deja l'energie de FLUCTUATION (event logfile retranche
   la moyenne avant de sommer), donc c'est exactement l'energie du sillage
   accumule. On en tire u'_sillage = sqrt(2 ke/3) puis
       beta_sillage = u'_sillage / u_inf
   et on injecte ce beta dans la loi petit-beta ajustee sur nos 3 points fiables
   (1 - c beta^2), pour predire de combien le sillage ralentit la bulle.

METHODE 2 -- mesure directe
   La derive de u_inf(t) = d(z_lab)/dt sur le plateau, pendant que le sillage
   s'accumule.

RESULTAT (lam7_frame)
---------------------
   beta_sillage sature a 0.045 (le sillage ne s'accumule PAS indefiniment :
       ke passe de 0.048 a ~0.45 vers t=220 puis plafonne)
   methode 1 : -0.76 %
   methode 2 : -0.8 %
   => les deux tombent d'accord a 0.04 point. Le biais de sillage vaut
      **environ -0,8 % sur u_inf, et il est borne**, soit un a deux ordres de
      grandeur sous les effets mesures (-3 % a -38 %).

   Sur les points turbulents, ajouter le sillage en quadrature
   (beta_eff = sqrt(beta^2 + beta_sillage^2)) deplace beta de +4.4 % au point le
   plus faible (beta=0.15) a +0.7 % au plus fort. C'est le point beta=0.15 qui
   est le plus sensible, ce qui est logique : c'est la que le sillage propre
   pese le plus par rapport a la turbulence imposee.

=> LIMITE A ASSUMER, PAS A CACHER : chiffree, bornee, et petite devant le signal.
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
# les 3 points d'ensemble fiables (n=3 membres) de la campagne weber
BETA_OK = np.array([0.15, 0.22, 0.31])
RATIO_OK = np.array([0.97, 0.82, 0.62])

st = load_table(RUN / "stats.dat", verbose=False)
fr = load_table(RUN / "frame.dat", verbose=False)
bu = load_table(RUN / "bubble.dat", verbose=False)

t_st, ke = st[:, 0], st[:, 2]
u_wake = np.sqrt(2 * ke / 3)
beta_wake = u_wake / U_INF

# vitesse d'ascension labo, sans convention de moyenne : d(frame_z + zc)/dt
zc_u = np.unwrap(bu[:, 4] * 2 * np.pi / L0) * L0 / (2 * np.pi)
z_lab = fr[:, 2] + (zc_u - zc_u[0])
w = np.gradient(z_lab, fr[:, 0])
w = np.convolve(np.pad(w, (7, 7), mode="edge"), np.ones(15) / 15, mode="valid")

# loi petit-beta ajustee (sans terme constant : ratio -> 1 quand beta -> 0)
c = np.sum((1 - RATIO_OK) * BETA_OK ** 2) / np.sum(BETA_OK ** 4)
b_end = beta_wake[-1]
pred = 100 * c * b_end ** 2

i0, i1 = np.argmin(abs(fr[:, 0] - 220)), np.argmin(abs(fr[:, 0] - 260))
mes = 100 * (w[i1] - w[i0]) / U_INF

print(f"beta_sillage (fin de run) = {b_end:.4f}   [ke sature a {ke[-1]:.3f}]")
print(f"loi ajustee : ratio = 1 - {c:.2f} beta^2")
print(f"  methode 1 (prediction par la courbe maitre) : {-pred:+.2f} %")
print(f"  methode 2 (derive mesuree de u_inf)         : {mes:+.2f} %")
print(f"  ecart entre les deux : {abs(pred + mes):.2f} point")
print()
for bb in BETA_OK.tolist() + [0.38]:
    print(f"  beta={bb:.2f} -> beta_eff={np.hypot(bb, b_end):.4f} "
          f"({100*(np.hypot(bb, b_end)/bb - 1):+.1f} %)")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))

a = ax[0]
a.plot(t_st, beta_wake, lw=1.8, color="C0")
a.axhline(b_end, color="C3", ls="--", lw=1.0,
          label=rf"saturation $\beta_{{sillage}}={b_end:.3f}$")
a.set_xlabel("t  [u.t.]"); a.set_ylabel(r"$\beta_{sillage}=u'_{sillage}/u_\infty$")
a.set_title("Le sillage propre s'accumule… puis SATURE\n(le biais est borné, pas cumulatif)")
a.legend(fontsize=8); a.grid(alpha=.3)

a = ax[1]
bb = np.linspace(0, 0.35, 200)
a.plot(bb, 1 - c * bb ** 2, color="0.5", lw=1.4,
       label=rf"loi ajustée $1-{c:.1f}\beta^2$")
a.plot(BETA_OK, RATIO_OK, "o", ms=7, color="C0", label="points d'ensemble (n=3)")
a.axvline(b_end, color="C3", ls="--", lw=1.0)
a.plot([b_end], [1 - c * b_end ** 2], "*", ms=14, color="C3",
       label=rf"sillage propre : $-{pred:.2f}\,\%$")
a.annotate(f"mesuré indépendamment : {mes:+.1f} %", xy=(b_end, 1 - c * b_end ** 2),
           xytext=(0.09, 0.80), fontsize=8, color="C3",
           arrowprops=dict(arrowstyle="->", color="C3", lw=1.0))
a.set_xlabel(r"$\beta=u'/u_\infty$"); a.set_ylabel(r"$\bar{u}_\infty/u_\infty$")
a.set_title("Le biais de sillage placé sur notre propre courbe maîtresse")
a.legend(fontsize=8); a.grid(alpha=.3)

fig.suptitle("Limite du domaine périodique : la bulle remonte dans son propre sillage — "
             f"effet chiffré à {-pred:.1f} %", fontsize=11)
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=140)
print(f"\n-> {OUT}")
