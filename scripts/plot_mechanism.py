#!/usr/bin/env python3
"""
plot_mechanism.py -- POURQUOI la turbulence ralentit la bulle ?
Depouille mechanism.dat (src/mechanism.c) : echantillonnage preferentiel (A) vs
trainee non lineaire (B).

CADRE
-----
Vitesses reconstruites en repere LABO depuis les colonnes STOCKEES (repere mobile)
via frame_uz = -uz_domaine :
  - montee LAB de la bulle :        W_lab = ub - uz_far
  - vitesse du fluide vue (cone) :  Uf_lab(r) = cone(r) - uz_far
Verifie en laminaire : W_lab ~ 12.2 = u_inf, Uf_lab(3D) ~ +0.53 = ecoulement induit.

DISCRIMINANT
------------
La reduction de vitesse (lab) vaut  dU = u_inf_lam - W_lab_turb.
  (A) Echantillonnage preferentiel : la bulle est PORTEE par du fluide qui descend
      -> Uf_lab(bulle) ~ -dU (negatif, du meme ordre que la reduction).
  (B) Trainee non lineaire : le fluide ambiant est en moyenne nul, la reduction
      vient de la trainee -> Uf_lab(bulle) ~ 0 (bien plus petit que dU).

MAIS Uf_lab contient l'ECOULEMENT INDUIT par la bulle (mesure en laminaire : +0.53
a 3D, proportionnel a la vitesse de glissement donc a la montee). On le retranche,
remis a l'echelle de la montee de chaque run :
  signal_preferentiel(r) = Uf_lab_turb(r) - Uf_lab_lam(r) * (W_lab_turb / u_inf_lam)
  -> ~0  => trainee non lineaire (B)
  -> ~ -dU => echantillonnage preferentiel (A)
  -> intermediaire => les deux contribuent (on chiffre la part de dU expliquee).

TEST COMPLEMENTAIRE (B) : la variance de la vitesse de glissement
  slip(t) = W_lab(t) - Uf_lab(t)  ; sigma_slip/‹slip› grand => trainee non lineaire
credible (fluctuations fortes de glissement).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os
from pathlib import Path

RUNS = [   # label, dossier, beta
    ("laminaire", "simulations/lam7_frame",            0.00),
    ("wt05",      "simulations/weber_fix/wt05_g4_bub",  0.15),
    ("wt21",      "simulations/weber_fix/wt21_g4_bub",  0.31),
    ("wt32",      "simulations/weber_fix/wt32_g4_bub",  0.38),
]
OUT = Path("scripts/figures/mechanism.png")
U_LAM = 12.27
T0 = 220.0
# colonnes de mechanism.dat : label z_c ub uz_far c1..c4 h1..h4 f1..f4 Vb
CONE = [4, 5, 6, 7]      # bandes cone r/D = [1,1.5][1.5,2][2,2.5][2.5,3]
RBAND = [1.25, 1.75, 2.25, 2.75]
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a84"
COL = {"laminaire": "#8a8a84", "wt05": "#2a78d6", "wt21": "#eb6834", "wt32": "#008300"}


def load(path):
    f = Path(path) / "mechanism.dat"
    if not f.is_file():
        return None
    rows = []
    for line in open(f):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        p = s.split()
        try:
            t = float(p[0]); vals = [float(x) for x in p[1:]]
        except ValueError:
            continue
        if len(vals) >= 16:
            rows.append([t] + vals)
    return np.array(rows) if rows else None


data = {}
for lab, d, beta in RUNS:
    a = load(d)
    if a is None:
        print(f"[attente] {lab} : mechanism.dat absent ({d})")
        continue
    t = a[:, 0]
    m = t >= T0
    if m.sum() < 3:
        m = t >= t.max() - 40
    ub = a[m, 2]; uz_far = a[m, 3]
    W_lab = ub - uz_far                       # montee LAB
    Uf = a[m, :][:, CONE] - uz_far[:, None]   # fluide vu (lab), 4 bandes cone
    data[lab] = dict(beta=beta, W=W_lab.mean(), Wstd=W_lab.std(),
                     Uf=Uf.mean(axis=0), Uf_t=Uf, W_t=W_lab, n=m.sum())

if "laminaire" not in data:
    print("\n>> Le run laminaire (baseline) est requis. Relancer quand present.")
    if not data:
        sys.exit(0)

# ---- verdict par run
print(f"\n{'run':10} {'beta':>5} {'W_lab':>7} {'dU':>6} {'Uf(3D)':>7} "
      f"{'induit*':>8} {'signal':>7} {'part dU':>8}")
lam = data.get("laminaire")
res = []
for lab, d in data.items():
    if lab == "laminaire":
        print(f"{lab:10} {d['beta']:5.2f} {d['W']:7.2f} {'--':>6} "
              f"{d['Uf'][-1]:+7.2f} {'(ref)':>8} {'--':>7} {'--':>8}")
        continue
    dU = U_LAM - d["W"]                        # reduction de vitesse (lab)
    Uf_far = d["Uf"][-1]                       # cone r=2.5-3D
    induced = lam["Uf"][-1] * (d["W"] / lam["W"])   # induit, remis a l'echelle
    signal = Uf_far - induced                  # echantillonnage preferentiel isole
    part = 100 * (-signal) / dU if dU > 0 else 0.
    res.append((lab, d["beta"], dU, Uf_far, induced, signal, part))
    print(f"{lab:10} {d['beta']:5.2f} {d['W']:7.2f} {dU:6.2f} {Uf_far:+7.2f} "
          f"{induced:+8.2f} {signal:+7.2f} {part:7.0f}%")

print("\nLecture : 'signal' = vitesse du fluide vue par la bulle une fois l'ecoulement")
print("induit retranche. ~0 => trainee non lineaire ; ~ -dU => echantillonnage preferentiel.")
print("'part dU' = fraction de la reduction expliquee par l'echantillonnage preferentiel.")

# ---- figure
fig, ax = plt.subplots(1, 3, figsize=(14, 4.6))
for a in ax:
    a.grid(alpha=.25, lw=.6, color=INK3); a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_color(INK3)
    a.tick_params(colors=INK2, labelsize=9)

# (a) fluide vu vs rayon, tous runs
a = ax[0]
for lab, d in data.items():
    a.plot(RBAND, d["Uf"], "o-", color=COL[lab], lw=1.6, ms=6,
           label=f"{lab} (β={d['beta']:.2f})")
a.axhline(0, color=INK3, lw=.8, ls=":")
a.set_xlabel(r"rayon de la coquille  $r/D$ (cône sup.)", fontsize=10, color=INK)
a.set_ylabel(r"$U_f^{lab}$ vue par la bulle", fontsize=10, color=INK)
a.set_title("Vitesse du fluide vue par la bulle", fontsize=11, color=INK, pad=8)
a.legend(fontsize=8, frameon=False, labelcolor=INK2)

# (b) signal preferentiel vs reduction
a = ax[1]
if res:
    betas = [r[1] for r in res]
    dUs = [r[2] for r in res]
    sigs = [-r[5] for r in res]        # -signal = downflow sampled
    a.plot(betas, dUs, "s-", color="#e34948", lw=1.8, ms=8,
           label="réduction totale $\\Delta U$")
    a.plot(betas, sigs, "o-", color="#4a3aa7", lw=1.8, ms=8,
           label="expliquée par échantillonnage")
    a.axhline(0, color=INK3, lw=.8)
    a.set_xlabel(r"$\beta$", fontsize=10, color=INK)
    a.set_ylabel("vitesse [u. code]", fontsize=10, color=INK)
    a.set_title("Échantillonnage préférentiel\nvs réduction totale", fontsize=11,
                color=INK, pad=8)
    a.legend(fontsize=8.5, frameon=False, labelcolor=INK2, loc="upper left")

# (c) fluctuations de glissement (test trainee non lineaire)
a = ax[2]
for lab, d in data.items():
    if lab == "laminaire":
        continue
    Uf_t = d["Uf_t"][:, -1]
    slip = d["W_t"] - Uf_t
    a.plot([d["beta"]], [slip.std() / abs(slip.mean())], "o", ms=10,
           color=COL[lab], mec="white", mew=1.5)
a.set_xlabel(r"$\beta$", fontsize=10, color=INK)
a.set_ylabel(r"$\sigma_{slip}/\langle slip\rangle$", fontsize=10, color=INK)
a.set_title("Fluctuations du glissement\n(fortes → traînée non linéaire)",
            fontsize=11, color=INK, pad=8)

fig.suptitle("Mécanisme du ralentissement : échantillonnage préférentiel "
             "vs traînée non linéaire", fontsize=12.5, color=INK, y=1.0)
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=150, facecolor="#fcfcfb")
print(f"\n-> {OUT}")
