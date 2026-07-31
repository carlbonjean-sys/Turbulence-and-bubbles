#!/usr/bin/env python3
"""
plot_mechanism.py -- Analyse du mécanisme de ralentissement (échantillonnage préférentiel vs traînée non linéaire).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
from pathlib import Path

RUNS = [
    ("laminaire", "simulations/lam7_frame",            0.00),
    ("wt05",      "simulations/betasweep/wt05_m0_bub",  0.15),
    ("wt21",      "simulations/betasweep/wt21_m0_bub",  0.31),
    ("wt32",      "simulations/betasweep/wt32_m0_bub",  0.38),
]
OUT = Path("scripts/figures/mechanism.png")
U_LAM = 12.27
T0 = 220.0

CONE = [4, 5, 6, 7]
RBAND = [1.25, 1.75, 2.25, 2.75]
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a84"
COL = {"laminaire": "#8a8a84", "wt05": "#2a78d6", "wt21": "#eb6834", "wt32": "#008300"}

def load(path):
    f = Path(path) / "mechanism.dat"
    if not f.is_file():
        return None
    d = []
    with open(f) as fp:
        for l in fp:
            if l.startswith("#"): continue
            p = l.split()
            if len(p) >= 17:
                d.append([float(x) for x in p[:17]])
    return np.array(d) if d else None

data = {}
for lab, p, b in RUNS:
    a = load(p)
    if a is not None and len(a) > 5:
        sel = a[:, 0] >= T0
        if sel.sum() > 2:
            data[lab] = dict(raw=a[sel], beta=b)

if not data or "laminaire" not in data:
    print("[mechanism] Donnees de mécanisme non trouvées dans betasweep. Ignoré.")
    sys.exit(0)

lam = data["laminaire"]["raw"]
w_lam = (lam[:, 2] - lam[:, 3]).mean()
uf_lam_r = np.array([(lam[:, col] - lam[:, 3]).mean() for col in CONE])

res = {}
for lab, d in data.items():
    raw = d["raw"]
    w_lab = (raw[:, 2] - raw[:, 3]).mean()
    uf_lab_r = np.array([(raw[:, col] - raw[:, 3]).mean() for col in CONE])
    du = U_LAM - w_lab
    sig_pref = uf_lab_r - uf_lam_r * (w_lab / w_lam)
    res[lab] = dict(beta=d["beta"], w_lab=w_lab, du=du, uf_lab_r=uf_lab_r, sig_pref=sig_pref)

plt.rcParams.update({'font.size': 11})
fig, ax = plt.subplots(1, 2, figsize=(13, 5))

for a in ax:
    a.grid(alpha=.25, lw=.6, color=INK3)
    a.set_axisbelow(True)
    for s in ("top", "right"): a.spines[s].set_visible(False)
    for s in ("left", "bottom"): a.spines[s].set_color(INK3)

a = ax[0]
for lab, r in res.items():
    a.plot(RBAND, r["uf_lab_r"], "o-", color=COL[lab], label=f"{lab} (beta={r['beta']})")
a.set_xlabel("r / D (distance radiale)")
a.set_ylabel("Champ moyen de vitesse fluide Uf_lab")
a.legend(frameon=True, facecolor="white")

a = ax[1]
for lab, r in res.items():
    if lab == "laminaire": continue
    a.plot(RBAND, r["sig_pref"], "s--", color=COL[lab], label=f"{lab} (signal préférentiel)")
a.set_xlabel("r / D")
a.set_ylabel("Signal d'entraînement préférentiel")
a.legend(frameon=True, facecolor="white")

plt.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT, dpi=300)
plt.close()
print("Figure du mécanisme générée avec succès dans scripts/figures/")
