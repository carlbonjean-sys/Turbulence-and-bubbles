#!/usr/bin/env python3
"""
analyse_weber_ensemble.py -- Dépouillement des simulations de remontee (moyennes d'ensemble).
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R0     = 8.0
D      = 2.0 * R0
G      = 4.0
V_THEO = 4.0/3.0*np.pi*R0**3
U_LAM  = 12.27
FR_FAC = U_LAM / np.sqrt(G*D)

CF = {"t": 0, "uz": 1, "z": 2, "ub": 3}
CB = {"t": 0, "volume": 1, "chi": 12}

POINTS = [("beta0150", 5, 0.15), ("beta0220", 11, 0.22), ("beta0310", 21, 0.31),
          ("beta0330", 25, 0.33), ("beta0380", 32, 0.38)]
MEMBERS = [0, 1, 2, 3, 4]

from dataio import load_table

def load(path):
    return load_table(path)

def plateau_uz(frame, t0):
    if frame is None or len(frame) < 3:
        return None
    t = frame[:, CF["t"]]
    m = t >= t0
    if m.sum() < 3:
        return None
    uz = frame[m, CF["uz"]]
    return dict(mean=float(np.mean(uz)), std=float(np.std(uz)),
                npts=int(m.sum()), tend=float(t[-1]))

def vol_check(bub, t0):
    if bub is None or len(bub) < 3:
        return None
    t = bub[:, CB["t"]]
    m = t >= t0
    if m.sum() < 1:
        m = np.ones(len(t), bool)
    return dict(vfrac_min=float(np.min(bub[m, CB["volume"]]))/V_THEO,
                chi_max=float(np.max(bub[m, CB["chi"]])))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--t0", type=float, default=220.0,
                    help="debut du plateau frame_uz")
    args = ap.parse_args()

    base = os.path.join(args.root, "simulations", "betasweep")
    figdir = os.path.join(args.root, "scripts", "figures")
    os.makedirs(figdir, exist_ok=True)

    rows = []
    for lab, ke, beta in POINTS:
        members = []
        for m in MEMBERS:
            d = os.path.join(base, f"{lab}_m{m}_bub")
            pl = plateau_uz(load(os.path.join(d, "frame.dat")), args.t0)
            vc = vol_check(load(os.path.join(d, "bubble.dat")), args.t0)
            if pl is not None:
                members.append(dict(m=m, uz=pl["mean"], std=pl["std"],
                                    npts=pl["npts"], tend=pl["tend"], vc=vc))
        vals = np.array([mm["uz"] for mm in members])
        n = len(vals)
        if n >= 1:
            mean = float(np.mean(vals))
            sd = float(np.std(vals, ddof=1)) if n >= 2 else float("nan")
            sem = sd/np.sqrt(n) if n >= 2 else float("nan")
            ratio = mean/U_LAM
            red = 100.0*(ratio-1.0)
        else:
            mean = sd = sem = ratio = red = float("nan")
        rows.append(dict(lab=lab, ke=ke, beta=beta, frp=beta*FR_FAC,
                         members=members, n=n, mean=mean, sd=sd, sem=sem,
                         ratio=ratio, red=red))

    csv = os.path.join(figdir, "weber_ensemble.csv")
    with open(csv, "w", encoding="utf-8") as f:
        f.write("Weber,KE,beta,Fr_prime,n_membres,uz_moy,uz_std,SEM,ratio,reduction_pct,membres_uz\n")
        for r in rows:
            mu = ";".join(f"m{mm['m']}:{mm['uz']:.3f}" for mm in r["members"])
            f.write(f"{r['lab']},{r['ke']},{r['beta']},{r['frp']:.3f},{r['n']},"
                    f"{r['mean']:.3f},{r['sd']:.3f},{r['sem']:.3f},{r['ratio']:.3f},"
                    f"{r['red']:.1f},{mu}\n")

    fig, ax = plt.subplots(figsize=(8, 6))
    bb = np.linspace(0, 0.42, 100)
    ax.plot(bb, 1.0 - bb**2, "k--", lw=1.2, alpha=0.7, label=r"$1-\beta^2$")
    ax.axhline(1.0, color="gray", lw=0.8, ls=":")

    for r in rows:
        if not np.isfinite(r["ratio"]):
            continue
        if r["n"] >= 2:
            ax.errorbar(r["beta"], r["ratio"], yerr=r["sem"]/U_LAM,
                        fmt="o", ms=9, color="C0", capsize=4, zorder=5)
        else:
            ax.plot(r["beta"], r["ratio"], "o", ms=9, mfc="white", mec="C1", mew=1.8, zorder=5)

    ax.set_xlabel(r"$\beta = u'/u_\infty$")
    ax.set_ylabel(r"$\bar u_{\infty,\mathrm{turb}} / u_\infty$")
    ax.set_xlim(0, 0.42)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = os.path.join(figdir, "weber_rise_prelim.png")
    plt.savefig(out, dpi=150)
    plt.close()

if __name__ == "__main__":
    main()
