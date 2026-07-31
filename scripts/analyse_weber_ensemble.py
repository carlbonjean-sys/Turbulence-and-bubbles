#!/usr/bin/env python3
# ============================================================================
#  analyse_weber_ensemble.py -- Depouillement de la campagne "weber_rise"
#  (effet de la turbulence sur la vitesse de montee, gravite periodique CORRIGEE)
# ----------------------------------------------------------------------------
#  Cadre : Bo=1, Ga~70 (g=4), balayage de l'intensite turbulente beta = u'/u_inf.
#  Chaque point Weber a plusieurs MEMBRES d'ensemble (positions d'injection
#  decalees, argv[11]=MEMBER) -> realisations turbulentes independantes du meme
#  precurseur. On POOLE les membres par Weber.
#
#  La vitesse d'ascension est le PLATEAU de frame_uz (u_inf dans le repere mobile
#  move_frame), moyenne sur t >= T0 (defaut 220). Ratio a la reference laminaire
#  u_inf_lam = 12.24 (run lam7_frame, job 38772).  ->  reduction = 1 - ratio.
#
#  Donnees attendues (rapatriees par le fetch tar depuis Kairos) :
#    simulations/weber_ensemble/wtXX_mN/frame.dat   (t frame_uz frame_z ub_repere)
#    simulations/weber_ensemble/wtXX_mN/bubble.dat  (controle volume/chi)
#
#  Sorties (scripts/figures/) : weber_rise_prelim.png (ratio vs beta, barres
#  d'erreur d'ensemble), weber_ensemble.csv, weber_ensemble_tab.tex.
#
#  Usage : python scripts/analyse_weber_ensemble.py [--root .] [--t0 220]
# ============================================================================
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- constantes du cas (unites code) ----
R0     = 8.0
D      = 2.0 * R0
G      = 4.0                       # |g| des runs weber (GRAV=4, Ga~70)
V_THEO = 4.0/3.0*np.pi*R0**3       # 2144.66
U_LAM  = 12.27   # u_inf laminaire de reference (lam7_frame, job 38772).
                 # MAJ 2026-07-22 : 12.24 -> 12.27. L'ancienne valeur etait le
                 # plateau de frame_uz SEUL, qui oublie le residuel de la bulle
                 # dans le repere. La bonne mesure est d(z_lab)/dt avec
                 # z_lab = frame_z + zc, qui ne depend d'aucune convention de
                 # moyenne (cf. scripts/plot_ascension_laminaire.py). Ecart 0.2 %.
FR_FAC = U_LAM / np.sqrt(G*D)      # Fr' = beta * U_LAM/sqrt(gD) = beta * 1.53

# colonnes
CF = {"t": 0, "uz": 1, "z": 2, "ub": 3}                  # frame.dat
CB = {"t": 0, "volume": 1, "chi": 12}                    # bubble.dat

# les 5 points Weber : (label, KE_TARGET, beta = u'/u_inf mesure sur precurseur)
POINTS = [("wt05", 5,  0.15), ("wt11", 11, 0.22), ("wt21", 21, 0.31),
          ("wt25", 25, 0.33), ("wt32", 32, 0.38)]
MEMBERS = [0, 1, 2, 3, 4]


from dataio import load_table


def load(path):
    """Chargeur commun ROBUSTE (scripts/dataio.py).

    L'ancienne version fixait le nombre de colonnes sur la 1re ligne et jetait
    SILENCIEUSEMENT les autres. Depuis la correction d'`event bubble`
    (2026-07-22), bubble.dat a 16 colonnes au lieu de 13 : un run repris avec le
    nouveau binaire produit un fichier a deux formats, dont toute la seconde
    moitie aurait disparu sans un mot. load_table tronque a la largeur commune
    et le signale sur stderr."""
    return load_table(path)


def plateau_uz(frame, t0):
    """Moyenne de frame_uz sur le plateau t>=t0. None si pas de point plateau."""
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
    """Fraction volumique min et chi max sur le plateau (non-fragmentation)."""
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
                    help="debut du plateau frame_uz (defaut 220)")
    args = ap.parse_args()

    base = os.path.join(args.root, "simulations", "weber_ensemble")
    figdir = os.path.join(args.root, "scripts", "figures")
    os.makedirs(figdir, exist_ok=True)

    print(f"\n{'='*66}\n  CAMPAGNE weber_rise -- ensemble (u_inf_lam = {U_LAM}, plateau t>={args.t0:.0f})\n{'='*66}")
    print(f"{'Weber':6} {'beta':>5} {'membres (frame_uz plateau)':32} "
          f"{'moy':>6} {'ratio':>6} {'red%':>6} {'SEM':>6}")

    rows = []
    for lab, ke, beta in POINTS:
        members = []
        for m in MEMBERS:
            d = os.path.join(base, f"{lab}_m{m}")
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
            red = 100.0*(ratio-1.0)   # changement signe ū/u_inf - 1 (<0 = ralentissement)
        else:
            mean = sd = sem = ratio = red = float("nan")
        rows.append(dict(lab=lab, ke=ke, beta=beta, frp=beta*FR_FAC,
                         members=members, n=n, mean=mean, sd=sd, sem=sem,
                         ratio=ratio, red=red))

        mstr = " ".join(f"m{mm['m']}={mm['uz']:.2f}" for mm in members) or "(aucun plateau)"
        sem_str = f"+/-{sem/U_LAM*100:.0f}%" if np.isfinite(sem) else " (n<2)"
        rat_str = f"{ratio:.2f}" if np.isfinite(ratio) else "  -"
        red_str = f"{red:+.0f}" if np.isfinite(red) else "  -"
        mean_str = f"{mean:.2f}" if np.isfinite(mean) else "  -"
        print(f"{lab:6} {beta:5.2f} {mstr:32} {mean_str:>6} {rat_str:>6} "
              f"{red_str:>6} {sem_str:>6}  [n={n}]")

    # --- controle volume / non-fragmentation ---
    print("\n  Controle non-fragmentation (min V/Vth, chi max sur plateau) :")
    for r in rows:
        for mm in r["members"]:
            if mm["vc"]:
                print(f"    {r['lab']} m{mm['m']}: V/Vth_min={mm['vc']['vfrac_min']*100:.0f}%  "
                      f"chi_max={mm['vc']['chi_max']:.2f}  (t_end={mm['tend']:.0f}, n_plat={mm['npts']})")

    # --- CSV ---
    csv = os.path.join(figdir, "weber_ensemble.csv")
    with open(csv, "w", encoding="utf-8") as f:
        f.write("Weber,KE,beta,Fr_prime,n_membres,uz_moy,uz_std,SEM,ratio,reduction_pct,membres_uz\n")
        for r in rows:
            mu = ";".join(f"m{mm['m']}:{mm['uz']:.3f}" for mm in r["members"])
            f.write(f"{r['lab']},{r['ke']},{r['beta']},{r['frp']:.3f},{r['n']},"
                    f"{r['mean']:.3f},{r['sd']:.3f},{r['sem']:.3f},{r['ratio']:.3f},"
                    f"{r['red']:.1f},{mu}\n")
    print(f"\n  CSV   -> {csv}")

    # --- lignes LaTeX pour sec9 (tab:weber_ensemble) ---
    tex = os.path.join(figdir, "weber_ensemble_tab.tex")
    with open(tex, "w", encoding="utf-8") as f:
        f.write("% Genere par analyse_weber_ensemble.py -- coller dans sec9.tex\n")
        f.write("% colonnes : Weber & beta & Fr' & n & u_inf_turb & ratio & reduction \\\\\n")
        for r in rows:
            if r["n"] >= 2:
                val = f"{r['mean']:.2f} $\\pm$ {r['sem']:.2f}"
                red = f"{r['red']:+.0f}\\%"
            elif r["n"] == 1:
                val = f"{r['mean']:.2f} (partiel)"
                red = f"{r['red']:+.0f}\\%"
            else:
                val = "--"; red = "--"
            rat = f"{r['ratio']:.2f}" if np.isfinite(r["ratio"]) else "--"
            f.write(f"        {r['lab']} & {r['beta']:.2f} & {r['frp']:.2f} & "
                    f"{r['n']} & {val} & {rat} & {red} \\\\\n")
    print(f"  LaTeX -> {tex}")

    # --- figure : ratio vs beta, barres d'erreur d'ensemble ---
    fig, ax = plt.subplots(figsize=(8, 6))
    bb = np.linspace(0, 0.42, 100)
    ax.plot(bb, 1.0 - bb**2, "k--", lw=1.2, alpha=0.7,
            label=r"loi petit-$\beta$ : $1-\beta^2$")
    ax.axhline(1.0, color="gray", lw=0.8, ls=":")

    for r in rows:
        if not np.isfinite(r["ratio"]):
            continue
        if r["n"] >= 2:
            ax.errorbar(r["beta"], r["ratio"], yerr=r["sem"]/U_LAM,
                        fmt="o", ms=9, color="C0", capsize=4, zorder=5)
            ax.annotate(f"n={r['n']}", (r["beta"], r["ratio"]),
                        textcoords="offset points", xytext=(8, 6), fontsize=8)
        else:
            ax.plot(r["beta"], r["ratio"], "o", ms=9, mfc="white",
                    mec="C1", mew=1.8, zorder=5)
            ax.annotate(f"n=1\n(partiel)", (r["beta"], r["ratio"]),
                        textcoords="offset points", xytext=(8, -18),
                        fontsize=8, color="C1")

    # legende propre (marqueurs factices)
    ax.plot([], [], "o", color="C0", label="ensemble ($n\\geq2$, barre = SEM)")
    ax.plot([], [], "o", mfc="white", mec="C1", mew=1.8, label="partiel ($n=1$)")

    ax.set_xlabel(r"$\beta = u'/u_\infty$")
    ax.set_ylabel(r"$\bar u_{\infty,\mathrm{turb}} / u_\infty$")
    ax.set_title(r"Effet de la turbulence sur la vitesse de montee ($Bo=1$, $Ga\simeq70$)")
    ax.set_xlim(0, 0.42)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="lower left")
    # axe secondaire Fr'
    axb = ax.secondary_xaxis("top", functions=(lambda b: b*FR_FAC,
                                               lambda fr: fr/FR_FAC))
    axb.set_xlabel(r"$Fr' = u'/\sqrt{gD}$")
    plt.tight_layout()
    out = os.path.join(figdir, "weber_rise_prelim.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  figure-> {out}")
    n_full = sum(1 for r in rows if r["n"] >= 2)
    print(f"\n  {n_full}/5 points a >=2 membres ; "
          f"{sum(1 for r in rows if r['n']>=1)}/5 points avec au moins 1 membre.\n")


if __name__ == "__main__":
    main()
