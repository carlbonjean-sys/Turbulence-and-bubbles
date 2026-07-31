#!/usr/bin/env python3
# ============================================================================
#  postprocess.py -- Courbes de convergence spatiale + statistiques turbulence
# ----------------------------------------------------------------------------
#  Lit :
#    simulations/precursor/stats.dat      (# t dissipation energy Reynolds)
#    simulations/lvl{7,8,9}/bubble.dat    (# t volume xc yc zc ux uy uz
#                                             width depth height d_eq chi)
#    simulations/lvl{7,8,9}/stats.dat     (turbulence PENDANT le run bulle,
#                                          pour tau = k/eps)
#
#  Produit (dans scripts/figures/) :
#    - rise_velocity.png       : vitesse de remontee u_z(t) pour chaque niveau
#    - cumulative_mean.png     : moyenne cumulee U_cum(t;t0) vs depuis injection (P6)
#    - t0_robustness.png       : U_b(t0) -> plateau = choix de t0 non arbitraire (P6)
#    - convergence.png         : U_b (regression depuis t0) vs maxlevel
#    - divergence_lvl7_lvl8.png: |z7-z8|(t), |uz7-uz8|(t) semilog (P7, chaos)
#    - aspect_ratio.png        : chi(t) et width/depth/height
#    - trajectoire.png         : montee + derive horizontale (deperiodisees)
#    - turbulence.png          : energie, dissipation, Reynolds du precurseur
#    - convergence_table.csv   : recap chiffre (U_b, C_D, Re_b, n_tau, traversees)
#
#  P6 (tuteurs) : la moyenne demarre a t0 = injection + transitoire, PAS a
#  l'injection. U_b = pente de la regression lineaire de z_deperiodise(t) sur
#  [t0, fin] (robuste), plus la moyenne cumulee pour visualiser la convergence.
#
#  Usage :
#    python scripts/postprocess.py [--root .] [--levels 7 8 9] [--tail 0.2]
#                                  [--t0-offset 8.0] [--t0 ABS]
# ============================================================================
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Colonnes de bubble.dat
COL = {"t": 0, "volume": 1, "xc": 2, "yc": 3, "zc": 4,
       "ux": 5, "uy": 6, "uz": 7, "width": 8, "depth": 9,
       "height": 10, "d_eq": 11, "chi": 12}

L0 = 120.0                                  # WIDTH (code units)
G = 1.0                                     # |G_Z|
DRHO_RHO = 1.0 - 1.0 / 850.0                # (rho1-rho2)/rho1
NU = 0.01 * (L0 / (2.0 * np.pi)) ** 2 / 2.0  # mu1/rho1 du solveur (~1.824)


from dataio import load_table


def load(path):
    """Chargeur commun ROBUSTE (scripts/dataio.py).

    ATTENTION, l'ancienne version se disait "robuste" mais fixait le nombre de
    colonnes sur la 1re ligne de donnees et jetait SILENCIEUSEMENT toutes les
    autres largeurs. Depuis la correction d'`event bubble` (2026-07-22),
    bubble.dat a 16 colonnes au lieu de 13 : un run repris avec le nouveau
    binaire aurait vu toute sa seconde moitie ignoree sans un mot. load_table
    tronque a la largeur commune (les colonnes 1-13 ont garde leur sens) et
    signale sur stderr des qu'il tronque ou ecarte quoi que ce soit."""
    if not os.path.isfile(path):
        print(f"  [skip] {path} introuvable")
        return None
    return load_table(path)


def tail_mean(t, y, frac):
    """Moyenne sur la derniere fraction temporelle (etat quasi-stationnaire)."""
    if len(t) == 0:
        return np.nan, np.nan
    t0 = t[-1] - frac * (t[-1] - t[0])
    mask = t >= t0
    return np.mean(y[mask]), np.std(y[mask])


def unwrap_periodic(p, period=L0):
    """Deperiodise une trajectoire (image minimale entre points successifs)."""
    d = np.diff(p)
    d -= period * np.round(d / period)
    return np.concatenate([[p[0]], p[0] + np.cumsum(d)])


def ub_regression(t, z, t0):
    """U_b = pente de la regression lineaire de z(t) sur [t0, fin] (P6)."""
    mask = t >= t0
    if mask.sum() < 10:
        return np.nan
    return np.polyfit(t[mask], z[mask], 1)[0]


def tau_turnover(stats, t0=None):
    """Temps de retournement tau = <k/eps> sur le run (stats.dat du niveau)."""
    if stats is None:
        return np.nan
    t, eps, k = stats[:, 0], stats[:, 1], stats[:, 2]
    mask = np.isfinite(eps) & (eps > 0)
    if t0 is not None:
        mask &= t >= t0
    if mask.sum() == 0:
        return np.nan
    return float(np.mean(k[mask] / eps[mask]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--levels", type=int, nargs="+", default=[7, 8, 9])
    ap.add_argument("--tail", type=float, default=0.2,
                    help="fraction finale pour les moyennes stationnaires")
    ap.add_argument("--t0-offset", type=float, default=8.0,
                    help="debut de moyenne = injection + offset (defaut 8 = "
                         "2*sqrt(D/g), transitoire de mise en mouvement)")
    ap.add_argument("--t0", type=float, default=None,
                    help="debut de moyenne ABSOLU (prioritaire sur --t0-offset)")
    args = ap.parse_args()

    figdir = os.path.join(args.root, "scripts", "figures")
    os.makedirs(figdir, exist_ok=True)

    bub, stt = {}, {}
    for lvl in args.levels:
        d = load(os.path.join(args.root, "simulations", f"lvl{lvl}", "bubble.dat"))
        if d is not None:
            bub[lvl] = d
        stt[lvl] = load(os.path.join(args.root, "simulations", f"lvl{lvl}",
                                     "stats.dat"))
    levels_sorted = sorted(bub.keys())

    # t0 : injection + transitoire (P6). Injection = 1er temps de bubble.dat.
    t0_abs = {}
    for lvl in levels_sorted:
        t_inj = bub[lvl][0, COL["t"]]
        t0_abs[lvl] = args.t0 if args.t0 is not None else t_inj + args.t0_offset

    # ---------------- 1) Vitesse de remontee u_z(t) ----------------
    plt.figure(figsize=(7, 5))
    for lvl, d in bub.items():
        plt.plot(d[:, COL["t"]], d[:, COL["uz"]], label=f"maxlevel {lvl}")
    if levels_sorted:
        plt.axvline(t0_abs[levels_sorted[0]], color="k", ls=":", lw=1,
                    label=r"$t_0$ (debut moyenne)")
    plt.xlabel("temps t")
    plt.ylabel(r"vitesse de remontee  $u_z$")
    plt.title("Vitesse de remontee de la bulle")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "rise_velocity.png"), dpi=150)
    plt.close()

    # ---------------- 2) P6 : moyenne cumulee U_cum(t; t0) ----------------
    #     U_cum(t;t0) = [z(t) - z(t0)] / (t - t0), z deperiodise.
    plt.figure(figsize=(7, 5))
    for lvl in levels_sorted:
        d = bub[lvl]
        t = d[:, COL["t"]]
        z = unwrap_periodic(d[:, COL["zc"]])
        t0 = t0_abs[lvl]
        # depuis l'injection (ancienne pratique, en pointille pour comparaison)
        m_inj = t > t[0] + 1e-9
        plt.plot(t[m_inj], (z[m_inj] - z[0]) / (t[m_inj] - t[0]),
                 ls=":", alpha=0.6, label=f"lvl{lvl} depuis injection")
        # depuis t0 (P6)
        i0 = np.searchsorted(t, t0)
        if i0 < len(t) - 10:
            tt, zz = t[i0:], z[i0:]
            m = tt > tt[0] + 1e-9
            plt.plot(tt[m], (zz[m] - zz[0]) / (tt[m] - tt[0]),
                     lw=2, label=f"lvl{lvl} depuis $t_0$")
    plt.xlabel("temps t")
    plt.ylabel(r"$U_{cum}(t) = [z(t)-z(t_0)]/(t-t_0)$")
    plt.title("Moyenne cumulee de la vitesse de remontee (P6)")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "cumulative_mean.png"), dpi=150)
    plt.close()

    # ---------------- 3) P6 : robustesse du choix de t0 ----------------
    plt.figure(figsize=(7, 5))
    for lvl in levels_sorted:
        d = bub[lvl]
        t = d[:, COL["t"]]
        z = unwrap_periodic(d[:, COL["zc"]])
        t0s = np.arange(t[0], t[0] + min(25.0, 0.5 * (t[-1] - t[0])), 0.5)
        ubs = [ub_regression(t, z, t0) for t0 in t0s]
        plt.plot(t0s - t[0], ubs, marker=".", ms=3, label=f"maxlevel {lvl}")
    plt.axvline(args.t0_offset, color="k", ls=":", lw=1, label="offset retenu")
    plt.xlabel(r"decalage $t_0 - t_{inj}$")
    plt.ylabel(r"$U_b$ (pente regression sur $[t_0,$ fin$]$)")
    plt.title(r"Robustesse de $U_b$ au choix de $t_0$ (P6)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "t0_robustness.png"), dpi=150)
    plt.close()

    # ---------------- 4) Convergence + table (U_b regression, C_D, Re_b,
    #                     duree en tau et traversees) ----------------
    rows = []
    for lvl in levels_sorted:
        d = bub[lvl]
        t = d[:, COL["t"]]
        z = unwrap_periodic(d[:, COL["zc"]])
        t0 = t0_abs[lvl]
        ub_reg = ub_regression(t, z, t0)
        ub_full = (z[-1] - z[0]) / (t[-1] - t[0])          # ancienne definition
        uz_tail, uz_std = tail_mean(t, d[:, COL["uz"]], args.tail)
        chi_m, _ = tail_mean(t, d[:, COL["chi"]], args.tail)
        deq_m, _ = tail_mean(t, d[:, COL["d_eq"]], args.tail)
        tau = tau_turnover(stt.get(lvl))
        T_win = t[-1] - t[0]
        n_tau = T_win / tau if np.isfinite(tau) and tau > 0 else np.nan
        i0 = np.searchsorted(t, t0)
        n_cross = (z[-1] - z[min(i0, len(z) - 1)]) / L0
        cd = (4.0 * G * deq_m * DRHO_RHO / (3.0 * ub_reg ** 2)
              if np.isfinite(ub_reg) and abs(ub_reg) > 1e-9 else np.nan)
        re_b = (ub_reg * deq_m / NU if np.isfinite(ub_reg) else np.nan)
        rows.append(dict(lvl=lvl, ub_reg=ub_reg, ub_full=ub_full,
                         uz_tail=uz_tail, uz_std=uz_std, chi=chi_m, deq=deq_m,
                         cd=cd, re_b=re_b, tau=tau, n_tau=n_tau,
                         n_cross=n_cross, T_win=T_win, t0=t0))

    if rows:
        plt.figure(figsize=(7, 5))
        plt.plot([r["lvl"] for r in rows], [r["ub_reg"] for r in rows],
                 marker="o", lw=1.5, label=r"$U_b$ regression depuis $t_0$ (P6)")
        plt.plot([r["lvl"] for r in rows], [r["ub_full"] for r in rows],
                 marker="s", ls="--", lw=1.2,
                 label=r"$U_b$ fenetre complete (ancien)")
        plt.xlabel("maxlevel")
        plt.ylabel(r"$U_b$")
        plt.title("Convergence spatiale de la vitesse de remontee")
        plt.xticks([r["lvl"] for r in rows])
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, "convergence.png"), dpi=150)
        plt.close()

        hdr = (" lvl |  U_b(reg,t0) | U_b(fenetre) | uz_tail |  chi  |  d_eq |"
               "   C_D  |  Re_b | tau  | n_tau | traversees")
        print("\n" + hdr)
        with open(os.path.join(figdir, "convergence_table.csv"), "w") as f:
            f.write("maxlevel,t0,Ub_regression,Ub_fenetre_complete,uz_tail,"
                    "uz_std,chi,d_eq,C_D,Re_b,tau,n_turnovers,n_traversees,"
                    "duree_fenetre\n")
            for r in rows:
                print(f"  {r['lvl']}  | {r['ub_reg']:12.5f} |"
                      f" {r['ub_full']:12.5f} | {r['uz_tail']:7.3f} |"
                      f" {r['chi']:.3f} | {r['deq']:.3f} | {r['cd']:6.2f} |"
                      f" {r['re_b']:5.1f} | {r['tau']:4.1f} |"
                      f" {r['n_tau']:5.2f} | {r['n_cross']:6.3f}")
                f.write(f"{r['lvl']},{r['t0']},{r['ub_reg']},{r['ub_full']},"
                        f"{r['uz_tail']},{r['uz_std']},{r['chi']},{r['deq']},"
                        f"{r['cd']},{r['re_b']},{r['tau']},{r['n_tau']},"
                        f"{r['n_cross']},{r['T_win']}\n")

    # ---------------- 5) P7 : divergence chaotique lvl7 vs lvl8 ----------------
    if 7 in bub and 8 in bub:
        d7, d8 = bub[7], bub[8]
        t7, t8 = d7[:, COL["t"]], d8[:, COL["t"]]
        z7 = unwrap_periodic(d7[:, COL["zc"]])
        z8 = unwrap_periodic(d8[:, COL["zc"]])
        ta = max(t7[0], t8[0])
        tb = min(t7[-1], t8[-1])
        tc = np.linspace(ta, tb, 2000)
        dz = np.abs(np.interp(tc, t7, z7) - np.interp(tc, t8, z8))
        duz = np.abs(np.interp(tc, t7, d7[:, COL["uz"]])
                     - np.interp(tc, t8, d8[:, COL["uz"]]))
        deq = np.nanmean(d8[:, COL["d_eq"]])

        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        ax[0].semilogy(tc - ta, np.maximum(dz, 1e-6))
        ax[0].axhline(deq, color="r", ls="--", lw=1,
                      label=f"diametre bulle ({deq:.1f})")
        ax[0].set_xlabel("temps depuis injection")
        ax[0].set_ylabel(r"$|z_7(t) - z_8(t)|$")
        ax[0].set_title("Divergence des trajectoires (position)")
        ax[0].grid(True, alpha=0.3, which="both")
        ax[0].legend()
        ax[1].semilogy(tc - ta, np.maximum(duz, 1e-6))
        ax[1].set_xlabel("temps depuis injection")
        ax[1].set_ylabel(r"$|u_{z,7}(t) - u_{z,8}(t)|$")
        ax[1].set_title("Divergence des trajectoires (vitesse)")
        ax[1].grid(True, alpha=0.3, which="both")
        fig.suptitle("P7 : meme precurseur, maillages differents -> les "
                     "trajectoires decorrelent (chaos), les statistiques "
                     "convergent", fontsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, "divergence_lvl7_lvl8.png"), dpi=150)
        plt.close()

    # ---------------- 6) Rapport d'aspect ----------------
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for lvl, d in bub.items():
        ax[0].plot(d[:, COL["t"]], d[:, COL["chi"]], label=f"maxlevel {lvl}")
    ax[0].set_xlabel("temps t")
    ax[0].set_ylabel(r"rapport d'aspect $\chi$")
    ax[0].set_title("Aspect ratio de la bulle")
    ax[0].grid(True, alpha=0.3)
    ax[0].legend()

    finest = bub.get(levels_sorted[-1]) if levels_sorted else None
    if finest is not None:
        ax[1].plot(finest[:, COL["t"]], finest[:, COL["width"]], label="width (x)")
        ax[1].plot(finest[:, COL["t"]], finest[:, COL["depth"]], label="depth (y)")
        ax[1].plot(finest[:, COL["t"]], finest[:, COL["height"]], label="height (z)")
        ax[1].set_xlabel("temps t")
        ax[1].set_ylabel("dimensions de la bulle")
        ax[1].set_title(f"Dimensions (maxlevel {levels_sorted[-1]})")
        ax[1].grid(True, alpha=0.3)
        ax[1].legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "aspect_ratio.png"), dpi=150)
    plt.close()

    # ---------------- 7) Trajectoire de la bulle (montee + derive) ----------------
    if finest is not None:
        t = finest[:, COL["t"]]
        z = unwrap_periodic(finest[:, COL["zc"]])
        fig = plt.figure(figsize=(12, 5))
        a1 = fig.add_subplot(1, 2, 1)
        a1.plot(t, z - z[0])
        a1.set_xlabel("temps t"); a1.set_ylabel("montee z(t) - z(0)")
        a1.set_title("Montee de la bulle"); a1.grid(True, alpha=0.3)
        a2 = fig.add_subplot(1, 2, 2)
        a2.plot(finest[:, COL["xc"]], finest[:, COL["yc"]], lw=1)
        a2.scatter([finest[0, COL["xc"]]], [finest[0, COL["yc"]]],
                   c="g", zorder=5, label="depart")
        a2.set_xlabel("x"); a2.set_ylabel("y"); a2.axis("equal")
        a2.set_title("Derive horizontale"); a2.grid(True, alpha=0.3); a2.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, "trajectoire.png"), dpi=150)
        plt.close()

    # ---------------- 8) Statistiques de turbulence (precurseur) ----------------
    st = load(os.path.join(args.root, "simulations", "precursor", "stats.dat"))
    if st is not None:
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        ax[0].plot(st[:, 0], st[:, 2]); ax[0].set_title("Energie cinetique k")
        ax[1].plot(st[:, 0], st[:, 1]); ax[1].set_title("Dissipation epsilon")
        ax[2].plot(st[:, 0], st[:, 3]); ax[2].set_title(r"$Re_\lambda$")
        for a in ax:
            a.set_xlabel("temps t")
            a.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, "turbulence.png"), dpi=150)
        plt.close()

    print(f"\n  Figures ecrites dans : {figdir}")


if __name__ == "__main__":
    main()
