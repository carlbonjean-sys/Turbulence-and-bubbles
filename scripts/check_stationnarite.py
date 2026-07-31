#!/usr/bin/env python3
# ============================================================================
#  check_stationnarite.py -- Atteint-on un REGIME PERMANENT (vitesse terminale) ?
#  Teste, sur frame.dat (u_inf = frame_uz du repere mobile) :
#    (1) pente residuelle du plateau (regression lineaire)  -> ~0 attendu
#    (2) stationnarite : moyenne 1re moitie vs 2e moitie du plateau
#    (3) ub = vitesse de la bulle DANS le repere -> ~0 attendu a l'equilibre
#    (4) derive de volume (bubble.dat) sur la fenetre -> diametre constant ?
#  Verdict par run : PERMANENT / QUASI / NON-CONVERGE.
# ============================================================================
import os, numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V_THEO = 4.0/3.0*np.pi*8.0**3
U_LAM  = 12.24

def load(p):
    if not os.path.isfile(p): return None
    r=[]
    for ln in open(p):
        s=ln.strip()
        if not s or s.startswith('#'): continue
        try: v=[float(x) for x in s.split()]
        except: continue
        r.append(v)
    return np.array(r) if r else None

def analyse(tag, frame_p, bub_p, t0):
    fr=load(frame_p)
    if fr is None:
        print(f"  {tag:14} : pas de frame.dat"); return
    t=fr[:,0]; uz=fr[:,1]; ub=fr[:,3]
    tend=t[-1]
    m=t>=t0
    if m.sum()<5:
        print(f"  {tag:14} : t_end={tend:.0f}  plateau t>={t0:.0f} vide (n={m.sum()}) -> NON CONVERGE (run court)")
        return
    tt=t[m]; uu=uz[m]; bb=ub[m]
    # (1) pente residuelle
    A=np.polyfit(tt,uu,1); slope=A[0]; mean=uu.mean()
    slope100=slope*100.0                 # variation sur 100 u.t.
    pct100=slope100/mean*100.0           # en % de la moyenne
    # (2) stationnarite : 1re vs 2e moitie
    half=tt[0]+(tt[-1]-tt[0])/2
    m1=uu[tt<half].mean(); m2=uu[tt>=half].mean()
    drift_halves=(m2-m1)/mean*100.0
    # (3) ub moyen (bulle dans repere) -> ~0
    ubm=bb.mean(); ubstd=bb.std()
    # (4) volume
    volinfo=""
    bub=load(bub_p)
    if bub is not None:
        tb=bub[:,0]; vol=bub[:,1]; mb=tb>=t0
        if mb.sum()>=5:
            vr=np.polyfit(tb[mb],vol[mb],1)[0]*100.0/V_THEO*100.0  # %Vth /100ut
            v0=vol[mb][0]/V_THEO*100; v1=vol[mb][-1]/V_THEO*100
            volinfo=f"vol {v0:.0f}->{v1:.0f}%Vth ({vr:+.1f}%/100ut)"
    # verdict
    ok_slope=abs(pct100)<5.0            # <5% de derive sur 100 u.t.
    ok_halves=abs(drift_halves)<4.0
    ok_ub=abs(ubm)<0.15*abs(mean) if mean else False
    verdict = "PERMANENT" if (ok_slope and ok_halves and ok_ub) else \
              ("QUASI" if (abs(pct100)<10 and abs(drift_halves)<8) else "NON-CONVERGE")
    print(f"  {tag:14} t[{t0:.0f},{tend:.0f}] n={m.sum():4d} | u_inf={mean:6.2f}"
          f" | pente={slope100:+6.2f}/100ut ({pct100:+5.1f}%) | 1/2v2/2={drift_halves:+5.1f}%"
          f" | ub={ubm:+.2f}±{ubstd:.2f} | {volinfo}  => {verdict}")
    return dict(tag=tag,mean=mean,pct100=pct100,drift=drift_halves,ub=ubm,verdict=verdict)

print("="*120)
print("  REGIME PERMANENT / VITESSE TERMINALE  (u_inf = plateau frame_uz ; criteres: |pente|<5%/100ut, |1/2-2/2|<4%, |ub|<0.15 u_inf)")
print("="*120)
print("\n  --- LAMINAIRE de reference (beta=0) ---")
analyse("laminaire", f"{ROOT}/simulations/lam7_frame/frame.dat",
        f"{ROOT}/simulations/lam7_frame/bubble.dat", 210)

print("\n  --- WEBER (turbulents), par membre ---")
for lab in ["wt05","wt11","wt21","wt25","wt32"]:
    for mm in [0,1,2]:
        d=f"{ROOT}/simulations/weber_ensemble/{lab}_m{mm}"
        if os.path.isdir(d):
            analyse(f"{lab}_m{mm}", f"{d}/frame.dat", f"{d}/bubble.dat", 220)

# momentum laminaire
print("\n  --- Garde-fou momentum (laminaire) : invariant pz_lab = pz+frame_uz ~ 0 ---")
mom=load(f"{ROOT}/simulations/lam7_frame/momentum.dat")
if mom is not None:
    # colonnes: t px py pz frame_uz pz_lab
    t=mom[:,0]; px=mom[:,1]; py=mom[:,2]; pz_lab=mom[:,5]
    m=t>=210
    print(f"    sur t>=210 : |px|max={np.abs(px[m]).max():.2e}  |py|max={np.abs(py[m]).max():.2e}"
          f"  |pz_lab|max={np.abs(pz_lab[m]).max():.2e}  (={np.abs(pz_lab[m]).max()/U_LAM*100:.2f}% de u_inf)")
    print(f"    derive pz_lab: {np.polyfit(t[m],pz_lab[m],1)[0]*100:+.2e} /100ut")
else:
    print("    momentum.dat absent")
