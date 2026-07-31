#!/usr/bin/env python3
"""
check_stationnarite.py -- Vérification du régime permanent (vitesse terminale).
"""
import os, numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V_THEO = 4.0/3.0*np.pi*8.0**3
U_LAM  = 12.27

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
        return
    t=fr[:,0]; uz=fr[:,1]; ub=fr[:,3]
    tend=t[-1]
    m=t>=t0
    if m.sum()<5:
        print(f"  {tag:14} : t_end={tend:.0f}  NON CONVERGE (run court)")
        return
    tt=t[m]; uu=uz[m]; bb=ub[m]
    A=np.polyfit(tt,uu,1); slope=A[0]; mean=uu.mean()
    slope100=slope*100.0
    pct100=slope100/mean*100.0
    half=tt[0]+(tt[-1]-tt[0])/2
    m1=uu[tt<half].mean(); m2=uu[tt>=half].mean()
    drift_halves=(m2-m1)/mean*100.0
    ubm=bb.mean(); ubstd=bb.std()
    volinfo=""
    bub=load(bub_p)
    if bub is not None:
        tb=bub[:,0]; vol=bub[:,1]; mb=tb>=t0
        if mb.sum()>=5:
            vr=np.polyfit(tb[mb],vol[mb],1)[0]*100.0/V_THEO*100.0
            v0=vol[mb][0]/V_THEO*100; v1=vol[mb][-1]/V_THEO*100
            volinfo=f"vol {v0:.0f}->{v1:.0f}%Vth"
    ok_slope=abs(pct100)<5.0
    ok_halves=abs(drift_halves)<4.0
    ok_ub=abs(ubm)<0.15*abs(mean) if mean else False
    verdict = "PERMANENT" if (ok_slope and ok_halves and ok_ub) else \
              ("QUASI" if (abs(pct100)<10 and abs(drift_halves)<8) else "NON-CONVERGE")
    print(f"  {tag:14} u_inf={mean:6.2f} | pente={slope100:+6.2f}/100ut | ub={ubm:+.2f} | {volinfo} => {verdict}")

print("--- Régime permanent (LAMINAIRE) ---")
analyse("laminaire", f"{ROOT}/simulations/lam7_frame/frame.dat",
        f"{ROOT}/simulations/lam7_frame/bubble.dat", 210)

print("\n--- Régime permanent (TURBULENCE - Betasweep) ---")
for lab in ["beta0150","beta0220","beta0310","beta0330","beta0380"]:
    for mm in range(5):
        d=f"{ROOT}/simulations/betasweep/{lab}_m{mm}_bub"
        if os.path.isdir(d):
            analyse(f"{lab}_m{mm}", f"{d}/frame.dat", f"{d}/bubble.dat", 220)
