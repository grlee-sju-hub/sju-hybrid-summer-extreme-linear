#!/usr/bin/env python3
"""Hybrid forecast of JA CDHD from the two concurrent (L0) track-A predictors:
   t2mna = JA 2-m temperature, tropical North Atlantic;  v10io = JA 10-m meridional wind, eastern Indian Ocean.
Both predictors are concurrent (L0), so at forecast time they are supplied by the C3S multi-model-ensemble
(MME) forecast, calibrated to the observed distribution. Two protocols, evaluated over the C3S hindcast
overlap 1993-2016:
  A) 1993-train : fit & leave-one-out within 1993-2016.
  B) 1979-train : fit beta on observed 1979-2024 (LOO 45-yr fits); predict 1993-2016 from the calibrated C3S MME.
Out: console + data/cdhd_trackA_hybrid.csv."""
import numpy as np, pandas as pd, xarray as xr
from pathlib import Path
D=Path("data")
def detr(x):
    x=np.asarray(x,float); n=len(x); t=np.arange(n)-(n-1)/2.; a=x-x.mean(); s=(t*t).sum(); return a-((t*a).sum()/s)*t
def corr(a,b): return float(np.corrcoef(a,b)[0,1])
def pploo(Xf,Xp,y):
    XtXi=np.linalg.inv(Xf.T@Xf); beta=XtXi@(Xf.T@y); e=y-Xf@beta
    h=np.einsum("ij,jk,ik->i",Xf,XtXi,Xf); cross=np.einsum("ij,jk,ik->i",Xp,XtXi,Xf)
    yloo=Xp@beta-cross*e/(1-h); return corr(yloo,y)
Y46=np.arange(1979,2025); E93=np.arange(1993,2017); pos={int(Y):k for k,Y in enumerate(Y46)}; PREDN=["t2mna","v10io"]

# predictand: CDHD (Korea-mean series)
d=xr.open_dataset(D/"ERA5_JA_CDHD_index_1979-2024.nc"); yc=d.year.values.astype(int); cd=d.korea_nday.values.astype(float); d.close()
ym={int(a):float(b) for a,b in zip(yc,cd)}
# observed track-A predictor indices (JA box means)
ta=pd.read_csv(D/"era5_JA_trackA_predictors_1979-2024.csv").set_index("year")
OBS={p:ta.reindex(Y46)[p].values.astype(float) for p in PREDN}
# C3S MME forecast of the L0 predictors (year -> value)
mme=pd.read_csv(D/"c3s_mme_L0_ja.csv").set_index("year")
MME={p:{int(Y):float(mme.loc[Y,p]) for Y in mme.index} for p in PREDN}

# ---- A) 1993-train : fit & LOO within 1993-2016 ----
i93=[pos[Y] for Y in E93]; yA=detr(np.array([ym[Y] for Y in E93]))
Xo=np.column_stack([np.ones(len(E93))]+[detr(OBS[p][i93]) for p in PREDN])
Xh=np.column_stack([np.ones(len(E93))]+[detr(np.array([MME[p][Y] for Y in E93])) for p in PREDN])
accA_obs=pploo(Xo,Xo,yA); accA_hyb=pploo(Xo,Xh,yA)

# ---- B) 1979-train, predict 1993-2016 (explicit LOO; C3S calibrated to observed-L0 distribution) ----
y46=detr(np.array([ym[Y] for Y in Y46])); Xf46={p:detr(OBS[p]) for p in PREDN}
cal={}
for p in PREDN:
    oL0=np.array([Xf46[p][pos[Y]] for Y in E93]); c=np.array([MME[p][Y] for Y in E93])
    cc=(c-c.mean())/c.std()*oL0.std()+oL0.mean(); cal[p]={int(Y):cc[k] for k,Y in enumerate(E93)}
yo=[]; yh=[]; yt=[]
for Y in E93:
    tr=Y46!=Y; Xtr=np.column_stack([np.ones(tr.sum())]+[Xf46[p][tr] for p in PREDN]); beta=np.linalg.lstsq(Xtr,y46[tr],rcond=None)[0]
    iY=pos[Y]; xo=np.array([1.0]+[Xf46[p][iY] for p in PREDN]); xh=np.array([1.0]+[cal[p][Y] for p in PREDN])
    yo.append(xo@beta); yh.append(xh@beta); yt.append(y46[iY])
accB_obs=corr(yo,yt); accB_hyb=corr(yh,yt)

res=pd.DataFrame([dict(index="CDHD", predictors="t2mna+v10io (track A, L0)",
    A_obs93=round(accA_obs,3), A_hyb93_C3S=round(accA_hyb,3),
    B_obs_tr79=round(accB_obs,3), B_hyb_tr79_C3S=round(accB_hyb,3))])
print(res.to_string(index=False),flush=True); res.to_csv(D/"cdhd_trackA_hybrid.csv",index=False)
print("[saved] cdhd_trackA_hybrid.csv  (CDHD track-A hybrid, eval 1993-2016)",flush=True)
