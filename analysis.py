import rasterio, numpy as np, json, warnings
from scipy import ndimage as ndi
warnings.filterwarnings("ignore")
D="data/"
dw=rasterio.open(D+"bujumbura_dw_built_and_counts.tif").read().astype(float)
b=dw[:10]; n=dw[10:]
wc20,wc21,el,sl=rasterio.open(D+"bujumbura_landcover_and_terrain.tif").read().astype(float)
gsrc=rasterio.open(D+"bujumbura_ghsl_worldpop_100m.tif"); g=gsrc.read().astype(float)
land=(wc21!=80)&(wc20!=80)
PX=4e-4  # km2 per 20 m px
R={}
def comp(idx):
    w=(n[idx]>=3).astype(float); return (b[idx]*w).sum(0)/np.maximum(w.sum(0),1)
E=comp([0,1,2]); L=comp([7,8,9]); M=comp([3,4,5])
R["pixels"]=dict(land_km2=float(land.sum()*PX))
# 1 calibration vs WC2020
cal={}
for t in (30,40,50,60,70,80):
    p=(M>=t)&land; r=(wc20==50)&land; tp=(p&r).sum()
    cal[t]=dict(km2=float(p.sum()*PX),prec=float(tp/p.sum()),rec=float(tp/r.sum()),f1=float(2*tp/(p.sum()+r.sum())))
R["calib"]=cal
# 2 yearly series
ser={}
for t in (50,70):
    ser[t]=[float((((b[i]>=t)&land&(n[i]>=3)).sum())*PX) for i in range(10)]
R["yearly"]=ser; R["nobs_median"]=[float(np.median(n[i])) for i in range(10)]
# 3 gain/loss at tau
out={}
for t in (50,70):
    Eb=(E>=t)&land; Lb=(L>=t)&land
    gain=Lb&~Eb; loss=Eb&~Lb
    # placebo: 2023 vs 2024-25 composite
    P1=comp([7]); P2=comp([8,9]); pb=(P2>=t)&~(P1>=t)&land; pl=(P1>=t)&~(P2>=t)&land
    out[t]=dict(early=float(Eb.sum()*PX),late=float(Lb.sum()*PX),gain=float(gain.sum()*PX),loss=float(loss.sum()*PX),
                placebo_gain=float(pb.sum()*PX),placebo_loss=float(pl.sum()*PX))
R["change"]=out
# main (tau 70) masks
t=70; Eb=(E>=t)&land; Lb=(L>=t)&land; gain=Lb&~Eb; loss=Eb&~Lb
# robust gain: also built in WorldCover 2021 (independent product)
gain_wc=gain&(wc21==50)
R["gain_confirmed_wc21_km2"]=float(gain_wc.sum()*PX); R["gain_share_wc21"]=float(gain_wc.sum()/gain.sum())
# components: size of gain patches
lab,nl=ndi.label(gain); sz=np.bincount(lab.ravel())[1:]*PX
R["patches"]=dict(n=int(nl),n_ge_1ha=int((sz>=0.01).sum()),area_ge_1ha=float(sz[sz>=0.01].sum()),largest=float(sz.max()))
# 4 terrain
cand=land&~Eb   # places where growth was possible
def rate(mask_bins,labels,var):
    rows=[]
    for (lo,hi),lab_ in zip(mask_bins,labels):
        m=cand&(var>=lo)&(var<hi)
        rows.append(dict(bin=lab_,cand_km2=float(m.sum()*PX),gain_km2=float((gain&m).sum()*PX),rate=float((gain&m).sum()/max(m.sum(),1)*100)))
    return rows
R["slope"]=rate([(0,5),(5,15),(15,25),(25,90)],["<5","5-15","15-25",">25"],sl)
R["elev"]=rate([(0,790),(790,900),(900,1200),(1200,1600),(1600,4000)],["<790","790-900","900-1200","1200-1600",">1600"],el)
# distance to early built
dist=ndi.distance_transform_edt(~Eb)*20.0
R["dist"]=rate([(0,100),(100,300),(300,1000),(1000,1e9)],["<100 m","100-300 m","300-1000 m",">1000 m"],dist)
# infill/edge/leapfrog share of gain area
R["gain_by_dist_share"]={k:float((gain&(dist>=lo)&(dist<hi)).sum()/gain.sum()) for k,(lo,hi) in {"<100":(0,100),"100-300":(100,300),"300-1000":(300,1000),">1000":(1000,1e9)}.items()}
# low-lying plain (<790 m, ~17 m above lake) and steep >15 deg exposure of ALL late built and of gain
low=el<790; steep=sl>=15
R["exposure"]=dict(late_built_low_pct=float((Lb&low).sum()/Lb.sum()*100),early_built_low_pct=float((Eb&low).sum()/Eb.sum()*100),
  gain_low_pct=float((gain&low).sum()/gain.sum()*100),late_built_steep_pct=float((Lb&steep).sum()/Lb.sum()*100),
  early_built_steep_pct=float((Eb&steep).sum()/Eb.sum()*100),gain_steep_pct=float((gain&steep).sum()/gain.sum()*100),
  gain_steep_km2=float((gain&steep).sum()*PX),gain_low_km2=float((gain&low).sum()*PX),
  land_low_pct=float(low[land].mean()*100),land_steep_pct=float(steep[land].mean()*100))
# Spatial-block bootstrap (1 km blocks) for gain share on steep/low
H,W=gain.shape; bs=50
by=np.arange(H)//bs; bx=np.arange(W)//bs; bid=(by[:,None]*(W//bs+2)+bx[None,:])
nb=bid.max()+1
def agg(m): return np.bincount(bid.ravel(),weights=m.ravel().astype(float),minlength=nb)
G=agg(gain); GS=agg(gain&steep); GL=agg(gain&low); LB=agg(Lb); LBS=agg(Lb&steep); LBL=agg(Lb&low); EB=agg(Eb); EBS=agg(Eb&steep); EBL=agg(Eb&low)
rng=np.random.default_rng(1); ok=np.where((G+LB+EB)>0)[0]; bt=[]
for _ in range(2000):
    s=rng.choice(ok,len(ok)); 
    bt.append([GS[s].sum()/G[s].sum()*100, EBS[s].sum()/EB[s].sum()*100, GL[s].sum()/G[s].sum()*100, EBL[s].sum()/EB[s].sum()*100, G[s].sum()*PX])
bt=np.array(bt); R["boot"]={k:[float(np.percentile(bt[:,i],2.5)),float(np.percentile(bt[:,i],97.5))] for i,k in enumerate(["gain_steep_pct","early_steep_pct","gain_low_pct","early_low_pct","gain_km2"])}
d_steep=bt[:,0]-bt[:,1]; d_low=bt[:,2]-bt[:,3]
R["boot"]["steep_diff_gain_minus_early"]=[float(np.percentile(d_steep,2.5)),float(np.percentile(d_steep,97.5)),float(d_steep.mean())]
R["boot"]["low_diff_gain_minus_early"]=[float(np.percentile(d_low,2.5)),float(np.percentile(d_low,97.5)),float(d_low.mean())]
# 5 cross-product at 100 m: aggregate DW fraction onto the GHSL grid
tr=gsrc.transform; gH,gW=g.shape[1:]
x0=D; # GHSL grid origin 753200,9637200 ; DW origin 753280..: find offset
dwsrc=rasterio.open(D+"bujumbura_dw_built_and_counts.tif")
def to100(mask):
    out=np.zeros((gH,gW)); 
    ox=round((dwsrc.transform.c-tr.c)/20); oy=round((tr.f-dwsrc.transform.f)/20)
    big=np.zeros((gH*5,gW*5)); h,w=mask.shape
    big[oy:oy+h,ox:ox+w]=mask[:min(h,gH*5-oy),:min(w,gW*5-ox)]
    return big.reshape(gH,5,gW,5).mean((1,3))
fE=to100(((E>=70)&land).astype(float)); fM=to100(((M>=70)&land).astype(float)); fL=to100(((L>=70)&land).astype(float))
fE17=to100(((comp([0,1])>=70)&land).astype(float)); fL21=to100(((comp([4,5])>=70)&land).astype(float))
g=np.nan_to_num(g,nan=0.0); ghsl15=g[1]/1e4; ghsl20=g[2]/1e4; ghsl00=g[0]/1e4
wp15=g[3]; wp20=g[4]
valid=to100(land.astype(float))>0.5
def corr(a,b_): 
    m=valid; return float(np.corrcoef(a[m],b_[m])[0,1])
R["x100"]=dict(corr_dwM_ghsl20=corr(fM,ghsl20),corr_dwE_ghsl15=corr(fE17,ghsl15),
  corr_change_dw_ghsl=corr(fL21-fE17,ghsl20-ghsl15),
  ghsl_gain_1520_km2=float(((ghsl20-ghsl15)[valid].sum())*1e-2),
  dw_gain_1621_km2=float(((fL21-fE17)[valid].sum())*1e-2))
def blk(a,k):
    h,w=(a.shape[0]//k)*k,(a.shape[1]//k)*k; return a[:h,:w].reshape(h//k,k,w//k,k).sum((1,3))
cc={}
for k in (1,5,10):
    a=blk(np.where(valid,fL21-fE17,0),k).ravel(); c=blk(np.where(valid,ghsl20-ghsl15,0),k).ravel(); p=blk(np.where(valid,wp20-wp15,0),k).ravel()
    cc[k]=dict(dw_vs_ghsl=float(np.corrcoef(a,c)[0,1]),dw_vs_pop=float(np.corrcoef(a,p)[0,1]),ghsl_vs_pop=float(np.corrcoef(c,p)[0,1]))
R["x_scale"]=cc
# 6 population
dp=wp20-wp15; 
R["pop"]=dict(total15=float(wp15[valid].sum()),total20=float(wp20[valid].sum()),growth_pct=float((wp20[valid].sum()/wp15[valid].sum()-1)*100),
  built_per_person_15_m2=float(ghsl15[valid].sum()*1e4/wp15[valid].sum()),built_per_person_20_m2=float(ghsl20[valid].sum()*1e4/wp20[valid].sum()),
  corr_dpop_dghsl=corr(dp,ghsl20-ghsl15),corr_dpop_ddw=corr(dp,fL21-fE17))
# share of added population in cells that were <10% built in 2015
low_b=(ghsl15<0.1)&valid
R["pop"]["added_pop_in_lowbuilt_2015_pct"]=float(dp[low_b].sum()/dp[valid].sum()*100)
# GHSL 2000-2020 totals
R["ghsl_tot"]=[float(x[valid].sum()*1e-2) for x in (ghsl00,ghsl15,ghsl20)]
# WorldCover 2020/21
R["wc_built"]=[float((wc20==50).sum()*PX),float((wc21==50).sum()*PX)]
np.savez_compressed("results/masks.npz",gain=gain,loss=loss,Eb=Eb,Lb=Lb,E=E.astype(np.float32),L=L.astype(np.float32),land=land,dist=dist.astype(np.float32),el=el.astype(np.float32),sl=sl.astype(np.float32),
  ghsl15=ghsl15,ghsl20=ghsl20,fE17=fE17,fL21=fL21,dp=dp)
json.dump(R,open("results/results.json","w"),indent=1)
print(json.dumps(R,indent=1))
