import rasterio, numpy as np, warnings
warnings.filterwarnings("ignore")
D="data/"
dw=rasterio.open(D+"bujumbura_dw_built_and_counts.tif").read().astype(float)
b=dw[:10]; n=dw[10:]
wc20,wc21,el,sl=rasterio.open(D+"bujumbura_landcover_and_terrain.tif").read().astype(float)
land=(wc21!=80)&(wc20!=80)
yrs=list(range(2016,2026))
# early/late composites (mean of yearly built prob, years with >=3 obs)
def comp(idx):
    w=(n[idx]>=3).astype(float); return (b[idx]*w).sum(0)/np.maximum(w.sum(0),1)
early=comp([0,1,2]); mid=comp([3,4,5]); late=comp([7,8,9])
print("early/late mean on land",early[land].mean(),late[land].mean())
# agreement of DW composite 2019-2021 vs WC2020 built
c2020=comp([3,4,5])
for t in (30,40,50,60,70):
    p=(c2020>=t)&land; r=(wc20==50)&land
    tp=(p&r).sum(); print("tau",t,"DW km2",p.sum()*4e-4,"prec",round(tp/p.sum(),3),"rec",round(tp/r.sum(),3),"F1",round(2*tp/(p.sum()+r.sum()),3))
# per-pixel noise: sd across years among pixels with early mean>70
m=(early>70)&land
print("stable built core px",m.sum(),"sd across years",np.nanmean(b[:, m].std(0)))
for t in (30,50,70):
    E=(early>=t)&land; L=(late>=t)&land
    print("tau",t,"early km2",E.sum()*4e-4,"late",L.sum()*4e-4,"gain",(L&~E).sum()*4e-4,"loss",(E&~L).sum()*4e-4)
g=rasterio.open(D+"bujumbura_ghsl_worldpop_100m.tif").read().astype(float)
print("ghsl band stats",[(float(np.nanmin(x)),float(np.nanmax(x))) for x in g])
