import numpy as np, json, rasterio, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":8,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":300})
R=json.load(open("results/results.json")); M=np.load("results/masks.npz")
D="data/"
FD="paper/figures/"
C1,C2,C3="#1f4e79","#c0504d","#6a8f3c"
# fig 1 series
fig,ax=plt.subplots(figsize=(3.4,2.4)); yrs=np.arange(2016,2026)
ax.plot(yrs,R["yearly"]["50"],"o-",c=C1,ms=3,label=r"DW built, $\tau$=50")
ax.plot(yrs,R["yearly"]["70"],"s-",c=C2,ms=3,label=r"DW built, $\tau$=70")
ax.plot([2020,2021],R["wc_built"],"^",c=C3,ms=5,label="WorldCover (v100, v200)")
ax.plot([2000,2015,2020],R["ghsl_tot"],"D",c="k",ms=3,label="GHSL built")
ax.set_xlim(1999,2026); ax.set_xlabel("Year"); ax.set_ylabel(r"Built area (km$^2$)"); ax.legend(frameon=False,fontsize=6.5,loc="upper left")
fig.tight_layout(); fig.savefig(FD+"fig_series.pdf"); plt.close()
# fig 2 maps
s2=rasterio.open(D+"bujumbura_s2_truecolour_2017_2025.tif").read().astype(float)
def rgb(a): 
    x=np.moveaxis(a,0,-1).astype(float); out=np.zeros_like(x)
    for k in range(3):
        lo,hi=np.percentile(x[...,k],(2,98)); out[...,k]=np.clip((x[...,k]-lo)/(hi-lo),0,1)
    return out**0.9
ext=[0,24.56,0,26.64]  # km (1228*20 m, 1332*20 m)
fig,axs=plt.subplots(1,3,figsize=(7.2,2.9))
for ax,a,t in zip(axs[:2],(s2[:3],s2[3:]),("Sentinel-2, 2017","Sentinel-2, 2025")):
    ax.imshow(rgb(a),extent=ext); ax.set_title(t,fontsize=8)
base=np.where(M["land"],0.88,0.7); ax=axs[2]; ax.imshow(base,cmap="gray",vmin=0,vmax=1,extent=ext)
ov=np.zeros(M["gain"].shape+(4,)); ov[M["Eb"]]=[.35,.35,.35,1]; ov[M["gain"]]=[.85,.1,.1,1]; ov[M["loss"]]=[.1,.3,.85,1]
ax.imshow(ov,extent=ext); ax.set_title("Early built (grey), gain (red), loss (blue)",fontsize=7.5)
for ax in axs: ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig(FD+"fig_maps.pdf"); plt.close()
# fig 3 terrain
fig,axs=plt.subplots(1,3,figsize=(7.2,2.1))
for ax,k,t in zip(axs,("slope","elev","dist"),("Slope (deg)","Elevation (m)","Distance to 2016-18 built-up")):
    rows=R[k]; ax.bar(range(len(rows)),[r["rate"] for r in rows],color=C1,width=.7)
    ax.set_xticks(range(len(rows))); ax.set_xticklabels([r["bin"] for r in rows],fontsize=6.5); ax.set_xlabel(t); 
    for i,r in enumerate(rows): ax.text(i,r["rate"],f'{r["rate"]:.1f}',ha="center",va="bottom",fontsize=6)
axs[0].set_ylabel("Share of unbuilt land\nthat became built (%)")
fig.tight_layout(); fig.savefig(FD+"fig_terrain.pdf"); plt.close()
# fig 4 agreement
v=M["land"].shape
fig,axs=plt.subplots(1,3,figsize=(7.2,2.3))
fE,fL,g15,g20=M["fE17"],M["fL21"],M["ghsl15"],M["ghsl20"]
m=(g20>0)|(fL>0)
axs[0].hexbin(g20[m],(fL)[m],gridsize=35,bins="log",cmap="Blues",mincnt=1); axs[0].set_xlabel("GHSL 2020 built fraction"); axs[0].set_ylabel("DW 2020-21 built fraction"); axs[0].set_title("Levels (100 m)",fontsize=8)
a=(fL-fE)[m]; c=(g20-g15)[m]
axs[1].hexbin(c,a,gridsize=35,bins="log",cmap="Blues",mincnt=1); axs[1].set_xlabel("GHSL change 2015-2020"); axs[1].set_ylabel("DW change 2016-17 to 2020-21"); axs[1].set_title("Change (100 m)",fontsize=8)
ks=["1","5","10"]; xs=np.arange(3); w=.27
for j,(key,lab,col) in enumerate((("dw_vs_ghsl","DW vs GHSL",C1),("dw_vs_pop","DW vs WorldPop",C2),("ghsl_vs_pop","GHSL vs WorldPop",C3))):
    axs[2].bar(xs+(j-1)*w,[R["x_scale"][k][key] for k in ks],w,color=col,label=lab)
axs[2].set_xticks(xs); axs[2].set_xticklabels(["100 m","500 m","1 km"]); axs[2].set_ylabel("Pearson r of change"); axs[2].axhline(0,c="k",lw=.5); axs[2].set_ylim(-0.2,0.6); axs[2].legend(frameon=False,fontsize=6,loc="upper right"); axs[2].set_title("Change agreement vs scale",fontsize=8)
fig.tight_layout(); fig.savefig(FD+"fig_agree.pdf"); plt.close()
print("ok")
