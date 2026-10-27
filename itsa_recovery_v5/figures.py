"""Publication plots of observations, fixed PRE target, fitted ITSA and residuals."""
import os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR",str(Path(__file__).resolve().parent/".runtime_matplotlib"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from scipy.stats import t as student_t
from model import design
from recovery import assess,LABELS,COLORS
from raster_outputs import qml

def save(fig,out,name,c):
    for ext in c.FIGURE_FORMATS:
        fig.savefig(out/f"{name}.{ext}",dpi=c.FIGURE_DPI,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def regional_graph(common,per_observation,obs,origin,settings,out,c):
    names=[o.name for o in obs]
    rows=[]
    for sample,o in zip(per_observation,obs):
        row=sample.table([o.name]).to_dict("records")[0] if sample.count else dict(observation=o.name)
        rows.append(row|dict(phase=o.phase,reference_date=o.reference.isoformat()))
    pd.DataFrame(rows).to_csv(out/"available_support_timeseries.csv",index=False)
    if not common.count:
        (out/"NO_COMMON_SUPPORT.md").write_text("No seven/all-date common support. No regional regression plotted.",
                                               encoding="utf-8")
        return
    table=common.table(names)
    table["phase"]=[o.phase for o in obs]
    table["reference_date"]=[o.reference.isoformat() for o in obs]
    table.to_csv(out/"common_support_timeseries.csv",index=False)
    y=table.median_sample.to_numpy()
    result=assess(y[:,None,None],obs,origin,settings)
    x,use,time=design(obs,origin)
    beta=result.model.beta[:,0,0]
    baseline=result.metrics["pre_baseline"].item()
    dates=mdates.date2num([o.reference for o in obs]); event=mdates.date2num(pd.Timestamp(origin))
    with plt.rc_context({"font.family":"serif","font.size":10,"pdf.fonttype":42,"svg.fonttype":"none"}):
        fig,ax=plt.subplots(figsize=(10,5.6),layout="constrained")
        for j,o in enumerate(obs):
            ax.vlines(dates[j],table.q25_sample.iloc[j],table.q75_sample.iloc[j],color=".7",lw=2)
            ax.plot(dates[j],y[j],marker="D" if o.phase=="event" else "o" if o.phase=="pre" else "s",
                color="#245a8d",mfc="white" if o.phase=="event" else "#245a8d",ms=6)
        for phase in ("pre","post"):
            which=np.array([o.phase==phase for o in obs])
            grid=np.linspace(dates[which].min(),event,100) if phase=="pre" else np.linspace(event,dates[which].max(),100)
            t=grid-event
            fitted=beta[0]+beta[1]*t+(beta[2]+beta[3]*t if phase=="post" else 0)
            ax.plot(grid,fitted,color="#245a8d",ls="-" if phase=="pre" else "--",lw=1.7,label=phase.upper()+" fit")
        grid=np.linspace(event,dates.max(),100)
        ax.plot(grid,beta[0]+beta[1]*(grid-event),color=".5",ls=":",label="PRE-trend counterfactual")
        ax.axhline(baseline,color="black",lw=1.1,label="Fixed observed PRE median")
        ax.axhline(baseline-settings.RETURN_TOLERANCE,color=".3",ls="-.",label="Baseline return threshold")
        ax.axvline(event,color="black",ls="--",lw=1,label="Karding / Noru interruption")
        ax.set(xlabel="Pair reference date (UTC; intervals in metadata)",ylabel="Sentinel-1 coherence",
               title="Exploratory ITSA-derived baseline recovery: unfilled primary inputs",ylim=(0,1))
        locator=mdates.AutoDateLocator(minticks=4,maxticks=9)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        ax.grid(axis="y",alpha=.15)
        ax.legend(fontsize=8,ncol=2,loc="upper center")
        save(fig,out,"ITSA_baseline_recovery",c)
        residual=y[use]-x[use]@beta
        fig,ax=plt.subplots(figsize=(9,3.5),layout="constrained")
        ax.axhline(0,color="black",lw=.8)
        ax.plot(dates[use],residual,"o-",color="#245a8d",lw=.8)
        ax.set(xlabel="Pair reference date (UTC)",ylabel="Coherence residual",
               title=f"Exploratory residual diagnostics (n={use.sum()}, df={int(result.model.df.item())})")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b %Y"))
        save(fig,out,"ITSA_residuals",c)
    pd.DataFrame({name:[value.item()] for name,value in result.metrics.items()}).to_csv(
        out/"regional_model_recovery_indicators.csv",index=False)
    xx=x[use]; dof=len(xx)-4
    covariance=np.linalg.pinv(xx.T@xx)*float(residual@residual)/dof
    se=np.sqrt(np.maximum(0,np.diag(covariance)))
    half=student_t.ppf(.975,dof)*se
    statistic=np.divide(beta,se,out=np.full(4,np.nan),where=se>0)
    pd.DataFrame(dict(term=["b0","b1","b2","b3"],estimate=beta,standard_error=se,
        ci95_lower=beta-half,ci95_upper=beta+half,p_value=2*student_t.sf(np.abs(statistic),dof),
        n=int(use.sum()),residual_df=dof)).to_csv(out/"regional_itsa_coefficients.csv",index=False)
    (out/"figure_caption.md").write_text(
        "Actual acquisition-pair data; points are seeded spatial-sample median estimates. "
        "IQR describes spatial variability, not confidence. The event diamond is excluded "
        "from OLS and retained for loss evidence. Recovery target is the observed PRE median, "
        "not the counterfactual. No spatial filling in the primary branch. "
        "Classical OLS uncertainty is conditional on unverified error independence. "
        "The aggregate fit is not a pixel- or zone-level recovery decision.",encoding="utf-8")

def recovery_maps(rasters,out,c):
    import rasterio
    from matplotlib.colors import ListedColormap,BoundaryNorm
    from matplotlib.patches import Patch
    from rasterio.enums import Resampling
    with rasterio.open(rasters/"ITSA_Recovery_Classes.tif") as src:
        h,w=max(1,round(src.height/8)),max(1,round(src.width/8))
        a=src.read(1,out_shape=(h,w),resampling=Resampling.nearest)
        bounds=src.bounds; crs=str(src.crs)
    fig,ax=plt.subplots(figsize=(8,8),layout="constrained")
    cmap=ListedColormap([np.array(COLORS[i])/255 for i in range(6)])
    ax.imshow(a,cmap=cmap,norm=BoundaryNorm(np.arange(-.5,6),6),interpolation="nearest",
              extent=(bounds.left,bounds.right,bounds.bottom,bounds.top))
    ax.set(title="Exploratory ITSA-derived baseline coherence recovery",xlabel=f"Easting ({crs})",ylabel="Northing")
    fig.legend(handles=[Patch(color=np.array(COLORS[i])/255,label=LABELS[i]) for i in range(1,6)],
               loc="outside lower center",ncol=2,fontsize=8,frameon=False)
    save(fig,out,"ITSA_recovery_map",c)
    with rasterio.open(rasters/"baseline_deficit_final.tif") as src:
        deficit=src.read(1,out_shape=(h,w),resampling=Resampling.bilinear,masked=True)
    if deficit.count():
        fig,ax=plt.subplots(figsize=(8,8),layout="constrained")
        im=ax.imshow(deficit,cmap="RdBu_r",vmin=-.2,vmax=.2,interpolation="bilinear",
                     extent=(bounds.left,bounds.right,bounds.bottom,bounds.top))
        ax.set(title="Modeled final deficit relative to fixed PRE baseline",xlabel=f"Easting ({crs})",ylabel="Northing")
        fig.colorbar(im,ax=ax,label="B - modeled POST coherence (display range only)",extend="both")
        save(fig,out,"ITSA_modeled_deficit_map",c)
