"""Detection conditions, exact sensitivity counts and bounded-sample distributions."""
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from recovery import COLORS,LABELS

def threshold_figure(table,out,formats=("png","pdf","svg"),dpi=600):
    """Show whole-boundary shares and rare recovery classes without hiding either."""
    import matplotlib.pyplot as plt
    names=list(table.scenario)
    x=np.arange(len(names))
    labels=[f"{int(name.split('_')[-1])/100:.2f}" for name in names]
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),layout="constrained")
    base=np.zeros(len(names))
    short={1:"Recovered",2:"Partial",3:"Not recovered / degraded",4:"No detectable loss",5:"Indeterminate"}
    for code in range(1,6):
        values=table[f"percent_{code}"].to_numpy()
        axes[0].bar(x,values,bottom=base,color=np.array(COLORS[code])/255,label=short[code])
        base+=values
    for code,marker in ((1,"o"),(2,"s"),(3,"^")):
        axes[1].plot(x,table[f"class_{code}"],marker=marker,color=np.array(COLORS[code])/255,label=short[code])
    axes[0].set(ylabel="Boundary pixels (%)",ylim=(0,100),title="All classification decisions")
    axes[1].set(ylabel="Pixel count",title="Recovery classes shown separately",ylim=(0,None))
    for ax in axes:
        ax.set_xticks(x,labels)
        ax.set_xlabel("Candidate minimum modeled loss (coherence)")
        ax.spines[["top","right"]].set_visible(False)
    axes[0].legend(fontsize=8,loc="upper center",bbox_to_anchor=(.5,-.23),ncol=2,frameon=False)
    axes[1].legend(fontsize=8,frameon=False)
    fig.suptitle("Threshold sensitivity: candidate thresholds require validation",fontsize=12)
    for ext in formats:
        fig.savefig(out/f"threshold_sensitivity.{ext}",dpi=dpi,bbox_inches="tight")
    plt.close(fig)

def condition_masks(stack,obs,result,c):
    inside=result.exclusion!=1
    pre=[i for i,o in enumerate(obs) if o.phase=="pre"]
    post=[i for i,o in enumerate(obs) if o.phase=="post"]
    event=next(i for i,o in enumerate(obs) if o.phase=="event")
    loss=result.metrics["fitted_level_loss"]
    se=result.metrics["level_change_se"]
    critical=student_t.ppf(1-c.DETECTION_ALPHA,np.maximum(result.model.df,1))
    masks={"sufficient_PRE":np.isfinite(stack[pre]).sum(axis=0)>=c.MIN_PRE,
        "finite_event":np.isfinite(stack[event]),
        "sufficient_POST":np.isfinite(stack[post]).sum(axis=0)>=c.MIN_POST,
        "model_available":result.model.valid,
        "unstable_model":(result.quality&512)>0,
        "model_lower_loss_above_threshold":loss>=result.metrics["loss_detection_boundary"]-1e-7,
        "model_upper_loss_below_threshold":loss+critical*se<c.LOSS_MIN-1e-7,
        "event_conflict":(result.quality&16)>0,
        "event_neutral_support":(result.quality&2048)>0,
        "ambiguous":(result.quality&8)>0,
        "all_zero_quarantined":(result.quality&8192)>0,
        "return_unobserved":result.exclusion==11,
        "operational_confident":result.metrics['evidence_status']==1,
        "operational_provisional":result.metrics['evidence_status']==2,
        "modeled_point_loss_above_effect_threshold":loss>=result.metrics['loss_effect_threshold']+getattr(c,'EFFECT_MARGIN',0)}
    from recovery import EXCLUSION
    masks.update({f'rejection_{code}_{label}':result.exclusion==code for code,label in EXCLUSION.items() if code>1})
    return masks

def write_plots(common,obs,origin,settings,summary,out,cfg):
    import matplotlib.pyplot as plt
    from recovery import assess
    if common.count:
        result=assess(common.values[:,None,:],obs,origin,settings)
        metrics=("fitted_level_loss","observed_event_change","post_slope","baseline_deficit_final",
                 "modeled_recovery_progress","level_change_se")
        fig,axes=plt.subplots(3,2,figsize=(10,10),layout="constrained")
        rows=[]
        for name,ax in zip(metrics,axes.ravel()):
            values=result.metrics[name].ravel()
            values=values[np.isfinite(values)]
            if not len(values):
                continue
            q=np.quantile(values,[.05,.25,.5,.75,.95])
            rows.append(dict(indicator=name,sample_count=len(values),p05=q[0],p25=q[1],
                             median=q[2],p75=q[3],p95=q[4]))
            ax.hist(values,bins=70,color="#245a8d")
            ax.axvline(0,color="black",lw=.8)
            ax.set(xlabel=name,ylabel="Spatial sample count")
        fig.suptitle("Model diagnostics on common-support spatial sample")
        for ext in cfg.FIGURE_FORMATS:
            fig.savefig(out/f"model_distributions.{ext}",dpi=cfg.FIGURE_DPI)
        plt.close(fig)
        pd.DataFrame(rows).to_csv(out/"model_distributions.csv",index=False)
    names=[name for name in summary.counts if name.startswith("loss_")]
    if names:
        table=pd.DataFrame([dict(scenario=name,
            **{f"class_{i}":int(v) for i,v in enumerate(summary.counts[name])},
            **{f"percent_{i}":100*v/summary.inside for i,v in enumerate(summary.counts[name])}) for name in names])
        threshold_figure(table,out,cfg.FIGURE_FORMATS,cfg.FIGURE_DPI)
