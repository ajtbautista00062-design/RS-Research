"""Exact boundary counts and V2 decision audit from existing exported coherence only."""
import sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(1,str(HERE/"preserved_v1"))
import json,argparse,time,warnings
from datetime import datetime,timezone
from types import SimpleNamespace
import numpy as np
import pandas as pd
import rasterio
from sarflow.temporal import Observation,timestamp
from sarflow.spatial import read_zones,mask_for,windows
from sarflow.sampling import Sample
from previous_classifier import assess
from v2_model import design,contrast_se
import config

NAMES=("fitted_level_loss","event_loss","post_slope","final_deficit","pre_noise",
       "pre_spread","v2_loss_threshold","v2_step_guard")

def observations(manifest):
    return [Observation(o["name"],Path(o["path"]),o["phase"],timestamp(o["start"]),timestamp(o["end"]),
            timestamp(o["start"]),timestamp(o["end"]),timestamp(o["reference"])) for o in manifest["observations"]]

def pre_noise(stack,obs):
    indices=[j for j,o in enumerate(obs) if o.phase=="pre"]
    t=np.array([(obs[j].reference-obs[indices[0]].reference).total_seconds()/86400 for j in indices])
    x=np.column_stack((np.ones(len(t)),t))
    data=stack[indices].reshape(len(indices),-1).astype(float)
    out=np.full(data.shape[1],np.nan)
    complete=np.isfinite(data).all(axis=0)
    if complete.any():
        y=data[:,complete]; residual=y-x@np.linalg.pinv(x)@y
        out[complete]=np.sqrt(np.sum(residual**2,axis=0)/(len(t)-2))
    return out.reshape(stack.shape[1:])

def run(source,out):
    manifest=json.loads((source/"run_manifest.json").read_text())
    obs=observations(manifest); c=SimpleNamespace(**manifest["thresholds"])
    raster=source/"Recovery Map Raster File"
    counts=np.zeros(6,dtype="int64"); funnel={}; study=0; mismatch=0; checked=0
    sample=Sample(50000,20220925,len(NAMES))
    controls=Sample(50000,20220926,1)
    started=time.monotonic()
    with rasterio.open(raster/"ITSA_Recovery_Classes.tif") as cls,rasterio.open(raster/"input_coherence.tif") as src:
        geoms=[g for _,g in read_zones(config.Study_Area,cls.crs,config.PROJECT_ROOT)]
        for number,window in enumerate(windows(cls.width,cls.height,256),1):
            inside=mask_for(geoms,(int(window.height),int(window.width)),
                rasterio.windows.transform(window,cls.transform),False)
            if not inside.any():
                continue
            codes=cls.read(1,window=window)
            counts+=np.bincount(codes[inside],minlength=6); study+=int(inside.sum())
            stack=src.read(window=window); stack[stack==src.nodata]=np.nan
            result=assess(stack,obs,manifest["event_origin"],c,inside)
            mismatch+=int((codes[inside]!=result.classes[inside]).sum())
            pre=[j for j,o in enumerate(obs) if o.phase=="pre"]
            post=[j for j,o in enumerate(obs) if o.phase=="post"]
            event=next(j for j,o in enumerate(obs) if o.phase=="event")
            b=result.metrics["pre_baseline"]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore",RuntimeWarning)
                spread=1.4826*np.nanmedian(np.abs(stack[pre]-np.nanmedian(stack[pre],axis=0)),axis=0)
            threshold=np.maximum(c.LOSS_MIN,c.LOSS_SPREAD_MULTIPLIER*np.maximum(spread,c.PRE_SPREAD_FLOOR))
            bg=np.maximum(c.BASELINE_GUARD_FLOOR,c.GUARD_MULTIPLIER*spread/np.sqrt(len(pre)))
            sg=c.GUARD_MULTIPLIER*contrast_se(result.model,np.array([0.,0.,1.,0.]))
            le=result.metrics["event_loss"]; lm=result.metrics["fitted_level_loss"]
            unstable=(result.exclusion==8)
            eligible=inside&result.model.valid&~unstable
            masks={"finite_three_pre":np.isfinite(stack[pre]).all(axis=0),
                   "finite_event":np.isfinite(stack[event]),
                   "finite_three_post":np.isfinite(stack[post]).all(axis=0),
                   "eligible_stable_model":eligible,"model_rejected_unstable":unstable,
                   "observed_loss_raw_above_H":eligible&(le>=threshold),
                   "model_loss_raw_above_H":eligible&(lm>=threshold),
                   "event_guarded_loss_above_H":eligible&(le-bg>=threshold-1e-7),
                   "model_guarded_loss_above_H":eligible&(lm-sg>=threshold-1e-7),
                   "both_guarded_loss":eligible&(le-bg>=threshold-1e-7)&(lm-sg>=threshold-1e-7),
                   "event_upper_below_H":eligible&(le+bg<threshold-1e-7),
                   "model_upper_below_H":eligible&(lm+sg<threshold-1e-7),
                   "both_upper_below_H_class4":eligible&(le+bg<threshold-1e-7)&(lm+sg<threshold-1e-7),
                   "loss_floor_controls_H":inside&(threshold<=c.LOSS_MIN+1e-7),
                   "spread_controls_H":inside&(threshold>c.LOSS_MIN+1e-7),
                   "all_seven_exact_zero":inside&np.all(stack==0,axis=0)}
            for key,mask in masks.items():
                funnel[key]=funnel.get(key,0)+int((inside&mask).sum())
            noise=pre_noise(stack,obs)
            vector=np.stack((lm,le,result.model.beta[1]+result.model.beta[3],
                result.metrics["baseline_deficit_final"],noise,spread,threshold,sg))
            valid=inside&np.isfinite(vector).all(axis=0)
            rr,cc=np.where(valid)
            if len(rr):
                sample.add(vector[:,valid],rr+int(window.row_off),cc+int(window.col_off))
            # PRE-only calibration set: finite PRE, nonzero temporal support; no POST/class selection.
            control=inside&np.isfinite(noise)&np.isfinite(stack[pre]).all(axis=0)&np.any(stack[pre]>0,axis=0)
            rr,cc=np.where(control)
            if len(rr):
                controls.add(noise[None,control],rr+int(window.row_off),cc+int(window.col_off))
            if checked<100:
                x,use,_=design(obs,manifest["event_origin"])
                rr,cc=np.where(result.model.valid&inside)
                for r,col in zip(rr[:min(5,100-checked)],cc[:min(5,100-checked)]):
                    independent=np.linalg.lstsq(x[use],stack[use,r,col].astype(float),rcond=None)[0]
                    if not np.allclose(independent,result.model.beta[:,r,col],atol=1e-7):
                        raise RuntimeError("Independent cached real-pixel check failed.")
                    checked+=1
            if number%100==0:
                print(f"Existing-output diagnostic: tile {number}, study pixels {study}, {time.monotonic()-started:.1f}s",flush=True)
    if mismatch:
        raise RuntimeError(f"V2 reproduction differs at {mismatch} pixels; inspect source configuration.")
    pd.DataFrame([dict(code=i,count=int(n),percent_inside=100*n/study) for i,n in enumerate(counts)]
        ).to_csv(out/"v2_class_counts_independent.csv",index=False)
    pd.DataFrame([dict(condition=k,accepted=v,rejected=study-v,denominator="all boundary pixels")
                  for k,v in funnel.items()]).to_csv(out/"v2_decision_conditions.csv",index=False)
    quantiles=[]
    for j,name in enumerate(NAMES):
        q=np.quantile(sample.values[j],[.05,.25,.5,.75,.95])
        quantiles.append(dict(indicator=name,p05=q[0],p25=q[1],median=q[2],p75=q[3],p95=q[4],
                              sample_pixels=sample.values.shape[1]))
    pd.DataFrame(quantiles).to_csv(out/"indicator_distributions.csv",index=False)
    calibration=dict(method="median pooled PRE-only detrended residual noise proxy",
        calibrated_floor=float(np.median(controls.values)),sample_pixels=controls.values.shape[1],
        population_pixels=controls.count,validated=False,independent_ground_truth_available=False,
        limitation="Three PRE observations give one residual degree of freedom; this is not validated false-alarm calibration")
    (out/"calibration.json").write_text(json.dumps(calibration,indent=2),encoding="utf-8")
    info=dict(source=str(source.resolve()),study_pixels=study,class_counts=counts.tolist(),conditions=funnel,
              independently_checked_real_pixels=checked,v2_reproduction_mismatches=mismatch,calibration=calibration)
    (out/"audit.json").write_text(json.dumps(info,indent=2),encoding="utf-8")
    np.savez_compressed(out/"sampled_distributions.npz",values=sample.values,names=np.array(NAMES))
    import os
    os.environ.setdefault("MPLCONFIGDIR",str(HERE/".runtime_matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout="constrained")
    for ax,j in zip(axes.ravel(),range(4)):
        ax.hist(sample.values[j],bins=70,color="#245a8d")
        ax.axvline(0,color="black",lw=.7)
        if j in (0,1):
            ax.axvline(.10,color=".3",ls="--",label="Previous absolute floor")
            ax.legend(fontsize=8)
        ax.set(xlabel=NAMES[j],ylabel="Spatial sample count")
    fig.suptitle("Existing V2 output diagnostics (sampled distributions; exact class counts separate)")
    for ext in ("png","pdf","svg"):
        fig.savefig(out/f"v2_indicator_distributions.{ext}",dpi=300)
    plt.close(fig)
    return info

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,required=True); parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=False)
    result=run(args.source,args.output)
    print(json.dumps(result,indent=2))
