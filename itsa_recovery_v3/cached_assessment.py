"""Requested V3 threshold diagnostic using saved, audited V2 unfilled inputs only."""
import sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(1,str(HERE/"preserved_v1"))
import json,copy,argparse,time,logging
from contextlib import ExitStack
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from scipy.stats import t as student_t
from sarflow.spatial import read_zones,mask_for,windows
from sarflow.temporal import file_digest
from recovery import assess,METRICS,QUALITY,EXCLUSION
from model import contrast_se
from raster_outputs import Outputs
from validation import independent_check
from investigate_v2 import observations
import config

def evidence(result,c,inside):
    boundary=result.metrics["loss_detection_boundary"]
    loss=result.metrics["fitted_level_loss"]
    critical=student_t.ppf(1-c.DETECTION_ALPHA,np.maximum(result.model.df,1))
    se=result.metrics["level_change_se"]
    status=np.zeros(inside.shape,dtype="uint8")
    status[inside]=3
    status[inside&(loss>=boundary-1e-7)]=1
    status[inside&(loss+critical*se<c.LOSS_MIN-1e-7)]=2
    bad=inside&np.isin(result.exclusion,[2,3,4,5,7,8,9,10])
    status[bad]=4
    return status

def run(source,out,calibration):
    manifest=json.loads((source/"run_manifest.json").read_text())
    if file_digest(Path(config.Study_Area))!=manifest["boundary_sha256"]:
        raise ValueError("Boundary changed since V2; cached comparison would not be like-for-like.")
    obs=observations(manifest)
    settings=config.analysis_settings()
    settings.LOSS_MIN=float(calibration["calibrated_floor"])
    settings.EVENT_EVIDENCE_MIN=settings.LOSS_MIN
    settings.ZERO_UNRESOLVED=manifest["zero_semantics"]=="unresolved"
    thresholds=list(config.DETECTION_SWEEP)
    cands={f"loss_{value:.2f}".replace(".",""):copy.copy(settings) for value in thresholds}
    for case,value in zip(cands.values(),thresholds):
        case.LOSS_MIN=value
    folder=source/"Recovery Map Raster File"
    count={name:np.zeros(6,dtype="int64") for name in ["v2_previous","v3_provisional",*cands,"zero_excluded_unfilled","zero_excluded_filled"]}
    matrices={name:np.zeros((6,6),dtype="int64") for name in ["v3_provisional",*cands,"zero_excluded_unfilled","zero_excluded_filled","fill_only"]}
    conditions={}; total=0; changed=0; checks=0; started=time.monotonic()
    with ExitStack() as resources:
        src=resources.enter_context(rasterio.open(folder/"input_coherence.tif","r"))
        prior=resources.enter_context(rasterio.open(folder/"ITSA_Recovery_Classes.tif","r"))
        flags=resources.enter_context(rasterio.open(folder/"input_provenance.tif","r"))
        geoms=[g for _,g in read_zones(config.Study_Area,src.crs,config.PROJECT_ROOT)]
        tags=dict(mode="cached-input threshold diagnostic",source_run=str(source.resolve()),
            analysis_status="EXPLORATORY; calibration and acquisition comparability unvalidated",
            loss_calibration=json.dumps(calibration),event_origin=manifest["event_origin"])
        primary=Outputs(resources,out/"Recovery Map Raster File",src,obs,config,tags,full=True)
        outputs={name:Outputs(resources,out/"Threshold Sweep"/name,src,obs,config,
                     tags|dict(threshold=str(case.LOSS_MIN)),full=False) for name,case in cands.items()}
        for name in ("zero_excluded_unfilled","zero_excluded_filled"):
            outputs[name]=Outputs(resources,out/"Spatial Sensitivity"/name,src,obs,config,tags|dict(scenario=name),full=False)
        excluded_folder=source/"Sensitivity"/"zero_excluded_filled"
        corrected=resources.enter_context(rasterio.open(excluded_folder/"input_coherence.tif","r"))
        corrected_flags=resources.enter_context(rasterio.open(excluded_folder/"input_provenance.tif","r"))
        profile=prior.profile|dict(nodata=0)
        evidence_writers={name:resources.enter_context(rasterio.open(
            out/"Threshold Sweep"/name/"ITSA_Disturbance_Evidence.partial.tif","w",**profile)) for name in cands}
        for dst in evidence_writers.values():
            dst.update_tags(code_1="Modeled loss detected",code_2="Modeled upper loss below minimum change",
                            code_3="Ambiguous",code_4="Conflicting/rejected",code_0="Outside")
            dst.write_colormap(1,{0:(0,0,0,0),1:(213,94,0,255),2:(0,114,178,255),3:(150,150,150,255),4:(150,100,160,255)})
        sensitive=resources.enter_context(rasterio.open(out/"Threshold_Sensitivity.partial.tif","w",**profile))
        for block,window in enumerate(windows(src.width,src.height,256),1):
            shape=(int(window.height),int(window.width))
            inside=mask_for(geoms,shape,rasterio.windows.transform(window,src.transform),False)
            if not inside.any():
                primary.empty(shape,window)
                for output in outputs.values():
                    output.empty(shape,window)
                for dst in [sensitive,*evidence_writers.values()]:
                    dst.write(np.zeros(shape,dtype="uint8"),1,window=window)
                continue
            stack=src.read(window=window); stack[stack==src.nodata]=np.nan
            provenance=flags.read(window=window)
            if np.any(provenance==2):
                raise ValueError("Primary cached stack unexpectedly contains filled pixels.")
            old=prior.read(1,window=window)
            result=assess(stack,obs,manifest["event_origin"],settings,inside)
            if checks<100:
                checks+=independent_check(stack,obs,manifest["event_origin"],settings,result,min(5,100-checks))
            primary.write(result,stack,provenance,window)
            total+=int(inside.sum())
            count["v2_previous"]+=np.bincount(old[inside],minlength=6)
            alternatives={"v3_provisional":result.classes}
            for name,case in cands.items():
                alternative=assess(stack,obs,manifest["event_origin"],case,inside)
                outputs[name].write(alternative,stack,provenance,window)
                alternatives[name]=alternative.classes
                evidence_writers[name].write(evidence(alternative,case,inside),1,window=window)
            vals=corrected.read(window=window); vals[vals==corrected.nodata]=np.nan
            pf=corrected_flags.read(window=window)
            unfilled=vals.copy(); unfilled[pf==2]=np.nan
            unfilled_flags=pf.copy(); unfilled_flags[pf==2]=0
            spatial_results={}
            for name,arr,prov in (("zero_excluded_unfilled",unfilled,unfilled_flags),("zero_excluded_filled",vals,pf)):
                alternative=assess(arr,obs,manifest["event_origin"],settings,inside,(prov==2).any(axis=0))
                outputs[name].write(alternative,arr,prov,window)
                alternatives[name]=alternative.classes
                spatial_results[name]=alternative.classes
            for name,codes in alternatives.items():
                count[name]+=np.bincount(codes[inside],minlength=6)
                matrices[name]+=np.bincount(old[inside]*6+codes[inside],minlength=36).reshape(6,6)
            matrices["fill_only"]+=np.bincount(
                spatial_results["zero_excluded_unfilled"][inside]*6+spatial_results["zero_excluded_filled"][inside],
                minlength=36).reshape(6,6)
            variations=np.stack([alternatives[name] for name in cands])
            varies=inside&np.any(variations!=variations[0],axis=0)
            changed+=int(varies.sum())
            sensitive.write(np.where(inside,np.where(varies,2,1),0).astype("uint8"),1,window=window)
            pre=[i for i,o in enumerate(obs) if o.phase=="pre"]
            post=[i for i,o in enumerate(obs) if o.phase=="post"]
            event=next(i for i,o in enumerate(obs) if o.phase=="event")
            eligible=inside&result.model.valid&~np.isin(result.exclusion,[2,3,4,5,8,9,10])
            status=evidence(result,settings,inside)
            masks={"finite_PRE":np.isfinite(stack[pre]).all(axis=0),
                "finite_EVENT":np.isfinite(stack[event]),"finite_POST":np.isfinite(stack[post]).all(axis=0),
                "stable_supported_model":eligible,"credible_ITSA_loss":status==1,
                "bounded_model_no_loss":status==2,
                "event_conflict":(result.quality&16)>0,"neutral_event_support":(result.quality&2048)>0,
                "uncertainty":(result.quality&8)>0,"all_zero_quarantined":(result.quality&8192)>0,
                "unstable":(result.quality&512)>0,"outside_return_interval":result.exclusion==11}
            for key,value in masks.items():
                conditions[key]=conditions.get(key,0)+int((inside&value).sum())
            if block%100==0:
                print(f"Cached threshold diagnostic tile {block}, {total} boundary pixels, {time.monotonic()-started:.1f}s",flush=True)
        primary.close()
        for output in outputs.values():
            output.close()
        for dst in [sensitive,*evidence_writers.values()]:
            dst.close()
        for output in [primary,*outputs.values()]:
            output.verify(); output.finalize()
    for path in out.rglob("*.partial.tif"):
        path.rename(path.with_name(path.name.replace(".partial.tif",".tif")))
    # Independent frequency check from written classes, without requesting GDAL statistics.
    observed=np.zeros(6,dtype="int64")
    with rasterio.open(out/"Recovery Map Raster File"/"ITSA_Recovery_Classes.tif") as cls,rasterio.open(out/"Recovery Map Raster File"/"exclusion_reason.tif") as exclusion:
        for _,window in cls.block_windows(1):
            a=cls.read(1,window=window); inside=exclusion.read(1,window=window)!=1
            observed+=np.bincount(a[inside],minlength=6)
    if not np.array_equal(observed,count["v3_provisional"]):
        raise RuntimeError("Independent V3 class count verification failed.")
    rows=[]
    for name,codes in count.items():
        baseline=matrices.get(name)
        changes=int(baseline.sum()-np.trace(baseline)) if baseline is not None else 0
        rows.append(dict(scenario=name,study_pixels=total,**{f"class_{i}":int(n) for i,n in enumerate(codes)},
            **{f"percent_{i}":100*n/total for i,n in enumerate(codes)},changed_from_v2=changes,
            independent_evidence_agreement="not available",calibration_validated=False))
    pd.DataFrame(rows).to_csv(out/"threshold_class_counts.csv",index=False)
    pd.DataFrame([dict(comparison=name,old_class=i,new_class=j,pixels=int(n))
        for name,matrix in matrices.items() for i,row in enumerate(matrix) for j,n in enumerate(row)]
        ).to_csv(out/"transition_counts.csv",index=False)
    pd.DataFrame([dict(condition=name,accepted=n,rejected=total-n) for name,n in conditions.items()]
        ).to_csv(out/"v3_decision_conditions.csv",index=False)
    report=dict(source=str(source.resolve()),study_pixels=total,calibration=calibration,
        threshold_sensitive_pixels=changed,threshold_sensitive_percent=100*changed/total,
        independent_real_pixel_checks=checks,independent_class_count_check="passed",
        full_source_workflow_rerun=False,authoritative=False,conditions=conditions)
    (out/"assessment_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    from figures import recovery_maps
    graphs=out/"Diagnostic Plots"; graphs.mkdir()
    recovery_maps(out/"Recovery Map Raster File",graphs,config)
    from diagnostics import threshold_figure
    table=pd.DataFrame(rows)
    threshold_figure(table[table.scenario.isin(cands)],graphs,config.FIGURE_FORMATS,config.FIGURE_DPI)
    return report

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--calibration",type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    report=run(args.source,args.output,json.loads(args.calibration.read_text()))
    print(json.dumps(report,indent=2))
