"""Windowed pilot/full execution. Called only by the explicit run command."""
from contextlib import ExitStack
from pathlib import Path
import json,logging,time,copy
import numpy as np
import pandas as pd
import rasterio
from sarflow.spatial import inventory,read_zones,mask_for,windows
from sarflow.sampling import Sample
from sarflow.recovery import baseline_recovery as legacy_classifier
from sarflow.zonal import ZoneStats
from preparation import open_readers,read_window
from recovery import assess,METRICS,QUALITY,metric_unit
from raster_outputs import Outputs
from reporting import Summary
from validation import variants,independent_check
from previous_classifier import assess as assess_v2
from types import SimpleNamespace
from diagnostics import condition_masks,write_plots

def execute(cfg,obs,root,metadata,legacy_settings):
    reports=root/"Reports"; reports.mkdir()
    graphs=root/"Time Series Graphs"; graphs.mkdir()
    spatial_mode=cfg.CLASSIFICATION_MODE!="pixel_itsa"
    rasters=root/("Pixel ITSA Reference" if cfg.EXPORT_PIXEL_LEVEL_REFERENCE else "Internal Pixel Evidence") if spatial_mode else root/"Recovery Map Raster File"
    settings=cfg.analysis_settings()
    summary=Summary(); checks=[]; valid_models=0
    with ExitStack() as resources:
        srcs=[resources.enter_context(rasterio.open(o.path,"r")) for o in obs]
        ref=next(src for src,o in zip(srcs,obs) if o.phase=="event")
        if cfg.REFERENCE_RASTER:
            ref=resources.enter_context(rasterio.open(cfg.REFERENCE_RASTER,"r"))
        zones=read_zones(cfg.Study_Area,ref.crs,cfg.PROJECT_ROOT)
        geoms=[g for _,g in zones]
        primary_reader=cfg.reader_settings(exclude_zero=cfg.ZERO_SEMANTICS=="background")
        pd.DataFrame(inventory(srcs,obs,primary_reader,ref)).to_csv(reports/"source_inventory.csv",index=False)
        primary_vrts=open_readers(resources,srcs,ref,primary_reader)
        excluded_reader=cfg.reader_settings(exclude_zero=True)
        excluded_vrts=open_readers(resources,srcs,ref,excluded_reader) if cfg.RUN_SENSITIVITY else None
        filled_reader=cfg.reader_settings(exclude_zero=True,filled=True)
        primary=Outputs(resources,rasters,ref,obs,cfg,metadata,full=True)
        scenarios=variants(settings,obs) if cfg.RUN_SENSITIVITY else []
        extra_names=[v[0] for v in scenarios]+(["zero_excluded_unfilled","zero_excluded_filled"] if cfg.RUN_SENSITIVITY else [])
        scenario_settings={name:vars(variant) for name,variant,_ in scenarios}
        sensitivity={}
        for name in extra_names:
            tags=metadata|dict(scenario=name)
            tags["thresholds"]=json.dumps(scenario_settings.get(name,vars(settings)))
            tags["pair_reference"]="midpoint" if name=="midpoint_indexing" else cfg.REFERENCE_DATE
            tags["zero_policy"]="excluded as sensitivity assumption" if name.startswith("zero_excluded") else cfg.ZERO_SEMANTICS
            tags["primary_correction"]="nearest-small-hole sensitivity" if name=="zero_excluded_filled" else "none"
            sensitivity[name]=Outputs(resources,root/"Sensitivity"/name,ref,obs,cfg,tags,
                                      full=False,provenance=name=="zero_excluded_filled")
        previous_settings=SimpleNamespace(**(vars(settings)|dict(LOSS_MIN=.10)))
        previous_v2=find_previous_v2(cfg,ref,resources)
        if previous_v2:
            previous_raster,previous_manifest=previous_v2
            previous_settings=SimpleNamespace(**previous_manifest["thresholds"])
            metadata["previous_v2_reference"]=str(previous_raster.name)
        existing=find_existing(cfg,ref,resources)
        if existing:
            metadata["existing_legacy_reference"]=str(existing.name)
        summary.grid=ref.width*ref.height
        grid=dict(crs=str(ref.crs),transform=tuple(ref.transform),shape=ref.shape)
        common=Sample(cfg.SAMPLE_PIXELS,cfg.RANDOM_SEED,len(obs))
        samples=[Sample(cfg.SAMPLE_PIXELS,cfg.RANDOM_SEED+j,1) for j in range(len(obs))]
        numeric_metrics=tuple(name for name in METRICS if metric_unit(name)!='code')
        zonal={name:ZoneStats(name,numeric_metrics) for name,_ in zones}
        zeros=np.zeros(len(obs),dtype="int64")
        started=last=time.monotonic()
        for block,window in enumerate(windows(ref.width,ref.height,cfg.WINDOW_SIZE),1):
            shape=(int(window.height),int(window.width))
            inside=mask_for(geoms,shape,rasterio.windows.transform(window,ref.transform),cfg.ALL_TOUCHED)
            if not inside.any():
                primary.empty(shape,window)
                for output in sensitivity.values():
                    output.empty(shape,window)
                continue
            stack,flags=read_window(primary_vrts,srcs,ref,window,primary_reader,geoms)
            result=assess(stack,obs,cfg.EVENT_ORIGIN,settings,inside)
            summary.inside+=int(inside.sum())
            valid_models+=int((inside&result.model.valid).sum())
            zeros+=((stack==0)&inside).sum(axis=(1,2))
            summary.add("primary_unfilled",result.classes,inside)
            summary.add("guarded_rule_same_fit",np.nan_to_num(result.metrics['guarded_reference_class'],nan=0).astype('uint8'),inside,result.classes)
            summary.reasons+=np.bincount(result.exclusion[inside],minlength=max(summary.reasons.size,12))
            for bit in QUALITY:
                summary.quality[str(bit)]+=int(((result.quality[inside]&bit)>0).sum())
            remaining=cfg.REAL_PIXEL_CHECK_LIMIT-summary.real_checks
            if remaining>0:
                checked=independent_check(stack,obs,cfg.EVENT_ORIGIN,settings,result,min(5,remaining))
                rr,cc=np.where(result.model.valid)
                for r,col in zip(rr[:checked],cc[:checked]):
                    checks.append(dict(row=int(window.row_off)+int(r),column=int(window.col_off)+int(col),
                        coherence_values=json.dumps(stack[:,r,col].tolist()),
                        coefficients=json.dumps(result.model.beta[:,r,col].tolist())))
                summary.real_checks+=checked
            for key,value in condition_masks(stack,obs,result,settings).items():
                summary.conditions[key]+=int((inside&value).sum())
            primary.write(result,stack,flags,window)
            old_v2=assess_v2(stack,obs,cfg.EVENT_ORIGIN,previous_settings,inside)
            summary.add("previous_v2_rule_same_inputs",old_v2.classes,inside,result.classes)
            if previous_v2:
                summary.add("previous_v2_raster",previous_raster.read(1,window=window),inside,result.classes)
            old,_=legacy_classifier(stack,obs,legacy_settings)
            summary.add("legacy_rule_same_primary_inputs",old,inside,result.classes)
            if existing:
                previous=existing.read(1,window=window)
                if not np.isin(previous,range(6)).all():
                    raise RuntimeError("Existing reference contains unexpected class codes.")
                summary.add("existing_legacy_output_reference",previous,inside,result.classes)
            sweep_classes=[]
            for name,variant,vobs in scenarios:
                alternative=assess(stack,vobs,cfg.EVENT_ORIGIN,variant,inside)
                summary.add(name,alternative.classes,inside,result.classes)
                sensitivity[name].write(alternative,stack,flags,window)
                if name.startswith("loss_"):
                    sweep_classes.append(alternative.classes)
            if sweep_classes:
                summary.threshold_sensitive+=int((inside&np.any(
                    np.stack(sweep_classes)!=sweep_classes[0],axis=0)).sum())
            if cfg.RUN_SENSITIVITY:
                for name,reader in (("zero_excluded_unfilled",excluded_reader),
                                    ("zero_excluded_filled",filled_reader)):
                    a,provenance=read_window(excluded_vrts,srcs,ref,window,reader,geoms)
                    alternative=assess(a,obs,cfg.EVENT_ORIGIN,settings,inside,(provenance==2).any(axis=0))
                    alternative.metrics['estimated_observation_count']=np.where(inside,(provenance==2).sum(axis=0),np.nan).astype('float32')
                    summary.add(name,alternative.classes,inside,result.classes)
                    sensitivity[name].write(alternative,a,provenance,window)
                    if reader.FILL_SMALL_HOLES:
                        fills=(provenance==2).sum(axis=(1,2))
                        summary.filled=(np.zeros(len(obs),dtype="int64") if summary.filled is None else summary.filled)+fills
                        legacy,_=legacy_classifier(a,obs,legacy_settings)
                        summary.add("legacy_zero_excluded_filled",legacy,inside,result.classes)
            complete=np.isfinite(stack).all(axis=0)
            rr,cc=np.where(complete)
            if len(rr):
                common.add(stack[:,complete],rr+int(window.row_off),cc+int(window.col_off))
            for j,sample in enumerate(samples):
                rr,cc=np.where(np.isfinite(stack[j]))
                if len(rr):
                    sample.add(stack[j:j+1,np.isfinite(stack[j])],rr+int(window.row_off),cc+int(window.col_off))
            metric_stack=np.stack([result.metrics[name] for name in numeric_metrics])
            for name,g in zones:
                zonal[name].add(metric_stack,mask_for([g],shape,rasterio.windows.transform(window,ref.transform),
                                                    cfg.ALL_TOUCHED),result.classes)
            if time.monotonic()-last>10:
                logging.info("Window %d; inside pixels %d; elapsed %.1f min",block,summary.inside,(time.monotonic()-started)/60)
                last=time.monotonic()
        if not summary.inside:
            raise ValueError("Study boundary does not overlap reference grid.")
        for output in [primary,*sensitivity.values()]:
            output.close()
        logging.info("Rereading outputs for metadata, numerical encoding and class-code checks.")
        for output in [primary,*sensitivity.values()]:
            output.verify()
        verify_decisions(primary,ref,obs,cfg,settings,summary)
        for output in [primary,*sensitivity.values()]:
            output.finalize()
    metadata["zero_cells_inside_primary_by_observation"]=zeros.tolist()
    if summary.filled is not None:
        summary.filled=summary.filled.tolist()
    pd.DataFrame(checks).to_csv(reports/"independent_real_pixel_checks.csv",index=False)
    pd.DataFrame([r for z in zonal.values() for r in z.rows()]).to_csv(reports/"zonal_indicator_statistics.csv",index=False)
    pd.DataFrame([z.classes_row() for z in zonal.values()]).to_csv(reports/"zonal_class_counts.csv",index=False)
    summary.export(reports,metadata)
    from figures import regional_graph,recovery_maps
    regional_graph(common,samples,obs,cfg.EVENT_ORIGIN,settings,graphs,cfg)
    spatial_details={}
    if spatial_mode:
        from object_pipeline import execute as object_execute
        spatial_details=object_execute(cfg,obs,root,rasters,metadata)
        recovery_maps(root/"Recovery Map Raster File",graphs,cfg)
    else:
        recovery_maps(rasters,graphs,cfg)
    write_plots(common,obs,cfg.EVENT_ORIGIN,settings,summary,graphs,cfg)
    return dict(spatial_details,grid=grid,inside_pixels=summary.inside,independent_real_checks=summary.real_checks,
                eligible_model_count=valid_models,classified_nonzero_count=int(summary.counts["primary_unfilled"][1:].sum()))

def find_existing(cfg,ref,resources):
    candidates=[]
    for path in (cfg.PROJECT_ROOT/"outputs"/"Recovery Map Raster File").glob(
            "series_*/Infanta_Baseline_Recovery_Classes_*.tif"):
        state=cfg.PROJECT_ROOT/"outputs"/path.parent.name/"run_state.json"
        if ".partial." not in path.name and state.is_file() and json.loads(state.read_text())["status"]=="complete":
            candidates.append(path)
    candidates.sort()
    if not candidates:
        return None
    src=resources.enter_context(rasterio.open(candidates[-1],"r"))
    if src.crs!=ref.crs or src.transform!=ref.transform or src.shape!=ref.shape:
        logging.warning("Existing legacy raster grid differs; comparison omitted.")
        return None
    return src

def verify_decisions(output,ref,obs,cfg,settings,summary):
    """Reassess exported unfilled coherence to ensure files retain the decision inputs."""
    with ExitStack() as resources:
        sources={name:resources.enter_context(rasterio.open(path,"r")) for name,path in output.paths.items()}
        for window in windows(ref.width,ref.height,cfg.WINDOW_SIZE):
            reason=sources["exclusion_reason"].read(1,window=window)
            inside=reason!=1
            if not inside.any():
                continue
            data=sources["input_coherence"].read(window=window)
            data[data==cfg.NODATA]=np.nan
            expected=assess(data,obs,cfg.EVENT_ORIGIN,settings,inside)
            for name,array in {"ITSA_Recovery_Classes":expected.classes,"quality_flags":expected.quality,
                               "exclusion_reason":expected.exclusion}.items():
                if not np.array_equal(sources[name].read(1,window=window),array):
                    raise RuntimeError(f"Exported decisions differ: {name}")
            for name,a in expected.metrics.items():
                target=np.where(np.isfinite(a),a,cfg.NODATA)
                if not np.allclose(sources[name].read(1,window=window),target,atol=1e-6):
                    raise RuntimeError(f"Exported model indicator differs: {name}")


def find_previous_v2(cfg,ref,resources):
    candidates=list((cfg.PROJECT_ROOT/"Corrected_outputs").glob("itsa_v2_*/run_manifest.json"))
    candidates+=list((cfg.PROJECT_ROOT/"Output_Archive").glob("Corrected_outputs_*/itsa_v2_*/run_manifest.json"))
    for manifest_path in sorted(candidates,reverse=True):
        state=manifest_path.parent/"run_state.json"
        if not state.exists() or json.loads(state.read_text()).get("status")!="complete":
            continue
        path=manifest_path.parent/"Recovery Map Raster File"/"ITSA_Recovery_Classes.tif"
        src=resources.enter_context(rasterio.open(path,"r"))
        if src.crs==ref.crs and src.transform==ref.transform and src.shape==ref.shape:
            return src,json.loads(manifest_path.read_text())
    return None
