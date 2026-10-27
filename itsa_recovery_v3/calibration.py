"""PRE-only data-informed noise calibration, independent of recovery class frequencies."""
from contextlib import ExitStack
import numpy as np
import rasterio
from sarflow.spatial import read_zones,mask_for,windows
from sarflow.sampling import Sample
from preparation import open_readers,read_window
from model import design

def noise_proxy(pre,time):
    y=pre.reshape(len(pre),-1).astype(float)
    valid=np.isfinite(y).all(axis=0)
    result=np.full(y.shape[1],np.nan)
    x=np.column_stack((np.ones(len(pre)),time))
    if valid.any() and len(pre)>2 and np.linalg.matrix_rank(x)==2:
        response=y[:,valid]
        residual=response-x@np.linalg.pinv(x)@response
        result[valid]=np.sqrt(np.sum(residual**2,axis=0)/(len(pre)-2))
    return result.reshape(pre.shape[1:])

def from_sample(sample,error_floor):
    if not sample.count:
        raise ValueError("No finite nonzero PRE controls available for noise calibration.")
    value=max(error_floor,float(np.median(sample.values[0])))
    return dict(operational_floor=value,method="median PRE-only detrended temporal noise proxy",
        sample_pixels=sample.values.shape[1],population_pixels=sample.count,validated=False,
        independent_ground_truth_available=False,
        limitation="PRE residuals may contain real environmental variability; three PRE dates give df=1")

def calibrate(cfg,obs):
    preobs=[o for o in obs if o.phase=="pre"]
    _,_,t=design(preobs,cfg.EVENT_ORIGIN)
    sample=Sample(cfg.CALIBRATION_SAMPLE_PIXELS,cfg.RANDOM_SEED+1,1)
    with ExitStack() as resources:
        srcs=[resources.enter_context(rasterio.open(o.path,"r")) for o in preobs]
        event=next(o for o in obs if o.phase=="event")
        ref=resources.enter_context(rasterio.open(cfg.REFERENCE_RASTER or event.path,"r"))
        geoms=[g for _,g in read_zones(cfg.Study_Area,ref.crs,cfg.PROJECT_ROOT)]
        settings=cfg.reader_settings(exclude_zero=cfg.ZERO_SEMANTICS=="background")
        readers=open_readers(resources,srcs,ref,settings)
        for window in windows(ref.width,ref.height,cfg.WINDOW_SIZE):
            inside=mask_for(geoms,(int(window.height),int(window.width)),
                rasterio.windows.transform(window,ref.transform),cfg.ALL_TOUCHED)
            if not inside.any():
                continue
            values,_=read_window(readers,srcs,ref,window,settings,geoms)
            sigma=noise_proxy(values,t)
            valid=inside&np.isfinite(sigma)&np.any(values>0,axis=0)
            rr,cc=np.where(valid)
            if len(rr):
                sample.add(sigma[None,valid],rr+int(window.row_off),cc+int(window.col_off))
    return from_sample(sample,cfg.analysis_settings().ERROR_FLOOR)
