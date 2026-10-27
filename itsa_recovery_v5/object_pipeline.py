"""Three bounded passes: PRE segmentation, per-object ITSA, class projection.

The existing pixel workflow supplies an aligned unfilled disk stack; no temporal
stack is loaded for the whole study area. Object IDs use global seed positions.
"""
from contextlib import ExitStack
from pathlib import Path
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from object_segmentation import pre_features,segment_cell,stable_ids
from object_statistics import classify_region
from spatial_settings import validate
from spatial_io import (active_cells,expanded,read_values,ancillary,strata_for,writer)


def pixel_scores(classes,guards,c):
    scale=max(c.RETURN_TOLERANCE,c.MIN_PROGRESS,c.LOSS_MIN)
    return np.where(np.isin(classes,[1,2,3,4]),np.clip(1-np.nan_to_num(guards,nan=np.inf)/scale,0,1),0)


def region_row(ids,key,values,features,strata,window,pixel_area,obs,c,settings,classes,scores):
    mask=ids==key; rr,cc=np.where(mask)
    row=classify_region(values[:,mask],obs,c.EVENT_ORIGIN,settings,c,
        int(mask.sum())*pixel_area,classes[mask],scores[mask])
    row.update(region_id=int(key),pixels=int(mask.sum()),area_m2=int(mask.sum())*pixel_area,
        stratum=int(strata[mask][0]),row_min=int(rr.min()+window.row_off),
        row_max=int(rr.max()+window.row_off+1),col_min=int(cc.min()+window.col_off),
        col_max=int(cc.max()+window.col_off+1),
        feature_min=np.min(features[:,mask],axis=1).tolist(),
        feature_max=np.max(features[:,mask],axis=1).tolist())
    if c.ENABLE_THRESHOLD_SENSITIVITY_REPORT:
        from copy import copy
        candidates=[(f'loss_{int(v*100):03d}','LOSS_MIN',v) for v in c.DETECTION_SWEEP]
        candidates += [('return_003','RETURN_TOLERANCE',.03),('return_007','RETURN_TOLERANCE',.07)]
        for name,field,threshold in candidates:
            variant=copy(settings);setattr(variant,field,threshold)
            candidate=classify_region(values[:,mask],obs,c.EVENT_ORIGIN,variant,c,
                row['area_m2'],classes[mask],scores[mask])
            row[f'sensitivity_{name}']=candidate['class_code']
    return row


def segment_and_fit(resources,stack,exclusion,pixels,guards,barrier,builtup,c,obs,folder,tags):
    geometry=validate(c,stack); settings=c.analysis_settings()
    id_writer=writer(resources,folder,'Spatial_Object_ID',stack,'uint32',0,tags)
    pre=[i for i,o in enumerate(obs) if o.phase=='pre']
    records={}
    for window in active_cells(exclusion,geometry['cell_side'],c.WINDOW_SIZE):
        area_window=expanded(window,stack,c.PRE_NEIGHBORHOOD_RADIUS)
        values=read_values(stack,area_window)
        strata=strata_for(area_window,values.shape[1:],barrier,builtup)
        inside=exclusion.read(1,window=area_window)!=1
        strata[~inside]=0
        features=pre_features(values[pre],strata,c.PRE_NEIGHBORHOOD_RADIUS)
        r=int(window.row_off-area_window.row_off); col=int(window.col_off-area_window.col_off)
        cut=(slice(r,r+int(window.height)),slice(col,col+int(window.width)))
        values=values[(slice(None),*cut)]; features=features[(slice(None),*cut)]
        strata=strata[cut]
        valid=(np.isfinite(values[pre]).sum(axis=0)>=settings.MIN_PRE)&(strata>0)
        labels=segment_cell(features,valid,strata,geometry['target_pixels'],geometry['maximum_pixels'],
                            geometry['minimum_pixels'],np.asarray(c.PRE_FEATURE_TOLERANCES))
        ids=stable_ids(labels,int(window.row_off),int(window.col_off),stack.width)
        id_writer.write(ids,1,window=window)
        cls=pixels.read(1,window=window)
        scores=pixel_scores(cls,read_values(guards,window)[0],settings)
        for key in np.unique(ids[ids>0]):
            records[int(key)]=region_row(ids,key,values,features,strata,window,
                geometry['pixel_area'],obs,c,settings,cls,scores)
            if len(records)>c.MAX_SPATIAL_OBJECTS:
                raise ValueError('MAX_SPATIAL_OBJECTS exceeded; review region features/sizes before rerunning.')
    id_writer.close()
    if not records:
        raise ValueError('No PRE-supported spatial objects; verify masks and input coverage.')
    return records,geometry


def execute(c,obs,root,pixel_folder,tags):
    """Called only by an explicit source run, never on import."""
    folder=root/'Recovery Map Raster File'; folder.mkdir(exist_ok=False)
    tags=dict(tags,implementation_version='5',classification_mode=c.CLASSIFICATION_MODE,
        segmentation='PRE-only bounded region growing; fixed global analytical cells',
        confidence='operational evidence score, not probability',
        builtup_mask_verified=bool(c.VERIFIED_BUILTUP_MASK),
        spatial_configuration=json.dumps({k:getattr(c,k) for k in (
            'MIN_REGION_AREA_M2','TARGET_REGION_AREA_M2','MAX_REGION_AREA_M2',
            'PRE_FEATURE_TOLERANCES','OBJECT_MIN_VALID_FRACTION','CLASS_CONFIDENCE_THRESHOLD',
            'SPATIAL_BARRIER_RASTER','VERIFIED_BUILTUP_MASK','BUILTUP_MASK_JUSTIFICATION')}))
    with ExitStack() as resources:
        stack=resources.enter_context(rasterio.open(pixel_folder/'input_coherence.tif','r'))
        reason=resources.enter_context(rasterio.open(pixel_folder/'exclusion_reason.tif','r'))
        pixels=resources.enter_context(rasterio.open(pixel_folder/'ITSA_Recovery_Classes.tif','r'))
        guards=resources.enter_context(rasterio.open(pixel_folder/'final_guard_halfwidth.tif','r'))
        barrier=ancillary(resources,c.SPATIAL_BARRIER_RASTER,stack)
        builtup=ancillary(resources,c.VERIFIED_BUILTUP_MASK,stack)
        records,geometry=segment_and_fit(resources,stack,reason,pixels,guards,barrier,builtup,c,obs,folder,tags)
        labels=resources.enter_context(rasterio.open(folder/'Spatial_Object_ID.tif','r'))
        from object_regularization import prepare_regularized
        regular,mapping=prepare_regularized(c,records,labels,stack,pixels,guards,obs,geometry)
        from object_products import export_products,export_tables
        export_products(resources,c,folder,labels,stack,reason,pixels,barrier,builtup,records,regular,mapping,tags)
        export_tables(root,obs,records,regular,geometry,c)
    from object_quality import verify
    verification=verify(folder,records,regular,c)
    (root/'Reports'/'spatial_verification.json').write_text(json.dumps(verification,indent=2),encoding='utf-8')
    from object_figures import example_trajectories
    example_trajectories(records,obs,root/'Time Series Graphs',c)
    return dict(spatial_objects=len(records),spatial_geometry=geometry,
                spatial_classification_mode=c.CLASSIFICATION_MODE,
                spatial_verification=verification,
                spatial_authoritative=False)
