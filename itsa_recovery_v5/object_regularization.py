"""Raster-backed region adjacency and independently refitted pair proposals."""
import numpy as np
from rasterio.windows import Window
from spatial_io import grid_windows,expanded,read_values
from spatial_regularization import boundary_pairs,regularize
from object_statistics import classify_region


def prepare_regularized(c,records,labels,stack,pixels,guards,obs,geometry):
    if c.CLASSIFICATION_MODE!='spatially_regularized_itsa' or not c.ENABLE_SPATIAL_REGULARIZATION:
        return None,None
    pairs=set()
    for window in grid_windows(labels,c.WINDOW_SIZE):
        pairs.update(boundary_pairs(labels.read(1,window=expanded(window,labels,1))))
        if len(pairs)>4*c.MAX_SPATIAL_OBJECTS:
            raise ValueError('Spatial adjacency memory cap exceeded; review segmentation.')
    def refit(keys):
        parents=[records[key] for key in keys]
        r=min(p['row_min'] for p in parents); col=min(p['col_min'] for p in parents)
        stop_r=max(p['row_max'] for p in parents); stop_c=max(p['col_max'] for p in parents)
        window=Window(col,r,stop_c-col,stop_r-r)
        # Disjoint adjacent original cell pairs have a bounded support rectangle.
        if window.width*window.height>4*geometry['cell_side']**2:
            raise ValueError('Unexpectedly unbounded regularization proposal.')
        mask=np.isin(labels.read(1,window=window),keys)
        values=read_values(stack,window)[:,mask]
        cls=pixels.read(1,window=window)[mask]
        from object_pipeline import pixel_scores
        scores=pixel_scores(cls,read_values(guards,window)[0][mask],c.analysis_settings())
        row=classify_region(values,obs,c.EVENT_ORIGIN,c.analysis_settings(),c,
            int(mask.sum())*geometry['pixel_area'],cls,scores)
        row.update(region_id=min(keys),pixels=int(mask.sum()),area_m2=int(mask.sum())*geometry['pixel_area'],
            stratum=parents[0]['stratum'],row_min=r,row_max=stop_r,col_min=col,col_max=stop_c,
            feature_min=np.min([p['feature_min'] for p in parents],axis=0).tolist(),
            feature_max=np.max([p['feature_max'] for p in parents],axis=0).tolist(),
            merged_from=','.join(map(str,keys)))
        return row
    return regularize(records,pairs,c,refit)
