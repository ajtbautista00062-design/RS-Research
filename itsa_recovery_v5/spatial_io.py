"""Bounded spatial reads and GeoTIFF products. Original rasters are read-only."""
from pathlib import Path
import numpy as np
import rasterio
from rasterio.windows import Window
from recovery import COLORS
from raster_outputs import qml


def grid_windows(ref,side):
    for r in range(0,ref.height,side):
        for c in range(0,ref.width,side):
            yield Window(c,r,min(side,ref.width-c),min(side,ref.height-r))


def expanded(window,ref,halo):
    r=max(0,int(window.row_off)-halo); c=max(0,int(window.col_off)-halo)
    stop_r=min(ref.height,int(window.row_off+window.height)+halo)
    stop_c=min(ref.width,int(window.col_off+window.width)+halo)
    return Window(c,r,stop_c-c,stop_r-r)


def active_cells(exclusion,side,io_size):
    """Read only the 1-band boundary mask to schedule globally anchored cells."""
    cells=set()
    for window in grid_windows(exclusion,io_size):
        rr,cc=np.where(exclusion.read(1,window=window)!=1)
        if len(rr):
            keys=np.column_stack(((rr+int(window.row_off))//side,(cc+int(window.col_off))//side))
            cells.update(map(tuple,np.unique(keys,axis=0).tolist()))
    for r,c in sorted(cells):
        yield Window(c*side,r*side,min(side,exclusion.width-c*side),min(side,exclusion.height-r*side))


def read_values(src,window):
    return src.read(window=window,masked=True).astype('float32').filled(np.nan)


def ancillary(resources,path,ref):
    if not path:
        return None
    source=Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    src=resources.enter_context(rasterio.open(source,'r'))
    if src.crs!=ref.crs or src.shape!=ref.shape or src.transform!=ref.transform or src.count!=1:
        raise ValueError('Ancillary masks/barriers require one band on the exact analysis grid.')
    return src


def strata_for(window,shape,barrier,builtup):
    result=np.ones(shape,dtype='uint32')
    if barrier:
        a=barrier.read(1,window=window,masked=True).astype('float64').filled(np.nan)
        good=np.isfinite(a)&(a>0)
        if np.any(good&((a!=np.floor(a))|(a>=np.iinfo('uint32').max))):
            raise ValueError('Barrier raster must contain positive integer category IDs.')
        result=np.where(good,a,0).astype('uint32')
    if builtup:
        a=read_values(builtup,window)[0]
        result[~np.isfinite(a)|(a<=0)]=0
    return result


def writer(resources,folder,name,ref,dtype,nodata,tags):
    path=folder/f'{name}.tif'
    dst=resources.enter_context(rasterio.open(path,'w',driver='GTiff',width=ref.width,
        height=ref.height,crs=ref.crs,transform=ref.transform,count=1,dtype=dtype,nodata=nodata,
        tiled=True,blockxsize=256,blockysize=256,compress='deflate',BIGTIFF='IF_SAFER',
        predictor=3 if dtype=='float32' else 2))
    dst.set_band_description(1,name)
    dst.update_tags(**tags,product=name)
    if name.endswith('Classes'):
        dst.write_colormap(1,COLORS); qml(path)
    return dst


def project_field(ids,records,field,missing):
    """No full-grid lookup indexed by potentially sparse global object IDs."""
    keys,inverse=np.unique(ids,return_inverse=True)
    values=np.array([records[int(key)].get(field,missing) if key else missing for key in keys])
    return values[inverse].reshape(ids.shape)
