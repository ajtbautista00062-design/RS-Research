"""Read-only aligned inputs, explicit zero policy, optional separate spatial estimates."""
import numpy as np
import rasterio
from rasterio.windows import Window
from sarflow.spatial import aligned,mask_for
from sarflow.correction import fill_holes

def open_readers(resources,srcs,ref,c):
    from rasterio.enums import Resampling
    return ([resources.enter_context(aligned(src,ref,c)) for src in srcs],
            [resources.enter_context(aligned(src,ref,c,Resampling.nearest)) for src in srcs])

def read_stack(vrts,footprints,srcs,window,c):
    layers=[]
    for src,vrt,footprint in zip(srcs,vrts,footprints):
        a=vrt.read(c.BAND,window=window,masked=True).filled(np.nan).astype("float32")
        support=footprint.read(c.BAND,window=window,masked=True).filled(np.nan)
        invalid=~np.isfinite(a)|~np.isfinite(support)
        if c.ZERO_IS_INVALID:
            invalid|=(a==0)|(support==0)
        a[invalid]=np.nan
        a=(a*src.scales[c.BAND-1]+src.offsets[c.BAND-1])*c.COHERENCE_SCALE+c.COHERENCE_OFFSET
        a[~np.isfinite(a)]=np.nan
        finite=a[np.isfinite(a)]
        if finite.size and (finite.min()<-1e-6 or finite.max()>1+1e-6):
            raise ValueError(f"Coherence range invalid: {src.name}")
        layers.append(np.clip(a,0,1))
    return np.stack(layers)

def read_window(readers,srcs,ref,window,c,geometries):
    halo=c.FILL_MAX_HOLE_PIXELS+c.FILL_RADIUS_PIXELS+1 if c.FILL_SMALL_HOLES else 0
    r=max(0,int(window.row_off)-halo); col=max(0,int(window.col_off)-halo)
    endr=min(ref.height,int(window.row_off+window.height)+halo)
    endc=min(ref.width,int(window.col_off+window.width)+halo)
    expanded=Window(col,r,endc-col,endr-r)
    stack=read_stack(*readers,srcs,expanded,c)
    inside=mask_for(geometries,stack.shape[1:],
        rasterio.windows.transform(expanded,ref.transform),c.ALL_TOUCHED)
    stack[:,~inside]=np.nan
    flags=np.isfinite(stack).astype("uint8")
    if c.FILL_SMALL_HOLES:
        for j,a in enumerate(stack):
            stack[j],flags[j]=fill_holes(a,c)
        stack[:,~inside]=np.nan; flags[:,~inside]=0
    rr=int(window.row_off)-r; cc=int(window.col_off)-col
    cut=(slice(None),slice(rr,rr+int(window.height)),slice(cc,cc+int(window.width)))
    return stack[cut],flags[cut]
