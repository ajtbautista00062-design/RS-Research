"""Bounded nearest-neighbor filling of small enclosed gaps before temporal analysis."""
import numpy as np
from scipy.ndimage import label, distance_transform_edt
from rasterio.windows import Window
import rasterio
from .spatial import read_stack, mask_for

def validate(c):
    if (type(c.FILL_MAX_HOLE_PIXELS) is not int or c.FILL_MAX_HOLE_PIXELS < 1
            or type(c.FILL_RADIUS_PIXELS) is not int or c.FILL_RADIUS_PIXELS < 1):
        raise ValueError("Fill hole area and distance limits must be positive integer pixels.")

def fill_holes(layer, c):
    """Original valid values are unchanged; edge-connected/large gaps remain missing."""
    result = layer.copy()
    finite = np.isfinite(layer)
    provenance = finite.astype("uint8")  # 0 missing, 1 valid input, 2 spatial estimate
    if not c.FILL_SMALL_HOLES or finite.all() or not finite.any():
        return result, provenance
    # 8-connected missing regions; preserve any region connected to the window edge.
    labels, _ = label(~finite, structure=np.ones((3,3)))
    sizes = np.bincount(labels.ravel())
    allowed = sizes <= c.FILL_MAX_HOLE_PIXELS
    allowed[0] = False
    edge = np.unique(np.concatenate((labels[0],labels[-1],labels[:,0],labels[:,-1])))
    allowed[edge] = False
    candidates = allowed[labels]
    if candidates.any():
        distance, indices = distance_transform_edt(~finite, return_indices=True)
        # Fill an eligible hole completely or leave it completely unfilled.
        too_far = np.unique(labels[candidates & (distance > c.FILL_RADIUS_PIXELS)])
        allowed[too_far] = False
        take = allowed[labels]
        result[take] = layer[indices[0][take],indices[1][take]]
        provenance[take] = 2
    return result, provenance

def corrected_window(vrts, srcs, ref, window, c, footprints, geometries):
    # Any <=max-area connected component extends at most max-area pixels.
    # This halo makes classification independent of processing tile edges.
    halo = c.FILL_MAX_HOLE_PIXELS+c.FILL_RADIUS_PIXELS+1 if c.FILL_SMALL_HOLES else 0
    row=max(0,int(window.row_off)-halo); col=max(0,int(window.col_off)-halo)
    bottom=min(ref.height,int(window.row_off+window.height)+halo)
    right=min(ref.width,int(window.col_off+window.width)+halo)
    expanded=Window(col,row,right-col,bottom-row)
    stack=read_stack(vrts,srcs,expanded,c,footprints)
    inside=mask_for(geometries,stack.shape[1:],
                    rasterio.windows.transform(expanded,ref.transform),c.ALL_TOUCHED)
    stack[:,~inside]=np.nan
    corrected=[]; provenance=[]
    for layer in stack:
        data, flag=fill_holes(layer,c)
        # Geometry mask is authoritative, including KML interior rings.
        data[~inside]=np.nan
        flag[~inside]=0
        corrected.append(data); provenance.append(flag)
    rr=int(window.row_off)-row; cc=int(window.col_off)-col
    cut=(slice(None),slice(rr,rr+int(window.height)),slice(cc,cc+int(window.width)))
    return np.stack(corrected)[cut],np.stack(provenance)[cut]
