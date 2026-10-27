"""PRE-only compact region growing on globally anchored, bounded analytical cells.

Cell boundaries are part of the design, not arbitrary raster I/O tile boundaries.
No POST values, fitted outcomes or final classes enter segmentation.
"""
from collections import deque
import warnings
import numpy as np


def pre_features(pre, strata, radius=1):
    """Masked local PRE median, temporal MAD, local spatial MAD (coherence units)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        baseline = np.nanmedian(pre, axis=0)
        temporal = 1.4826*np.nanmedian(np.abs(pre-baseline), axis=0)
        h, w = strata.shape
        pad = radius
        values = np.pad(baseline, pad, constant_values=np.nan)
        groups = np.pad(strata, pad, constant_values=0)
        neighbors = []
        for dr in range(-pad, pad+1):
            for dc in range(-pad, pad+1):
                v = values[pad+dr:pad+dr+h, pad+dc:pad+dc+w]
                g = groups[pad+dr:pad+dr+h, pad+dc:pad+dc+w]
                neighbors.append(np.where((g == strata) & (strata > 0), v, np.nan))
        local = np.nanmedian(neighbors, axis=0)
        spatial = 1.4826*np.nanmedian(np.abs(np.asarray(neighbors)-local), axis=0)
    return np.stack((local, temporal, spatial))


def segment_cell(features, valid, strata, target_pixels, maximum_pixels,
                 minimum_pixels, tolerances):
    """Four-connected growing, bounded feature range; small compatible regions merge.

    Returns local labels. Small unmergeable objects survive for explicit MMU flags;
    no morphological filling, deletion, or class-driven merging is performed.
    """
    valid = valid & (strata > 0) & np.isfinite(features).all(axis=0)
    labels = np.zeros(valid.shape, dtype="uint32")
    h, w = valid.shape
    spacing = max(1, int(np.sqrt(target_pixels)))
    seeds = [(r, c) for r in range(spacing//2, h, spacing)
             for c in range(spacing//2, w, spacing)]
    seeds += list(zip(*np.where(valid)))
    regions = {}; next_id = 0
    def adjacent(r, c):
        return ((r-1,c),(r+1,c),(r,c-1),(r,c+1))
    for sr, sc in seeds:
        if not valid[sr, sc] or labels[sr, sc]:
            continue
        next_id += 1
        group = strata[sr, sc]
        low = features[:,sr,sc].copy(); high = low.copy()
        pending = deque([(sr,sc)]); members = []
        while pending and len(members) < target_pixels:
            r, c = pending.popleft()
            if not (0 <= r < h and 0 <= c < w) or not valid[r,c] or labels[r,c] or strata[r,c] != group:
                continue
            new_low = np.minimum(low, features[:,r,c]); new_high = np.maximum(high, features[:,r,c])
            if np.any(new_high-new_low > tolerances):
                continue
            low, high = new_low, new_high
            labels[r,c] = next_id; members.append((r,c))
            pending.extend(adjacent(r,c))
        regions[next_id] = dict(members=members,low=low,high=high,group=group)
    for source in sorted(regions, key=lambda k:len(regions[k]["members"])):
        region = regions[source]
        if not region["members"] or len(region["members"]) >= minimum_pixels:
            continue
        candidates = set()
        for r,c in region["members"]:
            for rr,cc in adjacent(r,c):
                if 0 <= rr < h and 0 <= cc < w and labels[rr,cc] not in (0,source):
                    candidates.add(int(labels[rr,cc]))
        choices = []
        for target in candidates:
            other = regions[target]
            low=np.minimum(region["low"],other["low"]); high=np.maximum(region["high"],other["high"])
            if (region["group"]==other["group"] and
                    len(region["members"])+len(other["members"]) <= maximum_pixels and
                    np.all(high-low <= tolerances)):
                choices.append((float(np.max((high-low)/tolerances)),target,low,high))
        if choices:
            _,target,low,high = min(choices,key=lambda x:(x[0],x[1]))
            other=regions[target]
            for r,c in region["members"]:
                labels[r,c]=target
            other["members"].extend(region["members"])
            other["low"],other["high"]=low,high
            region["members"]=[]
    return labels


def stable_ids(labels, row_offset, col_offset, raster_width):
    """Unique seed-cell IDs tied to global pixel positions, independent of I/O tiles."""
    result=np.zeros_like(labels,dtype="uint32")
    for label in np.unique(labels[labels>0]):
        rr,cc=np.where(labels==label)
        value=(int(rr[0])+row_offset)*raster_width+int(cc[0])+col_offset+1
        if value >= np.iinfo("uint32").max:
            raise ValueError("Grid exceeds uint32 object identifier range.")
        result[rr,cc]=value
    return result
