"""Fixed PRE-only regions; identical members before and after the correction."""
from dataclasses import replace
import numpy as np
from .spatial_objects import cell_sizes, pre_change_only, segment, aggregate
from .recovery import assess, Assessment
from .baseline_recovery import assess as original_assess
from .temporal import day_offsets
from .itsa import median

# Upper bound on gathered values per vectorized median call (time x groups x size).
_GATHER_BUDGET = 1 << 24


def object_centers(stack, ids, keys, counts, c):
    """Per-object aggregate centers; equal to ``aggregate(stack, ids == gid, c)[0]``.

    Present member pixels are grouped ONCE by sorted object id instead of scanning
    the full tile mask for every id. Objects with the same pixel count are stacked
    into a compact (time, objects, size) array and reduced with one nanmedian.
    Members, finite counts, valid fraction, ``min(OBJECT_MIN_VALID_PIXELS, size)``
    and the median (valid zeros included, NaN ignored) are unchanged.
    """
    n_time = stack.shape[0]
    if not keys.size:
        return np.empty((n_time, 0))
    flat = ids.ravel()
    present = np.flatnonzero(flat > 0)
    order = np.argsort(flat[present], kind='stable')
    rows, cols = np.unravel_index(present[order], ids.shape)
    starts = np.concatenate(([0], np.cumsum(counts)[:-1]))
    out = None
    for size in np.unique(counts):
        size = int(size)
        groups = np.flatnonzero(counts == size)
        step = max(1, _GATHER_BUDGET // max(1, n_time * size))
        for first in range(0, groups.size, step):
            chunk = groups[first:first+step]
            take = starts[chunk][:, None] + np.arange(size)
            values = stack[:, rows[take], cols[take]]          # (time, objects, size)
            count = np.isfinite(values).sum(2)
            fraction = count / size
            center = median(values, axis=2)
            need = min(c.OBJECT_MIN_VALID_PIXELS, size)
            center[(count < need) | (fraction < c.OBJECT_MIN_VALID_FRACTION)] = np.nan
            if out is None:
                out = np.empty((n_time, keys.size), dtype=center.dtype)
            out[:, chunk] = center
    return out


def object_series(stack, inside, obs, c, area, row, col, reference_width):
    pre = np.array([o.phase == 'pre' for o in obs])
    cell,_ = cell_sizes(area,c)
    change = pre_change_only(stack[pre],day_offsets(obs,c.EVENT_ORIGIN)[pre])
    ids = segment(stack[pre],None,change,inside,np.ones(inside.shape,'uint8'),
                  cell,row,col,reference_width,c)
    keys,counts = np.unique(ids[ids>0],return_counts=True)
    series = object_centers(stack,ids,keys,counts,c) if keys.size else np.empty((len(obs),0))
    return ids,keys,counts,series


def project(batch, ids, keys, stack, inside, area, c):
    present = ids > 0
    index = np.searchsorted(keys,ids[present])
    def paint(a, fill, dtype):
        out = np.full(ids.shape,fill,dtype=dtype)
        out[present] = a[index]
        return out
    classes = paint(batch.classes,5,'uint8')
    absent = ~np.isfinite(stack).any(0)
    classes[absent | ~inside] = 0
    reasons = paint(batch.reasons,2,'uint8')
    reasons[absent] = 1
    reasons[~inside] = 11
    quality = paint(batch.quality,0,'uint16')
    counts = np.bincount(index,minlength=len(keys))
    if keys.size:
        small = counts*area < c.MIN_MAPPING_AREA_M2
        quality[present] |= np.where(small[index],4,0).astype('uint16')
    quality[~inside] = 65535
    metrics = {k:paint(v,np.nan,'float32') for k,v in batch.metrics.items()}
    gates = {k:paint(v,False,bool) for k,v in batch.gates.items()}
    return Assessment(classes,quality,reasons,metrics,None,gates)


def evaluate(stack,inside,obs,c,area,row,col,width):
    ids,keys,counts,series = object_series(stack,inside,obs,c,area,row,col,width)
    if not keys.size:
        # assess handles zero-column batches; avoid absent-pattern numerical work.
        empty = assess(stack,obs,c)
        old = original_assess(stack,obs,c)
        for result in (empty,old):
            result.reasons[~inside]=11
            result.quality[~inside]=65535
        return empty,old,ids,keys,counts,series
    current = assess(series,obs,c)
    original = original_assess(series,obs,c)
    return (project(current,ids,keys,stack,inside,area,c),
            project(original,ids,keys,stack,inside,area,c),ids,keys,counts,series)
