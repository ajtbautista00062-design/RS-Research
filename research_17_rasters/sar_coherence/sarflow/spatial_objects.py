"""PRE-only connected regions in fixed physical domains, independent of I/O tiles.

Hard ancillary labels are barriers. Fixed domains bound RAM and object extent;
they are an explicit segmentation constraint, not inferred landscape boundaries.
"""
import math
import numpy as np
from scipy.ndimage import label
from .itsa import median
from .recovery import assess, Assessment

CROSS = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])


def cell_sizes(pixel_area, c):
    target = max(1, math.floor(math.sqrt(c.TARGET_REGION_AREA_M2 / pixel_area)))
    domain = max(target, math.floor(math.sqrt(c.MAX_REGION_AREA_M2 / pixel_area)))
    domain = max(target, (domain // target) * target)
    return target, domain


def pre_change_only(pre, days):
    """Segmentation feature independent of event/POST availability and values."""
    y = pre.reshape(len(pre), -1)
    change = np.full(y.shape[1], np.nan)
    patterns, groups = np.unique(np.isfinite(y).T, axis=0, return_inverse=True)
    for j, pattern in enumerate(patterns):
        if pattern.sum() < 3:
            continue
        ix = np.flatnonzero(groups == j)
        x = np.column_stack((np.ones(pattern.sum()), days[pattern]-np.mean(days[pattern])))
        beta = np.linalg.lstsq(x, y[pattern][:, ix], rcond=None)[0]
        change[ix] = beta[1] * np.ptp(days[pattern])
    return change.reshape(pre.shape[1:])


def segment(pre, noise, pre_change, mask, barriers, cell, row, col, width, c):
    """Four-connected PRE-feature bins; globally unique IDs = minimum cell index+1.

    Object boundaries use only PRE features, external barriers and fixed spatial
    cells. No POST values or recovery labels enter segmentation.
    """
    baseline = median(pre)
    spread = 1.4826 * median(np.abs(pre - baseline))
    features = np.stack((baseline, spread, pre_change))
    valid = (mask & (barriers > 0) & (np.isfinite(pre).sum(0) >= c.ITSA_MIN_PRE)
             & np.isfinite(features).all(0))
    bins = np.floor(np.nan_to_num(features) /
                    np.asarray(c.PRE_FEATURE_BIN_WIDTHS)[:, None, None] + 1e-8).astype("int32")
    ids = np.zeros(mask.shape, dtype="uint32")
    height, wide = mask.shape
    for r in range(0, height, cell):
        for s in range(0, wide, cell):
            sl = np.s_[r:min(r+cell, height), s:min(s+cell, wide)]
            ok = valid[sl]
            if not ok.any():
                continue
            keys = np.column_stack((bins[:, sl[0], sl[1]][:, ok].T, barriers[sl][ok]))
            groups = np.unique(keys, axis=0)
            target = ids[sl]
            for key in groups:
                member = ok & (bins[:, sl[0], sl[1]] == key[:3, None, None]).all(0)
                member &= barriers[sl] == key[3]
                parts, number = label(member, CROSS)
                for part in range(1, number+1):
                    rr, cc = np.where(parts == part)
                    gid = int(np.min((row+r+rr) * width + col+s+cc)) + 1
                    if gid > np.iinfo("uint32").max:
                        raise ValueError("Analysis grid exceeds uint32 region-ID capacity.")
                    target[rr, cc] = gid
    return ids


def aggregate(stack, member, c):
    values = stack[:, member]
    count = np.isfinite(values).sum(1)
    fraction = count / member.sum()
    center = median(values, axis=1)
    spread = 1.4826 * median(np.abs(values-center[:, None]), axis=1)
    need = min(c.OBJECT_MIN_VALID_PIXELS, int(member.sum()))
    center[(count < need) | (fraction < c.OBJECT_MIN_VALID_FRACTION)] = np.nan
    return center, count, fraction, spread


def region_assessments(stack, ids, obs, c, pixel_area):
    records = {}
    for gid in np.unique(ids[ids > 0]):
        member = ids == gid
        series, counts, fractions, spread = aggregate(stack, member, c)
        records[int(gid)] = (None, series, counts, fractions, spread)
    if records:
        batch = assess(np.column_stack([v[1] for v in records.values()]), obs, c)
        for index, gid in enumerate(records):
            result = Assessment(batch.classes[index:index+1].copy(), batch.quality[index:index+1].copy(),
                                batch.reasons[index:index+1].copy(),
                                {k: v[index:index+1].copy() for k, v in batch.metrics.items()}, None,
                                {k: v[index:index+1].copy() for k, v in batch.gates.items()})
            if (ids == gid).sum() * pixel_area < c.MIN_REGION_AREA_M2:
                result.quality[:] |= 4
            records[gid] = (result, *records[gid][1:])
    return records


def neighbors(ids):
    pairs = set()
    for a, b in ((ids[:-1], ids[1:]), (ids[:, :-1], ids[:, 1:])):
        use = (a > 0) & (b > 0) & (a != b)
        for left, right in np.unique(np.column_stack((a[use], b[use])), axis=0):
            pairs.add(tuple(sorted((int(left), int(right)))))
    return sorted(pairs)


def regularize(stack, ids, records, barriers, obs, c, pixel_area):
    """Merge ONLY same-class low-confidence neighboring objects, then refit.

    No voting, class promotion, inferred confidence from pixel count, or NoData
    filling. The categorical map may be identical; object partitions can improve.
    A proposed merge must retain the class after exact pooled-pixel aggregation.
    """
    result_ids, result_records = ids.copy(), dict(records)
    edges = neighbors(ids)
    degree = {}
    for left, right in edges:
        degree[left] = degree.get(left, 0)+1
        degree[right] = degree.get(right, 0)+1
    for left, right in edges:
        if left not in result_records or right not in result_records:
            continue
        if max(degree[left], degree[right]) < c.REGULARIZATION_MIN_NEIGHBORS:
            continue
        a, b = result_records[left][0], result_records[right][0]
        code = int(a.classes[0])
        if code not in (1, 2, 3, 4) or code != int(b.classes[0]):
            continue
        # Retain small high-confidence regions as independent analytical units.
        if not (int(a.quality[0]) & 3 and int(b.quality[0]) & 3):
            continue
        member = (result_ids == left) | (result_ids == right)
        if member.sum() * pixel_area > c.MAX_REGION_AREA_M2:
            continue
        if len(np.unique(barriers[member])) != 1:
            continue
        pre = np.array([o.phase == "pre" for o in obs])
        left_pre, right_pre = result_records[left][1][pre], result_records[right][1][pre]
        common = np.isfinite(left_pre) & np.isfinite(right_pre)
        if (common.sum() < 3 or np.max(np.abs(left_pre[common]-right_pre[common])) >
                c.REGULARIZATION_FEATURE_TOLERANCE):
            continue
        keys = ("baseline", "modeled_level_loss", "modeled_progress", "baseline_deficit")
        if any(abs(float(a.metrics[key][0])-float(b.metrics[key][0])) >
               c.REGULARIZATION_FEATURE_TOLERANCE for key in keys):
            continue
        series, counts, fractions, spread = aggregate(stack, member, c)
        pooled = assess(series[:, None], obs, c)
        if int(pooled.classes[0]) != code:
            continue
        if member.sum() * pixel_area < c.MIN_REGION_AREA_M2:
            pooled.quality[:] |= 4
        result_ids[member] = left
        result_records[left] = (pooled, series, counts, fractions, spread)
        del result_records[right]
    return result_ids, result_records


def paint(ids, records, key, dtype="float32", nodata=np.nan):
    result = np.full(ids.shape, nodata, dtype=dtype)
    # Sort/map once instead of scanning the raster once per object.
    if not records:
        return result
    keys = np.array(sorted(records), dtype="uint32")
    def value(gid):
        fit = records[int(gid)][0]
        if key.startswith("gate:"):
            return fit.gates[key[5:]][0]
        return getattr(fit, key)[0] if key in ("classes", "quality", "reasons") else fit.metrics[key][0]
    values = np.array([value(gid) for gid in keys], dtype=dtype)
    use = ids > 0
    result[use] = values[np.searchsorted(keys, ids[use])]
    return result
