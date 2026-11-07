"""Multi-scale zonal ITSA/RRI recovery classification.

Why this module exists
----------------------
Single 14 m pixels of 12-day C-band coherence over vegetated Infanta are noisy
(median PRE residual noise ~0.03 with a median baseline ~0.15). Pixel decisions
are therefore speckled and most remain Indeterminate. This module classifies
*zones* instead of pixels:

1. Zones are square grid cells aligned to the reference grid, optionally split by
   barangay so a zone never mixes two barangays. The base cell is
   ``ZONAL_BASE_CELLS`` pixels (5 x 14 m = 70 m) and each further level doubles the
   side (70, 140, 280, 560 m). Levels nest exactly.
2. Each zone's per-date coherence is the median of its valid member pixels
   (multi-looking reduces estimator variance by roughly sqrt(n)).
3. The zone series is fitted with the same segmented ITSA OLS (``fit_itsa``) and
   classified with an *observation-anchored* decision tree (below).
4. Every pixel takes the class of the FINEST zone that resolves it. Pixels whose
   base zone stays Indeterminate are passed to the parent zone, and so on.
   The level used is exported for every pixel (``Zonal_Zone_Size_m``).
5. Outputs use the standard run layout: zonal ITSA files in
   ``B. ITSA/{Overall,Per Barangay}/{Raster Files,ITSA Graph,CSV}`` and zonal RRI
   files in ``A. RRI/{Raster Files,Graphs,CSV}/{Overall,Per Barangay}``; every file
   name starts with ``Zonal_``.

Decision tree (per zone; B = PRE median, E = event-straddling pair, H = adaptive
loss threshold, tol = RECOVERY_TOLERANCE, k = LATE_POST_CONFIRMATIONS)
- Indeterminate when: insufficient PRE/POST/event support, unstable fit, or the
  signal observability gate fails (B must exceed tol + H + margin and 3 x noise).
- Disturbed when the observed event-pair loss B - E >= H, or the modeled level
  loss -beta2 >= H confirmed by an observed early-POST deficit >= H.
- Disturbed zones: late = median of last k POST pairs, trough = min(E, first POST)
    late >= B - tol                         -> 1 Coherence recovered
    late - trough >= MIN_RECOVERY_IMPROVEMENT -> 2 Partially recovered
    otherwise                               -> 3 Not recovered / further degraded
- Not disturbed but early AND late POST deficits persist (early >= H, late > tol)
                                            -> 3 (post-event decline, evidence 3)
- Otherwise                                 -> 4 No detectable event-associated loss

Differences from the original pixel/object rules are deliberate and documented:
no indecision margins between adjacent classes, the observed event pair counts
as disturbance evidence (the original treated an event dip followed by a
post-event gain as "conflicting"), and the end state is the observed late median
(the modeled end state is kept as a quality check, flag 32).

Scientific limitations (exported in tags/CSV): event-pair decorrelation can be
caused by rain/wind/soil moisture rather than damage (flag 1 marks event-only
evidence); coherence magnitude bias at low true coherence is not corrected;
only 3 PRE and 3 POST pairs; zone medians mix land covers within a zone; no field
validation. Outputs are an exploratory zonal product, not accuracy-assessed.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.features import shapes
from rasterio.warp import transform_geom
from rasterio.windows import Window, transform as window_transform

from .itsa import fit_itsa, median
from .recovery import baseline_eligibility
from .rri import calculate, INTERPRETATION, STATUS
from .study_inputs import observations, effective_settings

CLASSES = {0: "NoData", 1: "Coherence recovered", 2: "Partially recovered",
           3: "Not recovered / further degraded", 4: "No detectable event-associated loss",
           5: "Indeterminate"}
REASONS = {0: "classified", 1: "no valid observations", 2: "insufficient PRE/POST/event support",
           3: "unstable or out-of-range model", 4: "baseline fails signal observability gate",
           12: "insufficient valid zone pixels at every zone size"}
EVIDENCE = {0: "not assessed", 1: "observed event-pair loss", 2: "modeled loss confirmed by early POST",
            3: "post-event decline without event-pair loss", 4: "no loss evidence"}
QUALITY = {1: "event-pair-only disturbance (weather-sensitive; early POST already near baseline)",
           2: "baseline below historical 0.40 observability gate",
           4: "resolved at a zone larger than the base zone",
           8: "late coherence exceeds baseline + tolerance",
           16: "source validity semantics unconfirmed",
           32: "modeled end state disagrees with observed late median"}

COLORS = {
    "class": {0: ("#000000", 0), 1: ("#1a9641", 255), 2: ("#a6d96a", 255), 3: ("#d7191c", 255),
              4: ("#2c7bb6", 255), 5: ("#bdbdbd", 255)},
    "rri": {0: ("#000000", 0), 1: ("#d7191c", 255), 2: ("#fdae61", 255), 3: ("#ffffbf", 255),
            4: ("#a6d96a", 255), 5: ("#1a9641", 255)},
    "status": {0: ("#000000", 0), 1: ("#1a9641", 255), 2: ("#2c7bb6", 255), 3: ("#fdae61", 255),
               4: ("#f46d43", 255), 5: ("#a50026", 255), 6: ("#bdbdbd", 255), 255: ("#000000", 0)},
    "evidence": {0: ("#bdbdbd", 255), 1: ("#d7191c", 255), 2: ("#fdae61", 255), 3: ("#7b3294", 255),
                 4: ("#2c7bb6", 255), 255: ("#000000", 0)},
    "reason": {0: ("#1a9641", 255), 1: ("#000000", 0), 2: ("#fdae61", 255), 3: ("#a50026", 255),
               4: ("#7b3294", 255), 12: ("#f4a582", 255), 255: ("#000000", 0)},
}
LABELS = {"class": CLASSES, "rri": {0: "NoData / RRI not valid", **{k: v for k, v in INTERPRETATION.items() if k}},
          "status": {**STATUS, 255: "NoData"}, "evidence": {**EVIDENCE, 255: "NoData"},
          "reason": {**REASONS, 255: "NoData"}}


def setting(c, name, default):
    return getattr(c, name, default)


# ---------------------------------------------------------------------------
# Classification of zone series
# ---------------------------------------------------------------------------
@dataclass
class ZoneResult:
    classes: np.ndarray
    reasons: np.ndarray
    evidence: np.ndarray
    quality: np.ndarray
    metrics: dict


def classify_series(y, obs, c):
    """Classify (time, n) zone series. Returns per-zone arrays."""
    y = np.asarray(y, dtype=float)
    n = y.shape[1]
    phase = np.array([o.phase for o in obs])
    pre, post = phase == "pre", phase == "post"
    event = np.flatnonzero(phase == "event")
    if len(event) != 1 or post.sum() < 3:
        raise ValueError("Zonal ITSA requires one event pair and at least three POST pairs.")
    fit = fit_itsa(y, obs, c)
    B, b = fit.baseline, fit.beta
    margin = c.EFFECT_DECISION_MARGIN
    tol = c.RECOVERY_TOLERANCE
    k = c.LATE_POST_CONFIRMATIONS
    noise = np.where(np.isfinite(fit.noise), fit.noise, np.nan)
    H = np.maximum(c.EVENT_LOSS_MINIMUM, c.PRE_RESIDUAL_NOISE_MULTIPLIER * np.fmax(noise, c.ITSA_NOISE_FLOOR))
    p = y[post]
    E = y[event[0]]
    early, late = median(p[:k]), median(p[-k:])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        trough = np.fmin(E, p[0])
    slope = b[1] + b[3]
    from .temporal import day_offsets
    days = day_offsets(obs, c.EVENT_ORIGIN)
    modeled_final = b[0] + b[2] + slope * days[post].max()
    predictions = b[0] + b[2] + slope * days[post, None]

    support = ((fit.pre_count >= c.ITSA_MIN_PRE) & (fit.post_count >= c.ITSA_MIN_POST)
               & np.isfinite(E) & np.isfinite(p[-k:]).all(0))
    stable = (np.isfinite(b).all(0) & (fit.rmse <= c.MAX_ITSA_RMSE)
              & ((predictions >= -0.05) & (predictions <= 1.05)).all(0))
    observable, absolute, _, _ = baseline_eligibility(B, fit, H, tol, c)
    usable = support & stable & observable

    event_loss = B - E
    modeled_loss = -b[2]
    early_deficit = B - early
    late_deficit = B - late
    disturbed_event = usable & (event_loss >= H)
    disturbed_model = usable & ~disturbed_event & (modeled_loss >= H) & (early_deficit >= H)
    disturbed = disturbed_event | disturbed_model
    recovered = disturbed & (late >= B - tol)
    partial = disturbed & ~recovered & (late - trough >= c.MIN_RECOVERY_IMPROVEMENT)
    not_recovered = disturbed & ~recovered & ~partial
    decline = usable & ~disturbed & (early_deficit >= H) & (late_deficit > tol)
    no_loss = usable & ~disturbed & ~decline

    classes = np.full(n, 5, "uint8")
    reasons = np.zeros(n, "uint8")
    evidence = np.zeros(n, "uint8")
    for mask, code in ((recovered, 1), (partial, 2), (not_recovered, 3), (decline, 3), (no_loss, 4)):
        classes[mask] = code
    evidence[disturbed_event], evidence[disturbed_model] = 1, 2
    evidence[decline], evidence[no_loss] = 3, 4
    reasons[~support] = 2
    reasons[support & ~stable] = 3
    reasons[support & stable & ~observable] = 4
    absent = ~np.isfinite(y).any(0)
    classes[absent], reasons[absent] = 0, 1

    quality = np.zeros(n, "uint16")
    quality[disturbed_event & (early_deficit < tol)] |= 1
    quality[np.isfinite(B) & ~absolute] |= 2
    quality[usable & (late > B + tol)] |= 8
    if not c.SOURCE_VALIDITY_CONFIRMED:
        quality[:] |= 16
    quality[usable & (np.abs(modeled_final - late) > c.OBSERVATION_DISAGREEMENT_TOLERANCE)] |= 32

    value, rri_class, rri_status, _, _ = calculate(
        B, trough, late, disturbed, support, stable, np.ones(n, bool), H, no_loss=classes == 4)
    rri_status[absent] = 0
    metrics = dict(baseline=B, event_value=E, early_post=early, late_post=late, trough=trough,
                   observed_event_loss=event_loss, modeled_level_loss=modeled_loss,
                   early_deficit=early_deficit, late_deficit=late_deficit, loss_threshold=H,
                   modeled_final=modeled_final, rmse=fit.rmse, pre_noise=fit.noise,
                   recovery_index=value, rri_class=rri_class, rri_status=rri_status)
    return ZoneResult(classes, reasons, evidence, quality, metrics)


# ---------------------------------------------------------------------------
# Multi-scale zoning
# ---------------------------------------------------------------------------
@dataclass
class ZonalProduct:
    classes: np.ndarray
    reasons: np.ndarray
    evidence: np.ndarray
    quality: np.ndarray
    zone_id: np.ndarray
    zone_size_m: np.ndarray
    rri: np.ndarray
    rri_class: np.ndarray
    rri_status: np.ndarray
    zones: pd.DataFrame
    levels: pd.DataFrame


def grouped_median(values, inverse, groups):
    frame = pd.DataFrame(values.T)
    frame["_g"] = inverse
    grouped = frame.groupby("_g", sort=True)
    med = np.array(grouped.median().reindex(range(groups)).to_numpy(dtype=float).T)
    cnt = np.array(grouped.count().reindex(range(groups)).fillna(0).to_numpy(dtype=float).T)
    return med, cnt


def zonal_classification(stack, inside, owner, obs, c, pixel_size, base_cells=None, levels=None):
    """Assign every inside pixel the class of its finest resolving zone."""
    base = int(base_cells or setting(c, "ZONAL_BASE_CELLS", 5))
    nlev = int(levels or setting(c, "ZONAL_LEVELS", 4))
    min_pixels = int(setting(c, "ZONAL_MIN_VALID_PIXELS", 9))
    min_fraction = float(setting(c, "ZONAL_MIN_VALID_FRACTION", 0.5))
    split = bool(setting(c, "ZONAL_RESPECT_BARANGAYS", True))
    T, Hh, Ww = stack.shape
    rows, cols = np.nonzero(inside)
    vals = stack[:, rows, cols].astype("float64")
    own = owner[rows, cols].astype("int64") if split else np.zeros(rows.size, "int64")
    npx = rows.size
    out = dict(classes=np.full(npx, 5, "uint8"), reasons=np.full(npx, 12, "uint8"),
               evidence=np.zeros(npx, "uint8"), quality=np.zeros(npx, "uint16"),
               zone=np.zeros(npx, "int64"), size=np.zeros(npx, "uint16"),
               rri=np.full(npx, np.nan, "float32"), rri_class=np.zeros(npx, "uint8"),
               rri_status=np.full(npx, 6, "uint8"))
    resolved = np.zeros(npx, bool)
    tables, level_rows = [], []
    for level in range(nlev):
        b = base * 2 ** level
        last = level == nlev - 1
        nbx = -(-Ww // b)
        key = ((rows // b) * nbx + (cols // b)) * (int(own.max()) + 1) + own
        uniq, inverse = np.unique(key, return_inverse=True)
        size = np.bincount(inverse, minlength=uniq.size)
        med, cnt = grouped_median(vals, inverse, uniq.size)
        need = np.minimum(min_pixels, size) if last else np.full(uniq.size, min_pixels)
        ok = (cnt >= need) & (cnt >= min_fraction * size)
        med[~ok] = np.nan
        result = classify_series(med, obs, c)
        has_obs = np.isfinite(med).any(0)
        zone_class = np.where(has_obs, result.classes, 5)
        zone_reason = np.where(has_obs, result.reasons, 12).astype("uint8")
        zone_ids = (level + 1) * 10_000_000 + np.arange(uniq.size)
        decisive = np.isin(zone_class, (1, 2, 3, 4))
        take = ~resolved & (decisive[inverse] | last)
        z = inverse[take]
        out["classes"][take] = zone_class[z]
        out["reasons"][take] = zone_reason[z]
        out["evidence"][take] = result.evidence[z]
        out["quality"][take] = result.quality[z] | (4 if level else 0)
        out["zone"][take] = zone_ids[z]
        out["size"][take] = int(round(b * pixel_size))
        out["rri"][take] = result.metrics["recovery_index"][z]
        out["rri_class"][take] = result.metrics["rri_class"][z]
        out["rri_status"][take] = np.where(has_obs[z], result.metrics["rri_status"][z], 3)
        resolved |= take & decisive[inverse]
        assigned = np.bincount(z, minlength=uniq.size)
        level_rows.append(dict(level=level, zone_side_m=round(b * pixel_size, 1), zone_cells=b,
                               zones=int(uniq.size), pixels_assigned=int(take.sum()),
                               pixels_resolved=int((take & decisive[inverse]).sum())))
        used = assigned > 0
        frame = pd.DataFrame(dict(zone_id=zone_ids[used], level=level, zone_side_m=round(b * pixel_size, 1),
                                  barangay_index=(uniq % (int(own.max()) + 1))[used],
                                  member_pixels=size[used], pixels_assigned=assigned[used],
                                  min_valid_pixels=cnt.min(0)[used], class_code=zone_class[used],
                                  class_label=[CLASSES[int(v)] for v in zone_class[used]],
                                  reason=zone_reason[used], evidence=result.evidence[used],
                                  quality=result.quality[used]))
        for name, arr in result.metrics.items():
            frame[name] = np.asarray(arr)[used]
        tables.append(frame)
    absent = ~np.isfinite(vals).any(0)
    out["classes"][absent], out["reasons"][absent], out["rri_status"][absent] = 0, 1, 0

    def paint(a, fill, dtype):
        img = np.full((Hh, Ww), fill, dtype)
        img[rows, cols] = a
        return img
    rri_img = paint(out["rri"], np.nan, "float32")
    return ZonalProduct(paint(out["classes"], 0, "uint8"), paint(out["reasons"], 255, "uint8"),
                        paint(out["evidence"], 255, "uint8"), paint(out["quality"], 65535, "uint16"),
                        paint(out["zone"], 0, "int32"), paint(out["size"], 0, "uint16"), rri_img,
                        paint(out["rri_class"], 0, "uint8"), paint(out["rri_status"], 0, "uint8"),
                        pd.concat(tables, ignore_index=True), pd.DataFrame(level_rows))


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------
def qml(kind):
    rows = "\n".join(
        f'        <paletteEntry value="{v}" color="{col}" alpha="{a}" label="{v} {LABELS[kind].get(v, "")}"/>'
        for v, (col, a) in COLORS[kind].items())
    return f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="paletted" band="1" opacity="1" alphaBand="-1" nodataColor="">
      <rasterTransparency/>
      <colorPalette>
{rows}
      </colorPalette>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation saturation="0" grayscaleMode="0" invertColors="0" colorizeOn="0" colorizeRed="255" colorizeGreen="128" colorizeBlue="128" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
"""


RRI_QML = """<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="singlebandpseudocolor" band="1" opacity="1" alphaBand="-1" classificationMin="-0.5" classificationMax="1.5" nodataColor="">
      <rasterTransparency/>
      <rastershader>
        <colorrampshader colorRampType="INTERPOLATED" classificationMode="1" clip="0" minimumValue="-0.5" maximumValue="1.5">
          <item value="-0.5" color="#a50026" alpha="255" label="-0.5 further degradation"/>
          <item value="0" color="#fdae61" alpha="255" label="0 no recovery"/>
          <item value="0.5" color="#ffffbf" alpha="255" label="0.5 half recovered"/>
          <item value="1" color="#66bd63" alpha="255" label="1 back to baseline"/>
          <item value="1.5" color="#006837" alpha="255" label="1.5 above baseline"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
"""

PRODUCTS = {  # name: (attribute, dtype, nodata, style kind)
    "Zonal_Recovery_Classes": ("classes", "uint8", 0, "class"),
    "Zonal_RRI_Interpretation_Classes": ("rri_class", "uint8", 0, "rri"),
    "Zonal_RRI_Validity_Status": ("rri_status", "uint8", 255, "status"),
    "Zonal_RRI": ("rri", "float32", -9999.0, "rri_value"),
    "Zonal_Disturbance_Evidence": ("evidence", "uint8", 255, "evidence"),
    "Zonal_Indeterminate_Reasons": ("reasons", "uint8", 255, "reason"),
    "Zonal_Quality_Flags": ("quality", "uint16", 65535, None),
    "Zonal_Zone_Size_m": ("zone_size_m", "uint16", 0, None),
    "Zonal_Zone_ID": ("zone_id", "int32", 0, None),
}


def write_raster(path, array, dtype, nodata, profile, kind, member, tags):
    p = dict(driver="GTiff", height=array.shape[0], width=array.shape[1], count=1, dtype=dtype,
             crs=profile["crs"], transform=profile["transform"], nodata=nodata, compress="deflate",
             tiled=True, blockxsize=256, blockysize=256)
    a = array.astype(dtype, copy=True)
    if np.issubdtype(a.dtype, np.floating):
        a[~np.isfinite(a)] = nodata
    a[~member] = nodata
    if kind in COLORS:
        p["photometric"] = "palette"
    with rasterio.open(path, "w", **p) as dst:
        if kind in COLORS:  # colour table must precede pixel data in tiled GeoTIFFs
            cmap = {v: tuple(int(col[i:i + 2], 16) for i in (1, 3, 5)) + (alpha,)
                    for v, (col, alpha) in COLORS[kind].items() if v < 256}
            dst.write_colormap(1, cmap)
        dst.write(a, 1)
        dst.update_tags(**tags)
    if kind in COLORS:
        path.with_suffix(".qml").write_text(qml(kind), encoding="utf-8")
    elif kind == "rri_value":
        path.with_suffix(".qml").write_text(RRI_QML, encoding="utf-8")


def class_table(classes, member, area_m2, labels):
    v = classes[member]
    total = int(v.size)
    rows = []
    for code, label in labels.items():
        n = int((v == code).sum())
        rows.append(dict(code=code, label=label, pixels=n, percent=100 * n / total if total else np.nan,
                         hectares=n * area_m2 / 1e4))
    return pd.DataFrame(rows)


def target_row(zone, original, zonal, member):
    o = np.bincount(original[member], minlength=6)[:6]
    z = np.bincount(zonal[member], minlength=6)[:6]
    i0, i1, c0, c1 = int(o[5]), int(z[5]), int(o[1:4].sum()), int(z[1:4].sum())
    d0, d1 = int(o[1:5].sum()), int(z[1:5].sum())
    red = 100 * (i0 - i1) / i0 if i0 else np.nan
    inc = 100 * (c1 - c0) / c0 if c0 else np.nan
    inc_all = 100 * (d1 - d0) / d0 if d0 else np.nan
    return dict(zone=zone, population=int(member.sum()), I0_original_indeterminate=i0, I1_zonal_indeterminate=i1,
                indeterminate_reduction_percent=red, target_indeterminate_reduction_90_met=bool(red >= 90) if i0 else None,
                C0_original_classes_1to3=c0, C1_zonal_classes_1to3=c1, recovery_class_increase_percent=inc,
                target_recovery_increase_90_met=bool(inc >= 90) if c0 else None,
                D0_original_classes_1to4=d0, D1_zonal_classes_1to4=d1, all_class_increase_percent=inc_all,
                **{f"original_{k}": int(o[k]) for k in range(6)}, **{f"zonal_{k}": int(z[k]) for k in range(6)})


def class_map(classes, member, path, title, kind="class"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    rgba = np.zeros((6, 4))
    for v, (col, a) in COLORS[kind].items():
        rgba[v] = [int(col[i:i + 2], 16) / 255 for i in (1, 3, 5)] + [a / 255]
    img = np.where(member, classes, 0)
    fig, ax = plt.subplots(figsize=(9, 8), layout="constrained")
    ax.imshow(img, cmap=ListedColormap(rgba), vmin=-0.5, vmax=5.5, interpolation="nearest")
    ax.set_axis_off()
    ax.set_title(title)
    ax.legend(handles=[Patch(color=rgba[k], label=LABELS[kind][k]) for k in range(1, 6)], loc="lower left",
              fontsize=8, frameon=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def comparison_chart(original, zonal, member, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    a = np.bincount(original[member], minlength=6)[:6]
    b = np.bincount(zonal[member], minlength=6)[:6]
    fig, ax = plt.subplots(figsize=(10, 4.5), layout="constrained")
    ix = np.arange(1, 6)
    ax.bar(ix - .2, a[1:], .4, label="Original pixel/object ITSA (0.40 gate)", color=".65")
    ax.bar(ix + .2, b[1:], .4, label="Multi-scale zonal ITSA", color="#245a8d")
    ax.set_xticks(ix, ["Recovered", "Partial", "Not recovered", "No detectable\nloss", "Indeterminate"])
    ax.set(title=title, ylabel="Study pixels")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def polygons(zone_id, member, profile, zones, names, path):
    """Zone polygons (WGS84 GeoJSON) carrying zone attributes."""
    attributes = zones.set_index("zone_id")
    features = []
    for geom, value in shapes(zone_id, mask=member & (zone_id > 0), transform=profile["transform"],
                              connectivity=4):
        zid = int(value)
        if zid not in attributes.index:
            continue
        r = attributes.loc[zid]
        props = dict(zone_id=zid, zone_side_m=float(r.zone_side_m), barangay=names.get(int(r.barangay_index), ""),
                     class_code=int(r.class_code), class_label=str(r.class_label), evidence=int(r.evidence),
                     evidence_label=EVIDENCE.get(int(r.evidence), ""), quality=int(r.quality),
                     baseline=round(float(r.baseline), 4) if np.isfinite(r.baseline) else None,
                     event_value=round(float(r.event_value), 4) if np.isfinite(r.event_value) else None,
                     late_post=round(float(r.late_post), 4) if np.isfinite(r.late_post) else None,
                     rri=round(float(r.recovery_index), 4) if np.isfinite(r.recovery_index) else None,
                     rri_class=int(r.rri_class))
        features.append(dict(type="Feature", properties=props,
                             geometry=transform_geom(profile["crs"], "EPSG:4326", geom)))
    path.write_text(json.dumps(dict(type="FeatureCollection", features=features)), encoding="utf-8")
    return len(features)


# ---------------------------------------------------------------------------
# Driver: outputs go into the standard run layout
#   A. RRI/{Raster Files,Graphs,CSV}/{Overall,Per Barangay/<name_code>}
#   B. ITSA/{Overall,Per Barangay/<name_code>}/{Raster Files,ITSA Graph,CSV}
# Every zonal file name starts with "Zonal_" so it never replaces a pixel/object
# product of the same run.
# ---------------------------------------------------------------------------
ITSA_PRODUCTS = ("Zonal_Recovery_Classes", "Zonal_Disturbance_Evidence", "Zonal_Indeterminate_Reasons",
                 "Zonal_Quality_Flags", "Zonal_Zone_Size_m", "Zonal_Zone_ID")
RRI_PRODUCTS = ("Zonal_RRI", "Zonal_RRI_Interpretation_Classes", "Zonal_RRI_Validity_Status")


def existing_zonal_files(run_root):
    return [p for method in ("A. RRI", "B. ITSA") for p in (Path(run_root) / method).rglob("Zonal_*")]


def run_zonal(run_root, c, replace=False):
    """Classify zones for a completed run and write Zonal_* files into A. RRI / B. ITSA."""
    from .study_layout import folders
    run_root = Path(run_root)
    rasters = run_root / "B. ITSA" / "Overall" / "Raster Files"
    previous = existing_zonal_files(run_root)
    if previous and not replace:
        raise FileExistsError(f"{len(previous)} Zonal_* files already exist in {run_root}; "
                              "pass replace=True (run_zonal.py --replace) to regenerate them.")
    for p in previous:
        p.unlink()
    with rasterio.open(rasters / "Aligned_Coherence_Observations.tif") as src:
        stack = src.read().astype("float32")
        if src.nodata is not None:
            stack[stack == src.nodata] = np.nan
        profile = dict(crs=src.crs, transform=src.transform)
        pixel = abs(src.transform.a)
        area = abs(src.transform.a * src.transform.e)
    with rasterio.open(rasters / "Study_Mask.tif") as src:
        inside = src.read(1) == 1
    with rasterio.open(rasters / "Barangay_Owner.tif") as src:
        owner = src.read(1).astype("int64")
    with rasterio.open(rasters / "Original_ITSA_Recovery_Classes.tif") as src:
        original = src.read(1).astype("uint8")
    stack[:, ~inside] = np.nan
    obs = observations(c)
    c = effective_settings(c, obs)
    if len(obs) != stack.shape[0]:
        raise ValueError("Observation count does not match the aligned stack bands.")
    from .analysis_boundaries import load_barangays
    barangays, _ = load_barangays(c.BARANGAY_BOUNDARY, profile["crs"], c)
    names = {0: "", **{i: b.name for i, b in enumerate(barangays, 1)}}

    product = zonal_classification(stack, inside, owner, obs, c, pixel)
    tags = dict(method="multi-scale zonal ITSA (observation-anchored); see sarflow/zonal_recovery.py",
                zone_base_m=str(round(setting(c, "ZONAL_BASE_CELLS", 5) * pixel, 1)),
                zone_levels=str(setting(c, "ZONAL_LEVELS", 4)),
                caveat="Exploratory zonal product; event-pair decorrelation may reflect weather; "
                       "coherence estimator bias uncorrected; no field validation.",
                classes=json.dumps(CLASSES), quality_bits=json.dumps(QUALITY))

    def export_zone(folder, name, window, member):
        itsa, rri = folders(run_root, "ITSA", folder), folders(run_root, "RRI", folder)
        sl = np.s_[int(window.row_off):int(window.row_off + window.height),
                   int(window.col_off):int(window.col_off + window.width)]
        prof = dict(crs=profile["crs"], transform=window_transform(window, profile["transform"]))
        for pname, (attr, dtype, nodata, kind) in PRODUCTS.items():
            target = itsa.rasters if pname in ITSA_PRODUCTS else rri.rasters
            write_raster(target / f"{pname}.tif", getattr(product, attr)[sl], dtype, nodata, prof, kind,
                         member, dict(tags, zone=name))
        cls = product.classes[sl]
        # B. ITSA CSV
        class_table(cls, member, area, CLASSES).to_csv(itsa.csv / "Zonal_Class_Area_Statistics.csv", index=False)
        class_table(product.evidence[sl], member, area, EVIDENCE).to_csv(
            itsa.csv / "Zonal_Disturbance_Evidence_Statistics.csv", index=False)
        ind = member & (cls == 5)
        class_table(product.reasons[sl], ind, area, REASONS).to_csv(
            itsa.csv / "Zonal_Indeterminate_Reasons.csv", index=False)
        sizes = product.zone_size_m[sl][member]
        pd.DataFrame([dict(zone_side_m=int(s), pixels=int((sizes == s).sum()),
                           percent=100 * (sizes == s).mean()) for s in np.unique(sizes)]).to_csv(
            itsa.csv / "Zonal_Resolution_Levels.csv", index=False)
        row = target_row(name, original[sl], cls, member)
        pd.DataFrame([row]).to_csv(itsa.csv / "Zonal_Target_Assessment.csv", index=False)
        tr = np.bincount(original[sl][member].astype(int) * 6 + cls[member], minlength=36).reshape(6, 6)
        pd.DataFrame(tr, index=[f"original_{CLASSES[i]}" for i in range(6)],
                     columns=[f"zonal_{CLASSES[j]}" for j in range(6)]).to_csv(
            itsa.csv / "Zonal_Original_to_Zonal_Transitions.csv")
        # A. RRI CSV
        class_table(product.rri_class[sl], member, area, LABELS["rri"]).to_csv(
            rri.csv / "Zonal_RRI_Class_Statistics.csv", index=False)
        class_table(product.rri_status[sl], member, area, STATUS).to_csv(
            rri.csv / "Zonal_RRI_Validity_Statistics.csv", index=False)
        # Graphs
        class_map(cls, member, itsa.graphs / "Zonal_Recovery_Classes_Map.png",
                  f"{name}: multi-scale zonal recovery classes")
        comparison_chart(original[sl], cls, member, itsa.graphs / "Zonal_Original_vs_Zonal_Counts.png",
                         f"{name}: original vs zonal classification")
        class_map(product.rri_class[sl], member, rri.graphs / "Zonal_RRI_Interpretation_Map.png",
                  f"{name}: zonal RRI interpretation", kind="rri")
        return row

    full = Window(0, 0, inside.shape[1], inside.shape[0])
    rows = [export_zone(None, "Overall Infanta", full, inside)]
    overall = folders(run_root, "ITSA")
    product.levels.to_csv(overall.csv / "Zonal_Zone_Size_Levels.csv", index=False)
    zones = product.zones.copy()
    zones.insert(4, "barangay", zones.barangay_index.map(names))
    zones.to_csv(overall.csv / "Zonal_Zone_Attributes.csv", index=False)
    for i, record in enumerate(barangays, 1):
        rr, cc = np.nonzero(inside & (owner == i))
        if not rr.size:
            continue
        win = Window(int(cc.min()), int(rr.min()), int(np.ptp(cc) + 1), int(np.ptp(rr) + 1))
        member = (owner[int(win.row_off):int(win.row_off + win.height), int(win.col_off):int(win.col_off + win.width)] == i) \
            & inside[int(win.row_off):int(win.row_off + win.height), int(win.col_off):int(win.col_off + win.width)]
        rows.append(export_zone(record.folder, record.name, win, member))
    pd.DataFrame(rows).to_csv(overall.csv / "Zonal_Target_Assessment_All_Barangays.csv", index=False)
    n = polygons(product.zone_id, inside, profile, product.zones, names,
                 overall.rasters / "Zonal_Recovery_Zones.geojson")

    sens = []
    for base in setting(c, "ZONAL_SENSITIVITY_BASE_CELLS", (3, 7)):
        alt = zonal_classification(stack, inside, owner, obs, c, pixel, base_cells=base)
        sens.append(dict(scenario=f"base_{round(base * pixel)}m", **target_row("Overall Infanta", original, alt.classes, inside)))
    for name, change in (("tolerance_0.03", dict(RECOVERY_TOLERANCE=.03)), ("tolerance_0.07", dict(RECOVERY_TOLERANCE=.07)),
                         ("loss_min_0.05", dict(EVENT_LOSS_MINIMUM=.05)),
                         ("no_barangay_split", dict(ZONAL_RESPECT_BARANGAYS=False))):
        from types import SimpleNamespace
        alt_c = SimpleNamespace(**{**{k: getattr(c, k) for k in dir(c) if not k.startswith("_")}, **change})
        alt = zonal_classification(stack, inside, owner, obs, alt_c, pixel)
        sens.append(dict(scenario=name, **target_row("Overall Infanta", original, alt.classes, inside)))
    pd.DataFrame(sens).to_csv(overall.csv / "Zonal_Sensitivity.csv", index=False)
    (overall.csv / "Zonal_Method_README.txt").write_text(
        (__doc__ or "") + f"\nZone polygons: B. ITSA/Overall/Raster Files/Zonal_Recovery_Zones.geojson "
        f"({n} features, WGS84).\nQuality bits: " + json.dumps(QUALITY, indent=1) + "\n", encoding="utf-8")
    return run_root, rows[0]
