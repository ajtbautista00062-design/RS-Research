"""Validated explicit observations and source-validity-constrained aligned reading."""
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
import re
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.features import geometry_window
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window, transform
from rasterio.windows import bounds, from_bounds
from rasterio.warp import transform_bounds
from .temporal import Observation, configured_pair, local_path, timestamp
from .spatial import read_zones, mask_for, inventory
from .analysis_boundaries import load_barangays, owners


@dataclass
class StudyInputs:
    observations: list[Observation]
    stack: np.ndarray
    inside: np.ndarray
    ownership: np.ndarray
    barangays: list
    profile: dict
    source_window: Window
    inventory: list[dict]
    boundary_audit: dict


def input_path(value, root):
    """Input rasters are read-only: absolute paths anywhere, relative paths from root."""
    path = Path(value).expanduser()
    return (path if path.is_absolute() else Path(root) / path).resolve()


SLOTS = (*(f'Pre{i}' for i in range(1, 9)), 'During', *(f'Post{i}' for i in range(1, 9)))


def observations(c) -> list[Observation]:
    """Observations from config.py: up to 8 PRE, exactly 1 event-straddling, up to 8 POST.

    Empty Pre/Post slots are skipped (at least ``MIN_PRE_INPUTS``/``MIN_POST_INPUTS``,
    default 3, are required). Dates come from ``PAIR_DATES`` or, when blank, from one
    unambiguous YYYYMMDD_YYYYMMDD pair in the filename. Raster paths may be
    absolute (any drive) or relative to PROJECT_ROOT.
    """
    if set(c.PATHS) != set(SLOTS):
        raise ValueError('config.PATHS must define exactly Pre1..Pre8, During, Post1..Post8.')
    event_start, event_end = timestamp(c.EVENT_START), timestamp(c.EVENT_END)
    errors, records = [], []
    if not event_start <= timestamp(c.EVENT_ORIGIN) <= event_end:
        errors.append('EVENT_ORIGIN must lie inside EVENT_START..EVENT_END.')
    if c.REFERENCE_DATE not in ('end', 'midpoint'):
        errors.append('REFERENCE_DATE must be end or midpoint.')
    for slot in SLOTS:
        if not c.PATHS[slot]:
            if slot == 'During':
                errors.append('During/Event_Straddler: the event-straddling raster is required.')
            continue
        try:
            path = input_path(c.PATHS[slot], c.PROJECT_ROOT)
            if not path.is_file():
                raise ValueError(f'raster does not exist: {path}')
            master, slave = configured_pair(c, slot, path)
            start, end = sorted((master, slave))
            if start == end:
                raise ValueError('pair acquisitions must differ')
            phase = 'pre' if slot.startswith('Pre') else 'post' if slot.startswith('Post') else 'event'
            if phase == 'pre' and end >= event_start:
                raise ValueError('PRE pair must finish before the event window')
            if phase == 'post' and start <= event_end:
                raise ValueError('POST pair must start after the event window')
            if phase == 'event' and (end < event_start or start > event_end):
                raise ValueError('During pair must overlap the event window')
            reference = end if c.REFERENCE_DATE == 'end' else start + (end - start) / 2
            records.append(Observation(f'{slot}_{start:%Y%m%d}_{end:%Y%m%d}', path, phase,
                                       master, slave, start, end, reference))
        except (ValueError, TypeError, KeyError) as exc:
            errors.append(f'{slot}: {exc}')
    records.sort(key=lambda o: o.reference)
    counts = {p: sum(o.phase == p for o in records) for p in ('pre', 'event', 'post')}
    if not errors:
        if counts['pre'] < getattr(c, 'MIN_PRE_INPUTS', 3) or counts['post'] < getattr(c, 'MIN_POST_INPUTS', 3):
            errors.append(f'Need at least {getattr(c, "MIN_PRE_INPUTS", 3)} PRE and '
                          f'{getattr(c, "MIN_POST_INPUTS", 3)} POST rasters; got {counts}.')
        if len({o.path for o in records}) != len(records):
            errors.append('Duplicate input paths are not allowed.')
        if len({(o.start, o.end) for o in records}) != len(records):
            errors.append('Duplicate acquisition pairs are not allowed.')
        if not all(a.reference < b.reference and a.end < b.end for a, b in zip(records, records[1:])):
            errors.append('Observations must have distinct increasing reference and pair-end dates.')
        spans = [(o.end - o.start).total_seconds() / 86400 for o in records]
        if records and max(spans) / min(spans) > c.MAX_PAIR_BASELINE_RATIO and not c.ALLOW_BASELINE_MISMATCH:
            errors.append('Pair temporal baselines differ too much: harmonize pairs or set '
                          'ALLOW_BASELINE_MISMATCH with a rationale in COMPARABILITY_NOTES.')
    if getattr(c, 'REQUIRE_COMPARABILITY_CONFIRMATION', True) and (
            not c.ACQUISITION_COMPARABILITY_CONFIRMED or not c.COMPARABILITY_NOTES.strip()):
        errors.append('Set ACQUISITION_COMPARABILITY_CONFIRMED = True and describe track, orbit, '
                      'polarization, bursts and processing in COMPARABILITY_NOTES.')
    if errors:
        raise ValueError('Configuration errors:\n - ' + '\n - '.join(errors))
    return records


def effective_settings(c, obs):
    """Settings copy whose ITSA minimum valid-date counts fit the supplied inputs."""
    from types import SimpleNamespace
    counts = {ph: sum(o.phase == ph for o in obs) for ph in ('pre', 'post')}
    values = {k: getattr(c, k) for k in dir(c) if not k.startswith('_')}
    values.update(EXPECTED_PRE_COUNT=counts['pre'], EXPECTED_POST_COUNT=counts['post'],
                  ITSA_MIN_PRE=min(c.ITSA_MIN_PRE, counts['pre']),
                  ITSA_MIN_POST=min(c.ITSA_MIN_POST, counts['post']))
    return SimpleNamespace(**values)


def load(c) -> StudyInputs:
    obs = observations(c)
    with ExitStack() as context:
        sources = [context.enter_context(rasterio.open(o.path)) for o in obs]
        ref = sources[0]
        metadata = inventory(sources, obs, c)
        zones = read_zones(c.Study_Area, ref.crs, c.PROJECT_ROOT)
        barangays, audit = load_barangays(c.BARANGAY_BOUNDARY, ref.crs, c)
        if not barangays or audit['skipped'] or audit['overlaps']:
            raise ValueError(f'Boundary must have verified named nonoverlapping polygons: {audit}')
        crop = geometry_window(ref, [g for _,g in zones]).intersection(Window(0,0,ref.width,ref.height))
        aff = transform(crop, ref.transform)
        shape = (int(crop.height), int(crop.width))
        inside = mask_for([g for _,g in zones], shape, aff, False)
        owner = owners(barangays, shape, aff)
        owner[~inside] = 0
        if np.any(inside & (owner == 0)):
            raise ValueError('Study cells lack named barangay ownership.')
        layers = []
        for source_index,src in enumerate(sources):
            source_bounds=transform_bounds(ref.crs,src.crs,*bounds(crop,ref.transform),densify_pts=21)
            raw_window=from_bounds(*source_bounds,src.transform)
            left=int(np.floor(raw_window.col_off))-2
            top=int(np.floor(raw_window.row_off))-2
            right=int(np.ceil(raw_window.col_off+raw_window.width))+2
            bottom=int(np.ceil(raw_window.row_off+raw_window.height))+2
            raw_window=Window(left,top,right-left,bottom-top).intersection(Window(0,0,src.width,src.height))
            original=src.read(c.BAND,window=raw_window,masked=True).astype('float32').filled(np.nan)
            original=(original*src.scales[c.BAND-1]+src.offsets[c.BAND-1])*c.COHERENCE_SCALE+c.COHERENCE_OFFSET
            values=original[np.isfinite(original)]
            if values.size and (values.min() < -1e-6 or values.max() > 1+1e-6):
                raise ValueError(f'Invalid source magnitude before interpolation: {src.name}')
            metadata[source_index]['raw_validation_scope']='study bounding source window plus 2-pixel bilinear halo'
            metadata[source_index]['raw_validation_window']=str(raw_window)
            metadata[source_index]['raw_validated_pixels']=len(values)
            kwargs = dict(crs=ref.crs, transform=ref.transform, width=ref.width,
                          height=ref.height, dtype='float32', nodata=np.nan)
            continuous = context.enter_context(WarpedVRT(src, **kwargs, resampling=Resampling[c.RESAMPLING]))
            nearest = context.enter_context(WarpedVRT(src, **kwargs, resampling=Resampling.nearest))
            a = continuous.read(c.BAND, window=crop, masked=True).filled(np.nan)
            footprint = nearest.read(c.BAND, window=crop, masked=True).filled(np.nan)
            a[~np.isfinite(footprint) | ~inside] = np.nan
            a = (a*src.scales[c.BAND-1]+src.offsets[c.BAND-1])*c.COHERENCE_SCALE+c.COHERENCE_OFFSET
            finite = a[np.isfinite(a)]
            if finite.size and (finite.min() < -1e-6 or finite.max() > 1+1e-6):
                raise ValueError(f'Invalid coherence magnitude in {src.name}')
            layers.append(np.clip(a, 0, 1).astype('float32'))
        profile = ref.profile.copy()
        profile.update(height=shape[0], width=shape[1], transform=aff)
    return StudyInputs(obs, np.stack(layers), inside, owner, barangays, profile,
                       crop, metadata, audit)
