"""Exact plotting-input cache; no replacement sample or raster computation."""
import hashlib
import json
from pathlib import Path
import numpy as np

SETTINGS = ('MIN_PRE','MIN_POST','BASELINE_MIN','BASELINE_MAD_MAX','NOISE_FLOOR',
            'LOSS_ABSOLUTE','LOSS_MAD_MULTIPLIER','RECOVERY_TOLERANCE','CONFIRMATIONS',
            'RECOVERY_CONFIRMATIONS','ENDSTATE_POST_PAIRS','MIN_CONFIRM_SPAN_DAYS',
            'MAX_CONFIRM_GAP_DAYS','ALLOW_OVERLAPPING_CONFIRMATION_PAIRS','REFERENCE_DATE',
            'EVENT_ORIGIN','RESAMPLING','COHERENCE_SCALE','COHERENCE_OFFSET','BAND',
            'ZERO_IS_INVALID','Study_Area','ALL_TOUCHED','RANDOM_SEED','SAMPLE_PIXELS')


def signature(c, obs, metrics):
    inputs = {'settings': {key: getattr(c, key, None) for key in SETTINGS},
              'observations': [[o.name, str(o.path), o.phase, str(o.master), str(o.slave)] for o in obs],
              'kernel_sha256': hashlib.sha256(Path(metrics.__file__).read_bytes()).hexdigest()}
    return hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()


def write_cache(out, sample, zones, obs, c, fraction, metrics):
    names = list(zones)
    arrays = {'signature': np.array(signature(c, obs, metrics)),
              'zone_names': np.array(names, dtype='U'), 'fraction': fraction}
    for prefix, current in [('global', sample), *[(f'zone_{i}', zones[name]) for i,name in enumerate(names)]]:
        for key in ('values', 'sums', 'coordinates', 'keys'):
            arrays[f'{prefix}_{key}'] = getattr(current, key)
        arrays[f'{prefix}_count'] = np.array(current.count, dtype='int64')
    np.savez_compressed(out/'publication_plot_cache.npz', **arrays)


def read_cache(out, c, obs, metrics, sample_class):
    path = out/'publication_plot_cache.npz'
    if not path.is_file():
        raise FileNotFoundError('publication_plot_cache.npz missing: run this project’s analysis once')
    with np.load(path, allow_pickle=False) as archive:
        if archive['signature'].item() != signature(c, obs, metrics):
            raise ValueError('Plot cache no longer matches configured analytical settings; rerun analysis')
        def restore(prefix):
            sample = sample_class(c.SAMPLE_PIXELS, c.RANDOM_SEED)
            for key in ('values','sums','coordinates','keys'):
                setattr(sample, key, archive[f'{prefix}_{key}'].copy())
            sample.count = int(archive[f'{prefix}_count'])
            return sample
        sample = restore('global')
        zones = {str(name): restore(f'zone_{i}') for i,name in enumerate(archive['zone_names'])}
        fraction = archive['fraction'].copy()
    return sample, zones, fraction
