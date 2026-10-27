"""Visualization-only adapters; all analytical eligibility comes from existing kernels."""
from dataclasses import dataclass
import numpy as np


@dataclass
class PlotData:
    sample: object
    zones: dict
    obs: list
    c: object
    fraction: np.ndarray
    table: object
    metrics: object
    anomaly: np.ndarray
    groups: list
    post_indices: list
    normalized: np.ndarray
    threshold: np.ndarray
    eligible: np.ndarray


def prepare(sample, zones, obs, c, fraction, metrics_module):
    phases = [o.phase for o in obs]
    pre = [i for i, phase in enumerate(phases) if phase == 'pre']
    event = phases.index('event')
    post = [i for i, phase in enumerate(phases) if phase == 'post']
    if (len(pre), phases.count('event'), len(post)) not in ((1, 1, 1), (8, 1, 8)):
        raise ValueError('Only the existing 1/1/1 and 8/1/8 workflows are supported')
    table = sample.table([o.name for o in obs])
    if len(pre) == 1:
        metrics = metrics_module.calculate_fraction(sample.values, obs, c)
        eligible = np.isin(metrics.classes, [1, 2, 3])
    else:
        metrics = metrics_module.calculate(sample.values, obs, c)
        loss_cut = np.maximum(c.LOSS_ABSOLUTE,
                              c.LOSS_MAD_MULTIPLIER*np.maximum(metrics.spread, c.NOISE_FLOOR))
        eligible = (metrics.loss >= loss_cut) & np.isin(metrics.bands[3], [3, 4])
    normalized = ((sample.values[post][:, eligible]-sample.values[event, eligible])
                  / metrics.loss[eligible])
    threshold = 1-c.RECOVERY_TOLERANCE/metrics.loss[eligible]
    groups = [np.median(sample.values[pre], axis=0), sample.values[event],
              np.median(sample.values[post], axis=0)]
    return PlotData(sample, zones, obs, c, fraction, table, metrics,
                    sample.values-metrics.baseline, groups, post, normalized, threshold, eligible)
