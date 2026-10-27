"""Descriptive coherence measurements; classification requires an approved rule."""
from dataclasses import dataclass
import warnings
import numpy as np

BAND_NAMES = (
    "pre_event_median_coherence", "event_straddling_coherence",
    "latest_post_coherence", "event_minus_pre_median",
    "latest_post_minus_pre_median", "latest_post_minus_event",
    "valid_pre_observation_count",
)
from .recovery import LABELS as CLASS_LABELS


@dataclass
class Measurements:
    baseline: np.ndarray
    event: np.ndarray
    latest_post: np.ndarray
    pre_count: np.ndarray
    bands: np.ndarray


def measure(stack, obs) -> Measurements:
    pre_indices = [i for i, o in enumerate(obs) if o.phase == "pre"]
    event_index = next(i for i, o in enumerate(obs) if o.phase == "event")
    post_indices = [i for i, o in enumerate(obs) if o.phase == "post"]
    pre = stack[pre_indices]
    pre_count = np.isfinite(pre).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        baseline = np.nanmedian(pre, axis=0)
    event = stack[event_index]
    latest = stack[post_indices[-1]]  # latest actual pair; no substitution/interpolation
    count_band = np.where(pre_count > 0, pre_count, np.nan)
    bands = np.stack((baseline, event, latest, event-baseline,
                      latest-baseline, latest-event, count_band)).astype("float32")
    return Measurements(baseline, event, latest, pre_count, bands)

