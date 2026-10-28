"""Pure numerical recovery kernel. Status is separate from duration."""
from dataclasses import dataclass
from types import ModuleType
import warnings
import numpy as np
from .temporal import Observation, day_offsets

STATUS = {1: "invalid_or_unstable_baseline", 2: "no_detectable_event_loss",
          3: "observed_sustained_return", 4: "no_confirmed_return_by_last_observation",
          5: "insufficient_event_or_post_observations"}


@dataclass
class Metrics:
    baseline: np.ndarray
    spread: np.ndarray
    threshold: np.ndarray
    loss: np.ndarray
    bands: np.ndarray


@dataclass
class FractionMetrics:
    baseline: np.ndarray
    spread: np.ndarray
    threshold: np.ndarray
    loss: np.ndarray
    fraction: np.ndarray
    classes: np.ndarray


def calculate_fraction(stack: np.ndarray, obs: list[Observation], c: ModuleType,
                       tolerance: float | None = None,
                       confirmations: int | None = None) -> FractionMetrics:
    """End-state coherence recovery fraction and 4-class interpretation.

    Recovery fraction is intentionally not clipped. Classes 1-3 require valid
    values in all three latest post pairs; class 1 additionally requires both
    latest consecutive pairs to pass the baseline-minus-tolerance rule.
    Class 4 needs only a valid baseline and event pair. Class 0 is NoData.
    """
    tol = c.RECOVERY_TOLERANCE if tolerance is None else tolerance
    pre, event, post = stack[:8], stack[8], stack[9:]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        baseline = np.nanmedian(pre, axis=0)
        spread = 1.4826 * np.nanmedian(np.abs(pre - baseline), axis=0)
    baseline_ok = ((np.isfinite(pre).sum(axis=0) >= c.MIN_PRE)
                   & (baseline >= c.BASELINE_MIN) & (spread <= c.BASELINE_MAD_MAX))
    loss = baseline - event
    loss_limit = np.maximum(c.LOSS_ABSOLUTE,
                            c.LOSS_MAD_MULTIPLIER * np.maximum(spread, c.NOISE_FLOOR))
    detectable = np.isfinite(event) & (loss >= loss_limit)
    threshold = np.maximum(0, baseline - tol)
    classes = np.full(baseline.shape, -9999, dtype="int16")
    fraction = np.full(baseline.shape, np.nan, dtype="float32")
    valid_reference = baseline_ok & np.isfinite(event)
    classes[valid_reference & ~detectable] = 4

    # The requested end-state statistic is the median of exactly the latest
    # three post-event pairs. Do not silently substitute older dates if missing.
    tail_n = c.ENDSTATE_POST_PAIRS
    tail = post[-tail_n:]
    tail_complete = np.isfinite(tail).all(axis=0)
    eligible = valid_reference & detectable & tail_complete
    if tail.shape[0] != tail_n:
        eligible[:] = False
    if eligible.any():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            post_end = np.nanmedian(tail, axis=0)
        fraction[eligible] = ((post_end[eligible] - event[eligible]) / loss[eligible]).astype("float32")
        k = c.RECOVERY_CONFIRMATIONS if confirmations is None else confirmations
        latest_values = post[-k:]
        latest = obs[-k:]
        confirmation_ok = np.isfinite(latest_values).all(axis=0) & (latest_values >= threshold).all(axis=0)
        temporal_ok = True
        for left, right in zip(latest, latest[1:]):
            gap = (right.reference - left.reference).total_seconds() / 86400
            temporal_ok = temporal_ok and gap <= c.MAX_CONFIRM_GAP_DAYS
            temporal_ok = temporal_ok and (c.ALLOW_OVERLAPPING_CONFIRMATION_PAIRS or left.end <= right.start)
        recovered = eligible & confirmation_ok & temporal_ok
        partial = eligible & ~recovered & (fraction > 0)
        not_recovered = eligible & ~recovered & (fraction <= 0)
        classes[recovered] = 1
        classes[partial] = 2
        classes[not_recovered] = 3
    return FractionMetrics(baseline, spread, threshold, loss, fraction, classes)


def calculate(stack: np.ndarray, obs: list[Observation], c: ModuleType,
              tolerance: float | None = None, confirmations: int | None = None) -> Metrics:
    """stack=(17,...); NaN means unobserved. Never interpolate missing time steps."""
    tol = c.RECOVERY_TOLERANCE if tolerance is None else tolerance
    k = c.CONFIRMATIONS if confirmations is None else confirmations
    pre, event, post = stack[:8], stack[8], stack[9:]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        baseline = np.nanmedian(pre, axis=0)
        spread = 1.4826 * np.nanmedian(np.abs(pre - baseline), axis=0)
    eligible_baseline = ((np.isfinite(pre).sum(axis=0) >= c.MIN_PRE)
                         & (baseline >= c.BASELINE_MIN) & (spread <= c.BASELINE_MAD_MAX))
    loss = baseline - event
    loss_cut = np.maximum(c.LOSS_ABSOLUTE,
                          c.LOSS_MAD_MULTIPLIER * np.maximum(spread, c.NOISE_FLOOR))
    detectable = np.isfinite(event) & (loss >= loss_cut)
    threshold = np.maximum(0, baseline - tol)
    status = np.ones(baseline.shape, dtype="float32")
    status[eligible_baseline] = 5
    status[eligible_baseline & np.isfinite(event) & ~detectable] = 2
    eligible = eligible_baseline & detectable & (np.isfinite(post).sum(axis=0) >= c.MIN_POST)
    status[eligible] = 4
    days = day_offsets(obs, c.EVENT_ORIGIN)[9:]
    passing = np.isfinite(post) & (post >= threshold)
    # Terminal observed run: any later dip OR missing observation breaks persistence.
    # A large temporal gap or overlapping acquisition intervals also breaks the run.
    run = np.zeros(baseline.shape, dtype="int16")
    first = np.full(baseline.shape, 8, dtype="int16")
    active = np.ones(baseline.shape, dtype=bool)
    for j in range(7, -1, -1):
        if j < 7:
            gap_ok = days[j + 1] - days[j] <= c.MAX_CONFIRM_GAP_DAYS
            overlap_ok = (c.ALLOW_OVERLAPPING_CONFIRMATION_PAIRS
                          or obs[j + 9].end <= obs[j + 10].start)
            active &= gap_ok and overlap_ok
        active &= passing[j]
        run[active] += 1
        first[active] = j
    safe_first = np.minimum(first, 7)
    observed = (eligible & (run >= k)
                & ((days[-1] - days[safe_first]) >= c.MIN_CONFIRM_SPAN_DAYS))
    status[observed] = 3
    result = np.full((5, *baseline.shape), np.nan, dtype="float32")
    result[0][observed] = days[safe_first][observed]
    # Operational observation bracket, NOT a confidence interval for physical recovery.
    # Lower bound is preceding pair reference, or event origin if first post qualifies.
    lower_days = np.r_[0.0, days[:-1]]
    result[1][observed] = lower_days[safe_first][observed]
    confirm_index = np.minimum(safe_first + k - 1, 7)
    # Minimum span may require more than k observations to confirm.
    for j in range(8):
        need_later = (days[confirm_index] - days[safe_first]) < c.MIN_CONFIRM_SPAN_DAYS
        confirm_index = np.where(need_later, np.minimum(confirm_index + 1, 7), confirm_index)
    result[2][observed] = days[confirm_index][observed]
    result[3] = status
    valid_post = np.isfinite(post)
    last = np.max(np.where(valid_post, np.arange(8).reshape((8,) + (1,) * baseline.ndim), -1), axis=0)
    followup = eligible_baseline & detectable & (last >= 0)
    result[4][followup] = days[np.maximum(last, 0)][followup]
    return Metrics(baseline, spread, threshold, loss, result)
