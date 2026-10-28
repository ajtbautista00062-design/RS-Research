"""ITSA-dependent operational recovery; uncertainty is exported separately."""
from dataclasses import dataclass
import numpy as np
from .itsa import fit_itsa, median
from .temporal import day_offsets

CLASSES = {0: "NoData", 1: "Coherence recovered", 2: "Partially recovered",
           3: "Not recovered / further degraded", 4: "No detectable event-associated loss",
           5: "Indeterminate"}
REASONS = {0: "classified", 1: "structurally unavailable", 2: "insufficient PRE/POST",
           3: "unstable/nonfinite model", 4: "low baseline", 5: "conflicting disturbance",
           6: "loss near threshold or unsupported", 7: "late observations unavailable",
           8: "model/observation contradiction", 9: "ambiguous recovery",
           10: "spatial evidence conflict", 11: "outside analysis mask or ancillary exclusion"}
QUALITY = {1: "weak level-change SE evidence", 2: "wide modeled progress uncertainty",
           4: "below minimum mapping area", 8: "partial temporal support",
           16: "unresolved source semantics", 32: "spatially disputed"}


@dataclass
class Assessment:
    classes: np.ndarray
    quality: np.ndarray
    reasons: np.ndarray
    metrics: dict
    fit: object
    gates: dict | None = None


def assess(values, obs, c, loss_floor=None, tolerance=None):
    shape = values.shape[1:]
    y = np.asarray(values).reshape(len(obs), -1)
    fit = fit_itsa(y, obs, c)
    b, baseline = fit.beta, fit.baseline
    days = day_offsets(obs, c.EVENT_ORIGIN)
    post = np.flatnonzero([o.phase == "post" for o in obs])
    event = np.flatnonzero([o.phase == "event" for o in obs])
    if len(event) != 1 or len(post) < 3:
        raise ValueError("Require one event pair and at least three POST pairs.")
    floor = c.EVENT_LOSS_MINIMUM if loss_floor is None else loss_floor
    tol = c.RECOVERY_TOLERANCE if tolerance is None else tolerance
    k = c.LATE_POST_CONFIRMATIONS
    count = y.shape[1]
    classes = np.full(count, 5, dtype="uint8")
    reason = np.full(count, 9, dtype="uint8")
    quality = np.zeros(count, dtype="uint16")
    threshold = np.full(count, floor)
    if c.USE_ADAPTIVE_DISTURBANCE_THRESHOLD:
        threshold = np.maximum(threshold, c.PRE_RESIDUAL_NOISE_MULTIPLIER *
                               np.maximum(fit.noise, c.ITSA_NOISE_FLOOR))
    p = y[post]
    first_index = np.argmax(np.isfinite(p), axis=0)
    last_index = len(post) - 1 - np.argmax(np.isfinite(p[::-1]), axis=0)
    first_day, last_day = days[post[first_index]], days[post[last_index]]
    slope = b[1] + b[3]
    initial = b[0] + b[2] + slope * first_day
    final = b[0] + b[2] + slope * last_day
    progress = final - initial
    loss = -b[2]
    event_loss = baseline - y[event[0]]
    early = median(p[:k])
    late = median(p[-k:])
    deficit = baseline - final
    early_deficit = baseline - initial
    margin = c.EFFECT_DECISION_MARGIN
    support = (fit.pre_count >= c.ITSA_MIN_PRE) & (fit.post_count >= c.ITSA_MIN_POST)
    finite = (np.isfinite(b).all(0) & np.isfinite(threshold) &
              (fit.rmse <= c.MAX_ITSA_RMSE))
    predictions = b[0] + b[2] + slope * days[post, None]
    finite &= ((predictions >= -0.05) & (predictions <= 1.05)).all(0)
    usable = support & finite & (baseline >= c.ITSA_BASELINE_MINIMUM)
    reason[~support] = 2
    reason[support & ~finite] = 3
    reason[support & finite & (baseline < c.ITSA_BASELINE_MINIMUM)] = 4
    absent = ~np.isfinite(y).any(0)
    classes[absent], reason[absent] = 0, 1
    quality[(loss < c.UNCERTAINTY_MULTIPLIER * fit.level_se) | ~np.isfinite(fit.level_se)] |= 1
    # This is an uncertainty flag, never an eligibility gate or a p-value test.
    quality[fit.progress_se * c.UNCERTAINTY_MULTIPLIER > np.abs(progress)] |= 2
    quality[(fit.pre_count < sum(o.phase == "pre" for o in obs)) |
            (fit.post_count < len(post))] |= 8
    if not c.SOURCE_VALIDITY_CONFIRMED:
        quality[:] |= 16
    support_loss = c.EVENT_SUPPORT_MINIMUM
    event_positive = event_loss >= support_loss
    early_positive = (early_deficit >= support_loss) & (baseline - early >= support_loss)
    opposite = event_loss <= -support_loss
    credible = usable & (loss >= threshold + margin) & (event_positive | early_positive) & ~opposite
    no_loss = (usable & (loss <= threshold - margin) & ~event_positive &
               (early_deficit < support_loss - margin) & (baseline - early < support_loss - margin))
    classes[no_loss], reason[no_loss] = 4, 0
    reason[usable & ~no_loss] = 6
    conflict = usable & (((loss >= threshold + margin) & opposite) |
                         ((loss <= -threshold) & event_positive))
    reason[conflict] = 5
    credible &= ~conflict
    late_ok = np.isfinite(p[-k:]).all(0)
    temporal_ok = (days[post[-1]] - days[post[-k]] >= c.MIN_LATE_CONFIRM_SPAN_DAYS)
    for left, right in zip(post[-k:-1], post[-k+1:]):
        temporal_ok &= days[right] - days[left] <= c.MAX_CONFIRM_GAP_DAYS
        temporal_ok &= c.ALLOW_OVERLAPPING_CONFIRMATION_PAIRS or obs[left].end <= obs[right].start
    late_ok &= temporal_ok
    reason[credible & ~late_ok] = 7
    late_predictions = predictions[-k:]
    agreement = (np.abs(p[-k:] - late_predictions) <= c.OBSERVATION_DISAGREEMENT_TOLERANCE).all(0)
    ready = credible & late_ok & agreement
    reason[credible & late_ok & ~agreement] = 8
    threshold_return = baseline - tol
    modeled_return = (late_predictions >= threshold_return).all(0)
    observed_return = (p[-k:] >= threshold_return).all(0)
    recovered = ready & modeled_return & observed_return
    deficit_supported = (deficit > tol + margin) & (baseline - late > tol - margin)
    observed_progress = late - early
    partial = (ready & ~recovered & deficit_supported &
               (progress > c.MIN_RECOVERY_IMPROVEMENT + margin) &
               (observed_progress >= c.MIN_RECOVERY_IMPROVEMENT - c.EFFECT_DECISION_MARGIN))
    degraded = (ready & ~recovered & deficit_supported &
                (progress <= c.MIN_RECOVERY_IMPROVEMENT - margin) &
                (observed_progress <= c.MIN_RECOVERY_IMPROVEMENT + margin))
    for mask, code in ((recovered, 1), (partial, 2), (degraded, 3)):
        classes[mask], reason[mask] = code, 0
    crossing = np.full(count, np.nan)
    candidate = np.full(count, np.nan)
    np.divide(threshold_return - (b[0] + b[2]), slope, out=candidate, where=slope > 0)
    bracketed = recovered & (candidate >= first_day) & (candidate <= last_day)
    crossing[bracketed] = candidate[bracketed]
    fraction = np.full(count, np.nan)
    np.divide(progress, early_deficit, out=fraction, where=credible & (early_deficit > 0))
    metrics = dict(baseline=baseline, modeled_level_loss=loss, observed_event_loss=event_loss,
                   post_slope=slope, modeled_initial=initial, modeled_final=final,
                   baseline_deficit=deficit, modeled_progress=progress, recovery_fraction=fraction,
                   crossing_days=crossing, adaptive_loss_threshold=threshold, rmse=fit.rmse,
                   level_se=fit.level_se, progress_se=fit.progress_se, final_prediction_se=fit.final_se,
                   pre_noise=fit.noise, pre_count=fit.pre_count,
                   post_count=fit.post_count, residual_df=fit.df, design_condition=fit.condition,
                   beta0=b[0], beta1=b[1], beta2=b[2], beta3=b[3])
    gates = dict(temporal_support=support, stable_model=finite, usable_reference=usable,
                 modeled_loss_above_threshold=loss >= threshold + margin,
                 event_support=event_positive, early_post_support=early_positive,
                 disturbance_conflict=conflict, credible_loss=credible,
                 no_loss_rule=no_loss, late_support=late_ok, late_model_agreement=agreement,
                 modeled_return=modeled_return, observed_return=observed_return)
    return Assessment(classes.reshape(shape), quality.reshape(shape), reason.reshape(shape),
                      {key: a.reshape(shape) for key, a in metrics.items()}, fit,
                      {key: a.reshape(shape) for key, a in gates.items()})
