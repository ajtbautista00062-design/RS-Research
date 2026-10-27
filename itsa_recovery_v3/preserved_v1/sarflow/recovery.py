"""Operational baseline recovery; independent of the ITSA counterfactual."""
import warnings
import numpy as np

LABELS = {1: "Coherence recovered", 2: "Partially recovered",
          3: "Not recovered / further degraded", 4: "No detectable event loss",
          5: "Indeterminate - insufficient POST support"}
COLORS = {0: (0, 0, 0, 0), 1: (0, 158, 115, 255), 2: (230, 171, 2, 255),
          3: (213, 94, 0, 255), 4: (0, 114, 178, 255), 5: (150, 150, 150, 255)}

def settings(c):
    return dict(reference="median of valid PRE observations; not ITSA counterfactual",
                loss_floor=c.LOSS_MIN, spread_multiplier=c.LOSS_SPREAD_MULTIPLIER,
                spread_floor=c.PRE_SPREAD_FLOOR, return_tolerance=c.RETURN_TOLERANCE,
                minimum_pre=c.CLASS_MIN_PRE, minimum_post=c.CLASS_MIN_POST,
                end_window=c.END_POST_COUNT, confirmation=c.CONFIRM_POST_COUNT,
                fraction="(median of last configured POST pairs - event)/(baseline - event)",
                interpretation="Exploratory operational coherence classes, not significance or reconstruction",
                comparability_confirmed=c.ACQUISITION_COMPARABILITY_CONFIRMED)

def validate(c):
    for name in ("CLASS_MIN_PRE", "CLASS_MIN_POST", "END_POST_COUNT", "CONFIRM_POST_COUNT"):
        value = getattr(c, name)
        if type(value) is not int or value < 1:
            raise ValueError(f"{name} must be a positive integer.")
    if c.CONFIRM_POST_COUNT < 2 or c.END_POST_COUNT < c.CONFIRM_POST_COUNT:
        raise ValueError("Use at least two consecutive confirmations within the endpoint window.")
    if c.CLASS_MIN_PRE < 3:
        raise ValueError("Use at least three PRE observations for the exploratory baseline.")
    if not (0 < c.LOSS_MIN <= 1 and 0 <= c.RETURN_TOLERANCE < c.LOSS_MIN
            and 0 < c.PRE_SPREAD_FLOOR <= 1 and c.LOSS_SPREAD_MULTIPLIER > 0
            and np.isfinite(c.LOSS_SPREAD_MULTIPLIER)):
        raise ValueError("Invalid loss/spread/return settings; return tolerance must be smaller than minimum loss.")

def baseline_recovery(stack, obs, c):
    """No temporal gap filling, baseline cutoff, or clipping of recovery fraction."""
    pre = stack[[i for i, o in enumerate(obs) if o.phase == "pre"]].astype("float64")
    post = stack[[i for i, o in enumerate(obs) if o.phase == "post"]].astype("float64")
    event = stack[next(i for i, o in enumerate(obs) if o.phase == "event")].astype("float64")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        baseline = np.nanmedian(pre, axis=0)
        spread = 1.4826*np.nanmedian(np.abs(pre-baseline), axis=0)
        endpoint = np.nanmedian(post[-c.END_POST_COUNT:], axis=0)
    loss = baseline-event
    threshold = np.maximum(c.LOSS_MIN, c.LOSS_SPREAD_MULTIPLIER*np.maximum(spread, c.PRE_SPREAD_FLOOR))
    reference_valid = (np.isfinite(event) & np.isfinite(baseline)
                       & (np.isfinite(pre).sum(axis=0) >= c.CLASS_MIN_PRE))
    # Small tolerance only protects float32 representation at exact decision boundaries.
    eps = 1e-7
    detected = reference_valid & (loss >= threshold-eps)
    codes = np.zeros(event.shape, dtype="uint8")
    codes[reference_valid & ~detected] = 4
    endpoint_valid = (np.isfinite(post[-c.END_POST_COUNT:]).all(axis=0)
                      & (np.isfinite(post).sum(axis=0) >= c.CLASS_MIN_POST))
    eligible = detected & endpoint_valid
    codes[detected & ~endpoint_valid] = 5
    fraction = np.full(event.shape, np.nan, dtype="float64")
    np.divide(endpoint-event, loss, out=fraction, where=eligible)
    codes[eligible] = 3
    codes[eligible & (endpoint > event+eps)] = 2
    returned = (post[-c.CONFIRM_POST_COUNT:] >= baseline-c.RETURN_TOLERANCE-eps).all(axis=0)
    codes[eligible & returned] = 1
    return codes, fraction.astype("float32")
