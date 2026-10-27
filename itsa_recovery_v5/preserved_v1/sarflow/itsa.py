"""Four-parameter ITSA: event-straddling pair is displayed but excluded."""
import numpy as np
from .temporal import timestamp

COEFFICIENTS = ("pre_level_at_event", "pre_slope_per_day", "post_level_change_at_event",
                "post_slope_change_per_day")
DIAGNOSTICS = ("post_slope_per_day", "final_minus_pretrend_prediction",
               "valid_model_observation_count", "residual_rmse", "model_r_squared")
MODEL = "y = b0 + b1*t + b2*POST + b3*t*POST; event-straddling pair excluded"
LIMITATION = ("Pilot estimates are exploratory. Three PRE and three POST observations "
              "estimate four parameters with only two residual degrees of freedom. "
              "OLS uncertainty assumes independent errors; shared SAR acquisitions may violate "
              "this assumption. Full-series analysis and external validation remain necessary.")

def design(obs, c):
    phases = [o.phase for o in obs]
    if phases.count("pre") < 3 or phases.count("post") < 3 or phases.count("event") != 1:
        raise ValueError("Require at least 3 unique PRE, 1 event-straddling and 3 POST pairs.")
    t = np.array([(o.reference-timestamp(c.EVENT_ORIGIN)).total_seconds()/86400 for o in obs])
    post = np.array([p == "post" for p in phases], dtype=float)
    x = np.column_stack((np.ones(len(obs)), t, post, t*post))
    x[np.array(phases) == "event"] = 0  # exclusion marker, never a fitted event value
    if np.linalg.matrix_rank(x) != 4:
        raise ValueError("Rank-deficient ITSA design.")
    return x

def band_names(obs):
    return tuple(f"observed_{o.phase}_{o.start:%Y%m%d}_{o.end:%Y%m%d}" for o in obs) + COEFFICIENTS + DIAGNOSTICS

def fit_complete(values, x):
    use = x[:, 0] == 1
    return np.linalg.pinv(x[use]) @ values[use]

def raster_bands(stack, x):
    use = x[:, 0] == 1
    count = np.isfinite(stack[use]).sum(axis=0)
    complete = count == int(use.sum())
    result = np.full((9, *stack.shape[1:]), np.nan, dtype="float32")
    result[6] = np.where(count > 0, count, np.nan)
    if complete.any():
        y = stack[use][:, complete].astype("float64")
        beta = np.linalg.pinv(x[use]) @ y
        residual = y-x[use]@beta
        ss = np.sum(residual**2, axis=0)
        total = np.sum((y-y.mean(axis=0))**2, axis=0)
        result[:4, complete] = beta
        result[4, complete] = beta[1]+beta[3]
        result[5, complete] = y[-1]-(beta[0]+beta[1]*x[use][-1,1])
        result[7, complete] = np.sqrt(ss/use.sum())
        r2 = np.full_like(ss, np.nan)
        np.divide(ss, total, out=r2, where=total > 1e-20)
        result[8, complete] = 1-r2  # undefined for constant series
    return np.concatenate((stack, result), axis=0)
