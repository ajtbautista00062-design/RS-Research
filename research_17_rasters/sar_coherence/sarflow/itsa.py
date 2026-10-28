"""Actual-date segmented OLS, grouped by missingness; no temporal interpolation."""
from dataclasses import dataclass
import warnings
import numpy as np
from .temporal import day_offsets


@dataclass
class Fit:
    beta: np.ndarray
    level_se: np.ndarray
    rmse: np.ndarray
    df: np.ndarray
    condition: np.ndarray
    pre_count: np.ndarray
    post_count: np.ndarray
    baseline: np.ndarray
    noise: np.ndarray
    pre_change: np.ndarray
    progress_se: np.ndarray
    final_se: np.ndarray


def median(a, axis=0):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(a, axis=axis)


def fit_itsa(values, obs, c):
    """Coefficients in days centered on EVENT_ORIGIN; event pair excluded from OLS.

    Numerical scaling of time is undone before reporting coefficients. Classical
    SEs are descriptive conditional-on-OLS screens; shared SLCs may correlate errors.
    """
    y = np.asarray(values, dtype=float).reshape(len(obs), -1)
    pre = np.array([o.phase == "pre" for o in obs])
    post = np.array([o.phase == "post" for o in obs])
    use = pre | post
    days = day_offsets(obs, c.EVENT_ORIGIN)
    scale = max(float(np.max(np.abs(days[use]))), 1.)
    t = days / scale
    x = np.column_stack((np.ones(len(obs)), t, post, t * post))
    n = y.shape[1]
    beta = np.full((4, n), np.nan)
    se, rmse, df, condition, noise, change = [np.full(n, np.nan) for _ in range(6)]
    progress_se, final_se = np.full(n, np.nan), np.full(n, np.nan)
    valid = np.isfinite(y) & use[:, None]
    patterns, membership = np.unique(valid.T, axis=0, return_inverse=True)
    for group, pattern in enumerate(patterns):
        ix = np.flatnonzero(membership == group)
        a = x[pattern]
        degrees = int(pattern.sum()) - 4
        df[ix] = degrees
        if degrees < 1 or np.linalg.matrix_rank(a) < 4:
            continue
        cond = np.linalg.cond(a)
        condition[ix] = cond
        if cond > c.MAX_DESIGN_CONDITION:
            continue
        b = np.linalg.lstsq(a, y[pattern][:, ix], rcond=None)[0]
        residual = y[pattern][:, ix] - a @ b
        variance = np.sum(residual**2, axis=0) / degrees
        covariance = np.linalg.inv(a.T @ a)
        se[ix] = np.sqrt(np.maximum(0, variance * covariance[2, 2]))
        post_dates = t[pattern & post]
        if len(post_dates):
            delta = np.array([0., 1., 0., 1.]) * np.ptp(post_dates)
            end = np.array([1., post_dates.max(), 1., post_dates.max()])
            progress_se[ix] = np.sqrt(np.maximum(0, variance * (delta @ covariance @ delta)))
            final_se[ix] = np.sqrt(np.maximum(0, variance * (end @ covariance @ end)))
        rmse[ix] = np.sqrt(np.mean(residual**2, axis=0))
        beta[:, ix] = b / np.array([1., scale, 1., scale])[:, None]
        p = pattern & pre
        if p.sum() >= 3:
            residual_pre = y[p][:, ix] - x[p, :2] @ b[:2]
            noise[ix] = 1.4826 * median(np.abs(residual_pre - median(residual_pre)))
            change[ix] = b[1] * np.ptp(t[p])
    return Fit(beta, se, rmse, df, condition, np.isfinite(y[pre]).sum(0),
               np.isfinite(y[post]).sum(0), median(y[pre]), noise, change, progress_se, final_se)
