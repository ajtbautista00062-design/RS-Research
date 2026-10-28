"""Conservative spatial refinement: candidates must independently support a class.

No classification thresholds are lowered. A regional label may change only where
the member pixel's own valid ITSA trajectory already supports the proposed class.
Consequently ambiguous raw pixel fits never acquire a class by proximity alone.
"""
import numpy as np


SIGNATURES = ("baseline", "modeled_level_loss", "modeled_initial", "modeled_final",
              "modeled_progress", "baseline_deficit")
CHANGE_REASONS = {0: "unchanged", 1: "own pixel ITSA agrees with connected neighboring seed",
                  255: "NoData/outside"}


def evidence_score(result, c):
    """Minimum of completeness, precision and decision eligibility, in [0,1].

    Precision = 1/(1 + 2*max(level_SE, progress_SE, final_SE)/effect_scale).
    effect_scale=max(loss threshold, recovery tolerance, minimum improvement).
    This OLS-conditional heuristic is not a calibrated probability.
    """
    m = result.metrics
    completeness = np.minimum(m["pre_count"] / 8., m["post_count"] / 8.)
    uncertainty = np.maximum.reduce([m["level_se"], m["progress_se"], m["final_prediction_se"]])
    scale = np.maximum(m["adaptive_loss_threshold"], max(c.RECOVERY_TOLERANCE, c.MIN_RECOVERY_IMPROVEMENT))
    precision = 1. / (1. + c.UNCERTAINTY_MULTIPLIER * uncertainty / scale)
    score = np.minimum(completeness, precision)
    eligible = np.isin(result.classes, (1, 2, 3, 4)) & (result.reasons == 0)
    return np.where(eligible & np.isfinite(score), np.clip(score, 0, 1), 0).astype("float32")


def shift(a, dr, dc, fill=0):
    """Non-wrapping neighbor lookup, also valid for tiny test arrays."""
    out = np.full_like(a, fill)
    h, w = a.shape
    if abs(dr) >= h or abs(dc) >= w:
        return out
    out[max(0, dr):min(h, h+dr), max(0, dc):min(w, w+dc)] = a[
        max(0, -dr):min(h, h-dr), max(0, -dc):min(w, w-dc)]
    return out


def refine(original, original_score, raw, strata, c, pixel_size):
    """Synchronous bounded four-neighbor growth, no dependence on iteration order.

    All signature comparisons use unmodified pixel model metrics. Strata encode
    supplied land cover and physical barriers; zero and NoData block every path.
    Caller reads a radius-sized halo; only the core is written, removing tile seams.
    """
    result = original.copy()
    score = np.asarray(original_score).copy()
    reason = np.where(original > 0, 0, 255).astype("uint8")
    if not (c.ENABLE_SPATIAL_REGION_GROWTH and c.ENABLE_NEIGHBOR_CLASS_ASSIGNMENT):
        return result, score, reason
    threshold = c.EVIDENCE_PROTECTION_SCORE
    own_score = evidence_score(raw, c)
    seeds = np.isin(original, (1, 2, 3)) & (score >= threshold) & (strata > 0)
    protected = np.isin(original, (1, 2, 3, 4)) & (score >= threshold)
    candidate = ((original > 0) & (original != raw.classes) & (strata > 0) & ~protected &
                 np.isin(raw.classes, (1, 2, 3)) & (own_score >= threshold))
    active = seeds.copy()
    distance = np.where(seeds, 0., np.inf)
    radius = c.NEIGHBOR_RADIUS_M
    steps = int(np.ceil(radius / min(pixel_size)))
    tolerance = c.NEIGHBOR_SIGNATURE_TOLERANCE
    for _ in range(steps):
        votes = np.zeros(original.shape, dtype="uint8")
        possible = votes.copy()
        best = np.zeros(original.shape, dtype="float32")
        best_distance = np.full(original.shape, np.inf)
        for dr, dc, length in ((1, 0, pixel_size[1]), (-1, 0, pixel_size[1]),
                               (0, 1, pixel_size[0]), (0, -1, pixel_size[0])):
            linked = shift(active, dr, dc) & (strata == shift(strata, dr, dc))
            proposed_distance = shift(distance, dr, dc, np.inf) + length
            linked &= proposed_distance <= radius
            possible += linked
            delta = np.maximum.reduce([np.abs(raw.metrics[k] - shift(raw.metrics[k], dr, dc, np.nan))
                                       for k in SIGNATURES])
            similarity = np.clip(1. - delta / tolerance, 0., 1.)
            same = linked & (shift(result, dr, dc) == raw.classes) & np.isfinite(similarity)
            votes += same & (similarity >= threshold)
            local = np.minimum.reduce([own_score, similarity,
                                       np.clip(1.-.2*proposed_distance/radius, 0, 1)])
            best = np.maximum(best, np.where(same, local, 0))
            best_distance = np.minimum(best_distance, np.where(same & (local >= threshold), proposed_distance, np.inf))
        agreement = votes / np.maximum(possible, 1)
        decision_score = np.minimum(best, agreement)
        accepted = candidate & ~active & (decision_score >= threshold) & (agreement >= c.NEIGHBOR_MIN_AGREEMENT)
        if not accepted.any():
            break
        changed = accepted & (result != raw.classes)
        result[accepted] = raw.classes[accepted]
        score[accepted] = decision_score[accepted]
        reason[changed] = 1
        distance[accepted] = best_distance[accepted]
        active |= accepted
    return result, score, reason
