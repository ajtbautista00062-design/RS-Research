"""Exact bounded-region temporal medians and operational ITSA evidence scores."""
import warnings
import numpy as np
from recovery import assess


def aggregate(values, minimum_pixels, minimum_fraction):
    """values[date,pixel]; missing dates remain NaN, not interpolated observations."""
    counts=np.isfinite(values).sum(axis=1)
    fractions=counts/values.shape[1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        medians=np.nanmedian(values,axis=1)
        mad=1.4826*np.nanmedian(np.abs(values-medians[:,None]),axis=1)
    eligible=(counts>=minimum_pixels)&(fractions>=minimum_fraction)
    return np.where(eligible,medians,np.nan),counts,fractions,mad


def evidence_score(result,c):
    """0..1 reproducibility/guard-width score; explicitly NOT a calibrated probability.

    A finite model with narrower guards receives a higher score. Indeterminate
    and NoData always score zero; thresholds must be validated on independent data.
    """
    width=result.metrics["final_guard_halfwidth"]
    scale=max(c.RETURN_TOLERANCE,c.MIN_PROGRESS,c.LOSS_MIN)
    score=np.clip(1-np.nan_to_num(width,nan=np.inf)/scale,0,1)
    return np.where(np.isin(result.classes,[1,2,3,4]),score,0)


def classify_region(values, obs, origin, analysis, spatial, area, pixel_classes=None,
                    pixel_scores=None):
    # The absolute spatial support requirement cannot exceed actual membership.
    # Small objects retain their real series and explicit MMU/provisional flags.
    minimum=min(spatial.OBJECT_MIN_VALID_PIXELS,values.shape[1])
    series,counts,fractions,mad=aggregate(values,minimum,
                                        spatial.OBJECT_MIN_VALID_FRACTION)
    result=assess(series[:,None,None],obs,origin,analysis)
    code=int(result.classes.item()); score=float(evidence_score(result,analysis).item())
    flags=[]
    if area < spatial.MIN_REGION_AREA_M2:
        flags.append("below_mapping_unit")
    if np.any(np.isfinite(mad)&(mad>spatial.OBJECT_MAX_SPATIAL_MAD)):
        flags.append("heterogeneous_temporal_evidence")
        if getattr(spatial,'OBJECT_HETEROGENEITY_POLICY','strict')=='strict':
            code=5
    # Do not conceal small strong pixel evidence inside a spatial median.
    protected=False
    if pixel_classes is not None and pixel_scores is not None:
        strong=(pixel_scores>=spatial.CLASS_CONFIDENCE_THRESHOLD)&np.isin(pixel_classes,[1,2,3])
        if np.any(strong&(pixel_classes!=code)):
            code=5; flags.append("conflicting_high_confidence_pixel_reference")
        protected=bool(np.any(strong&(pixel_classes==code)))
    if code in (0,5):
        # A segmented object has PRE support: temporal insufficiency is indeterminate.
        code=5; score=0.
    decision_flags=(1 if 'below_mapping_unit' in flags else 0)|(2 if 'heterogeneous_temporal_evidence' in flags else 0)|(4 if 'conflicting_high_confidence_pixel_reference' in flags else 0)
    row=dict(class_code=code,raw_itsa_class=int(result.classes.item()),evidence_score=score,
        protected_high_confidence=protected,
        object_flags=";".join(flags),object_decision_flags=decision_flags,
        quality_flags=int(result.quality.item()),exclusion_reason=int(result.exclusion.item()),
        n_fit=int(result.model.n.item()),residual_df=float(result.model.df.item()),
        rmse=float(result.model.rmse.item()),b0=float(result.model.beta[0].item()),
        b1=float(result.model.beta[1].item()),b2=float(result.model.beta[2].item()),
        b3=float(result.model.beta[3].item()))
    if int(result.classes.item())==0:
        row['exclusion_reason']=20
        row['quality_flags']|=64
    if code==5 and int(result.classes.item()) not in (0,5):
        row['exclusion_reason']=22 if decision_flags&4 else 21
        row['quality_flags']|=32 if decision_flags&4 else 8
    row.update({name:float(value.item()) for name,value in result.metrics.items()})
    if code==5 and row['exclusion_reason'] in (2,4):
        row['exclusion_reason']=20
    if code in (1,2,3,4) and decision_flags&3:
        row['quality_flags']|=16384
        row['evidence_status']=2
    if code==5:
        row['evidence_status']=0
    for j in range(len(obs)):
        row.update({f"date_{j+1}_median":float(series[j]),f"date_{j+1}_count":int(counts[j]),
                    f"date_{j+1}_fraction":float(fractions[j]),f"date_{j+1}_MAD":float(mad[j])})
    return row
