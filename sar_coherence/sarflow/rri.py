"""Pure numerical ITSA final-observation RRI; no spatial propagation or clipping."""
import numpy as np

INTERPRETATION = {0:'NoData', 1:'Further coherence degradation', 2:'No measurable recovery',
                  3:'Partial coherence recovery', 4:'Return to modeled baseline',
                  5:'Coherence exceeding modeled baseline'}
STATUS = {0:'Outside study area / structurally unavailable', 1:'Valid RRI',
          2:'No detectable coherence disturbance', 3:'Insufficient observations',
          4:'Invalid or near-zero disturbance denominator', 5:'Invalid fitted model result',
          6:'Unresolved or contradictory temporal evidence'}


def calculate(baseline, event, final, detected, supported, model_valid, inside,
              loss_threshold, minimum=.01, zero_tolerance=.05, one_tolerance=.05,
              no_loss=None):
    """All inputs broadcast to a common shape; output undefined values are NaN.

    Caller supplies existing raw ITSA disturbance/support decisions, never an
    enhanced neighbor label. Model estimates at intervention are beta0 and
    beta0+beta2. Threshold comparison is strict, preventing unstable division.
    """
    if minimum <= 0 or not 0 <= zero_tolerance < .5 or not 0 <= one_tolerance < .5:
        raise ValueError('Require positive denominator floor and disjoint tolerances <0.5.')
    b,e,f = np.broadcast_arrays(np.asarray(baseline,float),np.asarray(event,float),np.asarray(final,float))
    denominator, recovery = b-e, f-e
    status = np.zeros(b.shape,dtype='uint8')
    status[np.asarray(inside,bool)] = 6
    support = np.asarray(inside,bool) & np.asarray(supported,bool)
    status[np.asarray(inside,bool) & ~np.asarray(supported,bool)] = 3
    finite = np.isfinite(b)&np.isfinite(e)&np.isfinite(f)&np.asarray(model_valid,bool)
    finite &= (b>=0)&(b<=1)&(e>=0)&(e<=1)&(f>=0)&(f<=1)
    status[support & ~finite] = 5
    eligible = support & finite
    if no_loss is not None:
        status[eligible & np.asarray(no_loss,bool)] = 2
    disturbed = eligible & np.asarray(detected,bool)
    good_denominator = np.isfinite(denominator)&np.isfinite(loss_threshold)&(denominator > np.maximum(minimum,loss_threshold))
    status[disturbed & ~good_denominator] = 4
    valid = disturbed & good_denominator
    value = np.full(b.shape,np.nan)
    np.divide(recovery,denominator,out=value,where=valid)
    valid &= np.isfinite(value)&(np.abs(value)<=np.finfo('float32').max)
    status[disturbed & good_denominator & ~valid] = 5
    status[valid] = 1
    value[~valid] = np.nan
    eps=16*np.finfo(float).eps*(1+np.abs(value))
    classes = np.zeros(b.shape,dtype='uint8')
    classes[valid & (value < -zero_tolerance-eps)] = 1
    classes[valid & (np.abs(value) <= zero_tolerance+eps)] = 2
    classes[valid & (value > zero_tolerance+eps) & (value < 1-one_tolerance-eps)] = 3
    classes[valid & (np.abs(value-1) <= one_tolerance+eps)] = 4
    classes[valid & (value > 1+one_tolerance+eps)] = 5
    return value,classes,status,np.where(valid,denominator,np.nan),np.where(valid,recovery,np.nan)
