"""Effect-size ITSA decisions with separate uncertainty; no p-value class gateway."""
import numpy as np
from operational_support import quantities,evidence_values


def add_evidence(result,stack,obs,origin,c,inside=None,filled=None):
    q=quantities(result,stack,obs,origin,c,inside)
    values,_=evidence_values(q,stack,c,result.classes)
    if filled is not None:
        values['estimated_observation_count']=np.where(filled,np.nan,0)
    result.metrics.update({k:np.where(q.inside,np.asarray(v,dtype='float32'),np.nan) for k,v in values.items()})
    return result


def supported(q,stack,c):
    any_data=np.isfinite(stack).any(axis=0)
    reason=np.ones(q.baseline.shape,'uint8')
    reason[q.inside&~any_data]=9
    pre=q.inside&any_data&(q.npre>=c.MIN_PRE)&np.isfinite(q.baseline)
    reason[q.inside&any_data&~pre]=2
    post=pre&(q.npost>=c.MIN_POST)
    reason[pre&~post]=4
    event_ok=np.isfinite(q.event_loss)|~np.full(q.baseline.shape,c.REQUIRE_EVENT)
    reason[post&~event_ok]=3
    late=np.isfinite(q.late_observed).all(axis=0)
    reason[post&event_ok&~late]=19
    fit=post&event_ok&late&q.model.valid
    reason[post&event_ok&late&~q.model.valid]=5
    physical=np.any((q.model.predicted[q.post]<-c.MAX_PREDICTION_EXCESS)|
                    (q.model.predicted[q.post]>1+c.MAX_PREDICTION_EXCESS),axis=0)
    finite=np.isfinite(q.model.beta).all(axis=0)&np.isfinite(q.model.predicted[q.post]).all(axis=0)
    reason[fit&(~finite|physical)]=18
    rmse=fit&(q.model.rmse>c.MODEL_RMSE_FAILURE)
    reason[rmse]=17
    zero=q.inside&np.all(stack==0,axis=0)&c.ZERO_UNRESOLVED
    reason[zero]=10
    eligible=fit&finite&~physical&~rmse&~zero
    reason[eligible]=16
    return eligible,reason,any_data


def disturbance(q,eligible,c,reason):
    loss=eligible&(q.loss>=q.threshold+c.EFFECT_MARGIN)
    absence=eligible&(q.loss<=q.threshold-c.EFFECT_MARGIN)
    event_support=q.event_loss>=c.EVENT_EVIDENCE_MIN
    event_opposes=q.event_loss<=-c.EVENT_EVIDENCE_MIN
    early_support=(q.early_deficit>=c.EARLY_LOSS_MIN)&(q.observed_early_deficit>=c.EARLY_LOSS_MIN)
    contrary=event_support|early_support
    reason[eligible&~loss&~absence]=12
    reason[loss&~(event_support|early_support)]=13
    reason[loss&event_opposes]=14
    conflict=absence&contrary
    reason[conflict]=24
    # Fixed-baseline disturbance without level-change agreement is unresolved.
    no_loss=absence&~contrary&(q.early_deficit<c.EARLY_LOSS_MIN)&(q.observed_early_deficit<c.EARLY_LOSS_MIN)
    reason[absence&~no_loss&~conflict]=14
    return loss&(event_support|early_support)&~event_opposes,no_loss,event_support


def recovery_masks(q,loss,c,reason):
    model_return=(q.late_model>=q.target).all(axis=0)
    observed_return=(q.late_observed>=q.target).all(axis=0)
    below=q.final<q.target-c.EFFECT_MARGIN
    obs_below=(q.late_median<q.target)&(q.late_observed[-1]<q.target+c.OBS_CONFIRM_TOLERANCE)
    contradictory=loss&(np.max(np.abs(q.late_observed-q.late_model),axis=0)>c.MAX_OBSERVATION_DISAGREEMENT)
    contradictory|=loss&((model_return&~observed_return)|
        (below&~obs_below))
    reason[contradictory]=15
    consistent=loss&~contradictory
    recovered=consistent&model_return&observed_return
    partial=(consistent&below&obs_below&(q.progress>c.MIN_PROGRESS+c.EFFECT_MARGIN)&
             (q.observed_progress>=c.MIN_PROGRESS-c.OBS_CONFIRM_TOLERANCE))
    failed=(consistent&below&obs_below&(q.progress<=c.MIN_PROGRESS-c.EFFECT_MARGIN)&
            (q.observed_progress<=c.MIN_PROGRESS+c.OBS_CONFIRM_TOLERANCE))
    reason[consistent&~(recovered|partial|failed)]=16
    return recovered,partial,failed


def crossing(q,stack,recovered):
    candidate=np.full(q.baseline.shape,np.nan)
    np.divide(q.target-q.model.beta[0]-q.model.beta[2],q.slope,out=candidate,where=q.model.valid&(q.slope>0))
    lower=np.full_like(candidate,np.nan);upper=lower.copy()
    for left,right in zip(q.post[:-1],q.post[1:]):
        bracket=((candidate>=q.time[left])&(candidate<=q.time[right])&
                 (stack[left]<q.target)&(stack[right]>=q.target))
        lower[bracket]=q.time[left];upper[bracket]=q.time[right]
    supported=recovered&np.isfinite(lower)&(candidate>=q.start_time)&(candidate<=q.time[q.post[-1]])
    status=np.where(supported,1,np.where(recovered&(q.first_model>=q.target)&(q.first_observed>=q.target),2,
        np.where(recovered,3,np.where(candidate>q.time[q.post[-1]],4,0))))
    return dict(candidate_crossing_days=candidate,recovery_duration_days=np.where(supported,candidate,np.nan),
        crossing_lower_days=np.where(supported,lower,np.nan),crossing_upper_days=np.where(supported,upper,np.nan),
        crossing_status=status)


def classify(result,stack,obs,origin,c,inside=None,filled=None):
    q=quantities(result,stack,obs,origin,c,inside)
    guarded=result.classes.copy()
    eligible,reason,any_data=supported(q,stack,c)
    loss,no_loss,event_support=disturbance(q,eligible,c,reason)
    recovered,partial,failed=recovery_masks(q,loss,c,reason)
    classes=np.where(q.inside&any_data,5,0).astype('uint8')
    for mask,code in ((no_loss,4),(recovered,1),(partial,2),(failed,3)):
        classes[mask]=code;reason[mask]=0
    values,wide=evidence_values(q,stack,c,guarded)
    decided=np.isin(classes,[1,2,3,4])
    provisional=decided&((guarded!=classes)|wide)
    values['evidence_status']=np.where(decided,np.where(provisional,2,1),0)
    values.update(pre_baseline=q.baseline,modeled_first_post=q.first_model,modeled_final_post=q.final,
        baseline_deficit_final=q.baseline-q.final,modeled_recovery_progress=q.progress,
        loss_detection_boundary=q.threshold,final_guard_halfwidth=q.guard[q.post[-1]])
    values.update(crossing(q,stack,recovered))
    if filled is not None:
        values['estimated_observation_count']=np.where(filled,np.nan,0)
    quality=result.quality & np.uint16(1|2|4|128|4096|8192)
    def flag(mask,bit):
        quality[q.inside&mask]|=bit
    flag(wide,32768);flag(provisional,16384)
    flag(np.isin(reason,[2,3,4,5,19]),64)
    flag(np.isin(reason,[14,24]),16);flag(reason==15,32)
    flag(np.isin(reason,[12,13,16]),8)
    flag((q.model.rmse>c.MAX_RMSE)|(q.spread>c.MAX_PRE_SPREAD)|np.isin(reason,[17,18]),512)
    flag(loss&~event_support,2048)
    flag(values['crossing_status']==2,1024);flag(values['crossing_status']==4,256)
    quality[~q.inside]=65535
    result.classes=classes;result.exclusion=reason;result.quality=quality
    result.metrics.update({k:np.where(q.inside&np.isfinite(v),v,np.nan).astype('float32') for k,v in values.items()})
    return result
