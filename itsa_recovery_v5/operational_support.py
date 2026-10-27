"""Shared fitted quantities and uncertainty, independent of operational class gates."""
from types import SimpleNamespace
import warnings
import numpy as np
from scipy.stats import t as student_t
from model import design,contrast_se


def quantities(result,stack,obs,origin,c,inside):
    pre=[j for j,o in enumerate(obs) if o.phase=='pre']
    post=[j for j,o in enumerate(obs) if o.phase=='post']
    event=next(j for j,o in enumerate(obs) if o.phase=='event')
    _,_,time=design(obs,origin)
    if inside is None:
        inside=np.ones(stack.shape[1:],bool)
    m=result.model
    first_index=np.argmax(np.isfinite(stack[post]),axis=0)
    first=lambda a:np.take_along_axis(a,first_index[None],axis=0)[0]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        baseline=np.nanmedian(stack[pre],axis=0)
        spread=1.4826*np.nanmedian(np.abs(stack[pre]-baseline),axis=0)
        late_median=np.nanmedian(stack[post[-c.CONFIRM_POST:]],axis=0)
    npre=np.isfinite(stack[pre]).sum(axis=0);npost=np.isfinite(stack[post]).sum(axis=0)
    start_time=np.asarray(time[post])[first_index]
    dt=time[post[-1]]-start_time
    slope=m.beta[1]+m.beta[3]
    critical=student_t.ppf(1-c.DETECTION_ALPHA,np.maximum(m.df,1))
    baseline_guard=np.maximum(c.BASELINE_GUARD_FLOOR,c.GUARD_MULTIPLIER*spread/np.sqrt(np.maximum(npre,1)))
    guard=critical*m.se+baseline_guard
    step_se=contrast_se(m,np.array([0.,0.,1.,0.]))
    progress_guard=critical*contrast_se(m,np.array([0.,1.,0.,1.]))*dt
    threshold=np.maximum(c.LOSS_MIN,getattr(c,'LOSS_NOISE_MULTIPLIER',.5)*m.pre_noise)
    if getattr(c,'EVIDENCE_MODE','guarded')=='guarded':
        threshold=np.full(stack.shape[1:],c.LOSS_MIN)
    target=baseline-c.RETURN_TOLERANCE
    late=post[-c.CONFIRM_POST:]
    first_model=first(m.predicted[post]);first_observed=first(stack[post])
    return SimpleNamespace(pre=pre,post=post,event=event,time=time,inside=inside,model=m,
        baseline=baseline,spread=spread,npre=npre,npost=npost,critical=critical,guard=guard,
        loss=-m.beta[2],step_se=step_se,step_guard=critical*step_se,threshold=threshold,
        first_model=first_model,first_observed=first_observed,start_time=start_time,
        final=m.predicted[post[-1]],slope=slope,progress=slope*dt,progress_guard=progress_guard,
        target=target,late=late,late_median=late_median,
        late_observed=stack[late],late_model=m.predicted[late],
        event_loss=baseline-stack[event],early_deficit=baseline-first_model,
        observed_early_deficit=baseline-first_observed,
        observed_progress=late_median-first_observed)


def evidence_values(q,stack,c,guarded_class):
    width=q.guard[q.post[-1]]
    wide=(((q.loss-q.step_guard<=q.threshold)&(q.loss+q.step_guard>=q.threshold))|
          (width>max(c.RETURN_TOLERANCE,c.MIN_PROGRESS)))
    event_code=np.where(~np.isfinite(q.event_loss),0,np.where(q.event_loss>=c.EVENT_EVIDENCE_MIN,1,
                        np.where(q.event_loss<=-c.EVENT_EVIDENCE_MIN,3,2)))
    values=dict(valid_pre_count=q.npre,valid_post_count=q.npost,fit_observation_count=q.model.n,
        residual_degrees_freedom=q.model.df,regression_rmse=q.model.rmse,
        regression_condition=q.model.condition,loss_effect_threshold=q.threshold,
        level_loss_guard_halfwidth=q.step_guard,progress_guard_halfwidth=q.progress_guard,
        early_post_deficit=q.early_deficit,event_agreement_code=event_code,
        late_confirmation_count=np.sum(q.late_observed>=q.target,axis=0),
        evidence_status=np.where(np.isin(guarded_class,[1,2,3,4]),1,0),
        guarded_reference_class=guarded_class,crossing_status=np.zeros(q.baseline.shape),
        estimated_observation_count=np.zeros(q.baseline.shape))
    return values,wide
