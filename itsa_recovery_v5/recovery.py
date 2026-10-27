"""ITSA-derived baseline recovery. Every assessed class requires an eligible OLS fit."""
from dataclasses import dataclass
import warnings
import numpy as np
from model import fit,design,contrast_se
from scipy.stats import t as student_t

LABELS={0:"NoData",1:"Coherence recovered",2:"Partially recovered",
        3:"Not recovered / further degraded",4:"No detectable event-associated coherence loss",
        5:"Indeterminate"}
COLORS={0:(0,0,0,0),1:(0,158,115,255),2:(230,171,2,255),3:(213,94,0,255),
        4:(0,114,178,255),5:(150,150,150,255)}
QUALITY={1:"Unresolved zero semantics",2:"Unverified acquisition comparability",
    4:"Fewer than 8 PRE or 8 POST dates",8:"Uncertainty/threshold ambiguity",
    16:"Conflicting event/level loss evidence",32:"Model-observation disagreement",
    64:"Insufficient support",128:"Spatial estimates used",256:"Crossing after last POST",
    512:"Unstable/poor model",1024:"Crossing before first POST",
    2048:"Event pair neutral rather than confirming fitted loss",
    4096:"Detection calibration unvalidated",8192:"Unverified all-zero sequence",
    16384:"Provisional operational class; not guard-supported",
    32768:"Wide effect/return uncertainty; reported separately"}
EXCLUSION={0:"Assessed (classes 1-4)",1:"Outside boundary",2:"Insufficient PRE",
    3:"Missing event",4:"Insufficient POST",5:"Model fit unavailable",
    6:"Uncertain decision",7:"Conflicting evidence",8:"Unstable model",
    9:"Structurally unavailable all observations",10:"Unverified all-zero sequence",
    11:"Return not bracketed within observed POST interval",
    12:"Modeled loss within operational decision margin",13:"Modeled loss lacks corroboration",
    14:"Event/early POST evidence contradicts modeled disturbance",15:"Late observations contradict fitted outcome",
    16:"Recovery improvement/return genuinely unresolved",17:"Hard RMSE failure",
    18:"Nonfinite or physically invalid fitted trajectory",19:"Missing required late POST confirmation",
    20:"Object spatial observation coverage insufficient",21:"Object heterogeneity rejected under strict policy",
    22:"Strong pixel evidence contradicts object fit",23:"Optional MMU deferral",
    24:"No-loss hypothesis contradicted by event/early POST observations"}
CORE_METRICS=("pre_baseline","modeled_first_post","modeled_final_post","baseline_deficit_final",
    "modeled_recovery_progress","candidate_crossing_days","recovery_duration_days",
    "crossing_lower_days","crossing_upper_days","final_guard_halfwidth","event_loss","fitted_level_loss","post_slope",
    "pre_noise","level_change_se","loss_detection_boundary","observed_event_change")
EVIDENCE_METRICS=("valid_pre_count","valid_post_count","fit_observation_count","residual_degrees_freedom",
    "regression_rmse","regression_condition","loss_effect_threshold","level_loss_guard_halfwidth",
    "progress_guard_halfwidth","early_post_deficit","event_agreement_code","late_confirmation_count",
    "evidence_status","guarded_reference_class","crossing_status","estimated_observation_count")
METRICS=CORE_METRICS+EVIDENCE_METRICS

def metric_unit(name):
    if name.endswith('_days'):
        return 'days since event origin'
    if name=='post_slope':
        return 'coherence/day'
    if name.endswith('_count') or name=='residual_degrees_freedom':
        return 'count'
    if name in ('event_agreement_code','evidence_status','guarded_reference_class','crossing_status'):
        return 'code'
    if name=='regression_condition':
        return 'dimensionless'
    return 'coherence'

EVIDENCE_TAGS={'evidence_status':{0:'unresolved/unavailable',1:'operational class agrees with guarded rule; conditional',2:'provisional operational estimate'},
    'event_agreement_code':{0:'missing event',1:'supports loss',2:'neutral',3:'opposes loss'},
    'crossing_status':{0:'no supported return',1:'observed-period bracketed crossing',
                       2:'baseline reached by first observed POST; timing censored',
                       3:'return supported but crossing unobservable',4:'candidate after final POST'}}

@dataclass
class Assessment:
    classes: np.ndarray
    quality: np.ndarray
    exclusion: np.ndarray
    metrics: dict
    model: object

def _assess_guarded(stack,obs,origin,c,inside=None,filled=None):
    shape=stack.shape[1:]
    if inside is None:
        inside=np.ones(shape,dtype=bool)
    pre=[i for i,o in enumerate(obs) if o.phase=="pre"]
    post=[i for i,o in enumerate(obs) if o.phase=="post"]
    event=next(i for i,o in enumerate(obs) if o.phase=="event")
    x,use,t=design(obs,origin)
    model=fit(stack,obs,origin,c)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        b=np.nanmedian(stack[pre],axis=0)
        spread=1.4826*np.nanmedian(np.abs(stack[pre]-b),axis=0)
    npre=np.isfinite(stack[pre]).sum(axis=0)
    npost=np.isfinite(stack[post]).sum(axis=0)
    classes=np.zeros(shape,dtype="uint8")
    quality=np.zeros(shape,dtype="uint16")
    reason=np.ones(shape,dtype="uint8")
    def flag(mask,bit):
        quality[inside&mask]|=bit
    flag(np.ones(shape,bool),1 if c.ZERO_UNRESOLVED else 0)
    flag(np.ones(shape,bool),2 if c.COMPARABILITY_UNVERIFIED else 0)
    flag(np.ones(shape,bool),4 if len(pre)<8 or len(post)<8 else 0)
    if filled is not None:
        flag(filled,128)
    any_data=np.isfinite(stack).any(axis=0)
    classes[inside&any_data]=5
    reason[inside&~any_data]=9
    reference=inside&(npre>=c.MIN_PRE)&np.isfinite(b)
    reason[inside&any_data&~reference]=2
    reason[reference&~np.isfinite(stack[event])]=3
    reference &= np.isfinite(stack[event])
    flag(inside&~reference,64)
    flag(np.ones(shape,bool),4096 if not c.LOSS_CALIBRATION_VALIDATED else 0)
    classes[reference]=5; reason[reference]=4
    support=reference&(npost>=c.MIN_POST)&np.isfinite(stack[post[-c.CONFIRM_POST:]]).all(axis=0)
    reason[support]=5
    eligible=support&model.valid
    physical=np.any((model.predicted[post]<-c.MAX_PREDICTION_EXCESS)|
                    (model.predicted[post]>1+c.MAX_PREDICTION_EXCESS),axis=0)
    unstable=eligible&((model.rmse>c.MAX_RMSE)|(spread>c.MAX_PRE_SPREAD)|physical)
    reason[unstable]=8; flag(unstable,512)
    eligible &= ~unstable
    flag(reference&~eligible&~unstable,64)
    baseline_guard=np.maximum(c.BASELINE_GUARD_FLOOR,
        c.GUARD_MULTIPLIER*spread/np.sqrt(np.maximum(npre,1)))
    step_loss=-model.beta[2]
    critical=student_t.ppf(1-c.DETECTION_ALPHA,np.maximum(model.df,1))
    step_se=contrast_se(model,np.array([0.,0.,1.,0.]))
    step_guard=critical*step_se
    threshold=np.full(shape,c.LOSS_MIN)
    event_loss=b-stack[event]
    if not c.USE_UNCERTAINTY_GUARD:
        baseline_guard=np.zeros(shape); step_guard=np.zeros(shape)
    eps=1e-7
    all_zero=inside&np.all(stack==0,axis=0)&c.ZERO_UNRESOLVED
    eligible &= ~all_zero
    flag(all_zero,8192); reason[all_zero]=10
    # Primary disturbance evidence is the modeled level drop, not a dual threshold gate.
    modeled_loss=eligible&(step_loss-step_guard>=threshold-eps)
    modeled_no_loss=eligible&(step_loss+step_guard<threshold-eps)
    event_support=event_loss-baseline_guard>=c.EVENT_EVIDENCE_MIN-eps
    event_opposes=event_loss+baseline_guard<-c.EVENT_EVIDENCE_MIN-eps
    # Do not call equal sub-threshold changes conflicting solely because of different cutoffs.
    event_against_no_loss=event_loss-baseline_guard>=np.maximum(threshold,c.EVENT_EVIDENCE_MIN)-eps
    disagreement=eligible&((modeled_loss&event_opposes)|(modeled_no_loss&event_against_no_loss))
    loss=modeled_loss&~event_opposes
    no_loss=modeled_no_loss&~event_against_no_loss
    classes[no_loss]=4; reason[no_loss]=0
    flag(disagreement,16); reason[disagreement]=7
    flag(loss&~event_support,2048)
    ambiguous=eligible&~loss&~no_loss&~disagreement
    flag(ambiguous,8); reason[ambiguous]=6
    late=post[-c.CONFIRM_POST:]
    guard=critical*model.se+baseline_guard
    if not c.USE_UNCERTAINTY_GUARD:
        guard=np.zeros_like(model.predicted)
    target=b-c.RETURN_TOLERANCE
    model_return=(model.predicted[late]-guard[late]>=target-eps).all(axis=0)
    model_below=model.predicted[post[-1]]+guard[post[-1]]<target-eps
    observed_return=(stack[late]>=target-eps).all(axis=0)
    observed_below=(stack[late]<target-eps).all(axis=0)
    conflict=loss&((model_return&~observed_return)|(model_below&~observed_below))
    flag(conflict,32); reason[conflict]=7
    recovered=loss&model_return&observed_return&~conflict
    dt=t[post[-1]]-t[post[0]]
    progress=(model.beta[1]+model.beta[3])*dt
    progress_guard=critical*contrast_se(model,np.array([0.,dt,0.,dt]))
    if not c.USE_UNCERTAINTY_GUARD:
        progress_guard=np.zeros(shape)
    partial=loss&model_below&observed_below&(progress-progress_guard>c.MIN_PROGRESS+eps)&~conflict
    failed=loss&model_below&observed_below&(progress+progress_guard<=c.MIN_PROGRESS+eps)&~conflict
    classes[recovered]=1; classes[partial]=2; classes[failed]=3
    decided=recovered|partial|failed
    reason[decided]=0
    uncertain=loss&~decided&~conflict
    flag(uncertain,8); reason[uncertain]=6
    first=model.predicted[post[0]]; final=model.predicted[post[-1]]
    slope=model.beta[1]+model.beta[3]; intercept=model.beta[0]+model.beta[2]
    candidate=np.full(shape,np.nan)
    np.divide(target-intercept,slope,out=candidate,where=eligible&(slope>0))
    flag(np.isfinite(candidate)&(candidate>t[post[-1]]),256)
    flag(np.isfinite(candidate)&(candidate<t[post[0]]),1024)
    lower=np.full(shape,np.nan); upper=lower.copy()
    for left,right in zip(post[:-1],post[1:]):
        bracket=(candidate>=t[left]-eps)&(candidate<=t[right]+eps)&(stack[left]<target)&(stack[right]>=target)
        lower[bracket]=t[left]; upper[bracket]=t[right]
    duration=np.where(recovered&np.isfinite(lower)&np.isfinite(stack[post[0]])&
                      (stack[post[0]]<target)&(candidate>=t[post[0]]-eps)&
                      (candidate<=t[post[-1]]+eps),candidate,np.nan)
    unsupported_return=recovered&~np.isfinite(duration)
    classes[unsupported_return]=5; reason[unsupported_return]=11
    flag(unsupported_return,8)
    lower[~np.isfinite(duration)]=np.nan; upper[~np.isfinite(duration)]=np.nan
    values=dict(zip(CORE_METRICS,(np.where(npre>=c.MIN_PRE,b,np.nan),first,final,b-final,progress,candidate,duration,
                lower,upper,guard[post[-1]],event_loss,step_loss,slope,model.pre_noise,step_se,
                threshold+step_guard,-event_loss)))
    for name,value in values.items():
        values[name]=np.where(inside&np.isfinite(value),value,np.nan).astype("float32")
    quality[~inside]=65535
    return Assessment(classes,quality,reason,values,model)


def assess(stack,obs,origin,c,inside=None,filled=None):
    """One fitted model; selectable operational or historical guarded decisions."""
    result=_assess_guarded(stack,obs,origin,c,inside,filled)
    from operational_recovery import classify,add_evidence
    mode=getattr(c,'EVIDENCE_MODE','guarded')
    if mode=='operational':
        return classify(result,stack,obs,origin,c,inside,filled)
    if mode!='guarded':
        raise ValueError('EVIDENCE_MODE must be operational or guarded.')
    return add_evidence(result,stack,obs,origin,c,inside,filled)
