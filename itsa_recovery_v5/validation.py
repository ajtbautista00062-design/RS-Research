"""Independent checks and explicitly configured sensitivity scenarios."""
from copy import copy
from dataclasses import replace
import numpy as np
from model import design

def independent_check(stack,obs,origin,c,result,limit):
    x,use,_=design(obs,origin)
    y=stack.reshape(len(obs),-1)
    cols=np.where(result.model.valid.ravel())[0][:limit]
    beta=result.model.beta.reshape(4,-1)
    predicted=result.model.predicted.reshape(len(obs),-1)
    for col in cols:
        valid=use&np.isfinite(y[:,col])
        b=np.linalg.lstsq(x[valid],y[valid,col],rcond=None)[0]
        if not np.allclose(beta[:,col],b,atol=1e-7,rtol=1e-6):
            raise RuntimeError("Independent real-pixel coefficient check failed.")
        if not np.allclose(predicted[valid,col],x[valid]@b,atol=1e-7,rtol=1e-6):
            raise RuntimeError("Independent real-pixel prediction check failed.")
    return len(cols)

def variants(c,obs):
    scenarios=[]
    for name,field,value in (("loss_001","LOSS_MIN",.01),("loss_002","LOSS_MIN",.02),("loss_003","LOSS_MIN",.03),
                             ("loss_005","LOSS_MIN",.05),("loss_010","LOSS_MIN",.10),
                             ("return_003","RETURN_TOLERANCE",.03),("return_007","RETURN_TOLERANCE",.07),
                             ("incomplete_support","REQUIRE_COMPLETE",False)):
        changed=copy(c); setattr(changed,field,value)
        scenarios.append((name,changed,obs))
    strict=copy(c); strict.MIN_PRE=4; strict.MIN_POST=4
    scenarios.append(("minimum_4_each_phase",strict,obs))
    mid=[replace(o,reference=o.start+(o.end-o.start)/2) for o in obs]
    scenarios.append(("midpoint_indexing",copy(c),mid))
    return scenarios

def validate_settings(c):
    if getattr(c,'EVIDENCE_MODE','guarded') not in ('operational','guarded'):
        raise ValueError('Invalid classification evidence mode.')
    if getattr(c,'EVIDENCE_MODE','guarded')=='operational':
        for key in ('EARLY_LOSS_MIN','EFFECT_MARGIN','OBS_CONFIRM_TOLERANCE',
                    'MAX_OBSERVATION_DISAGREEMENT','MODEL_RMSE_FAILURE'):
            if not np.isfinite(getattr(c,key)) or getattr(c,key)<=0:
                raise ValueError(f'Invalid operational setting: {key}')
        if not 0<=c.LOSS_NOISE_MULTIPLIER<=5 or c.EFFECT_MARGIN>min(c.LOSS_MIN,c.MIN_PROGRESS):
            raise ValueError('Noise multiplier or operational margin invalid.')
        if not isinstance(c.REQUIRE_EVENT,(bool,np.bool_)):
            raise ValueError('REQUIRE_EVENT must be boolean.')
    if c.MIN_PRE<3 or c.MIN_POST<3 or c.CONFIRM_POST<2:
        raise ValueError("Minimum 3 PRE/POST and 2 late confirmations required.")
    if not (0<c.RETURN_TOLERANCE<1 and 0<c.LOSS_MIN<=1 and 0<c.DETECTION_ALPHA<.5) or c.MIN_PROGRESS<=0:
        raise ValueError("Invalid loss/return/progress thresholds.")
    for name in ("PRE_SPREAD_FLOOR","MAX_RMSE","MAX_PRE_SPREAD","ERROR_FLOOR",
                 "BASELINE_GUARD_FLOOR","GUARD_MULTIPLIER","MAX_CONDITION"):
        if not np.isfinite(getattr(c,name)) or getattr(c,name)<=0:
            raise ValueError(f"Invalid setting: {name}")
