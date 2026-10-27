"""Window-vectorized segmented OLS; supports complete or explicitly relaxed support."""
from dataclasses import dataclass
import numpy as np
import pandas as pd

@dataclass
class Fit:
    beta: np.ndarray
    covariance: np.ndarray
    predicted: np.ndarray
    se: np.ndarray
    rmse: np.ndarray
    r_squared: np.ndarray
    df: np.ndarray
    n: np.ndarray
    valid: np.ndarray
    condition: np.ndarray
    pre_noise: np.ndarray

def design(obs, origin):
    t=np.array([(o.reference-pd.Timestamp(origin)).total_seconds()/86400 for o in obs])
    post=np.array([o.phase=="post" for o in obs],dtype=float)
    use=np.array([o.phase!="event" for o in obs])
    return np.column_stack((np.ones(len(obs)),t,post,t*post)),use,t

def fit(stack,obs,origin,c):
    x,use,_=design(obs,origin)
    shape=stack.shape[1:]; n_pixels=int(np.prod(shape))
    y=stack.reshape(len(obs),n_pixels).astype("float64")
    xx=x[use]; yy=y[use]
    good=np.isfinite(yy)
    beta=np.full((4,n_pixels),np.nan); cov=np.full((4,4,n_pixels),np.nan)
    pred=np.full_like(y,np.nan); se=np.full_like(y,np.nan)
    rmse=np.full(n_pixels,np.nan); r2=rmse.copy(); df=rmse.copy(); cond=rmse.copy()
    valid=np.zeros(n_pixels,dtype=bool)
    pre_noise=np.full(n_pixels,np.nan)
    phases=np.array([o.phase for o in obs])[use]
    enough=(good[phases=="pre"].sum(axis=0)>=c.MIN_PRE)&(good[phases=="post"].sum(axis=0)>=c.MIN_POST)
    if c.REQUIRE_COMPLETE:
        groups=[(np.ones(use.sum(),dtype=bool),np.where(enough&good.all(axis=0))[0])]
    else:
        indices=np.where(enough)[0]
        patterns,inverse=np.unique(good[:,indices].T,axis=0,return_inverse=True)
        groups=[(p,indices[inverse==i]) for i,p in enumerate(patterns)]
    for mask,cols in groups:
        if not len(cols):
            continue
        a=xx[mask]; dof=len(a)-4
        scales=np.linalg.norm(a,axis=0)
        condition=np.linalg.cond(a/np.where(scales>0,scales,1))
        cond[cols]=condition
        if dof<=0 or np.linalg.matrix_rank(a)!=4 or condition>c.MAX_CONDITION:
            continue
        response=yy[mask][:,cols]
        b=np.linalg.pinv(a)@response
        residual=response-a@b
        sse=np.sum(residual**2,axis=0)
        pre_rows=phases[mask]=="pre"
        pre_a=a[pre_rows,:2]
        pre_y=response[pre_rows]
        pre_residual=pre_y-pre_a@np.linalg.pinv(pre_a)@pre_y
        sigma=np.sqrt(np.sum(pre_residual**2,axis=0)/(pre_rows.sum()-2))
        pre_noise[cols]=sigma
        # PRE-only temporal noise protects against declaring a noisy series unchanged.
        variance=np.maximum.reduce((sse/dof,np.full(len(cols),c.ERROR_FLOOR**2),sigma**2))
        inv=np.linalg.pinv(a.T@a)
        beta[:,cols]=b
        cov[:,:,cols]=inv[:,:,None]*variance
        pred[:,cols]=x@b
        se[:,cols]=np.sqrt(np.maximum(0,np.einsum("ij,jk,ik->i",x,inv,x)[:,None]*variance))
        # Event has no instantaneous model prediction in outputs.
        pred[~use]=np.nan; se[~use]=np.nan
        rmse[cols]=np.sqrt(sse/len(a)); df[cols]=dof; valid[cols]=True
        total=np.sum((response-response.mean(axis=0))**2,axis=0)
        ratios=np.full(len(cols),np.nan)
        np.divide(sse,total,out=ratios,where=total>1e-20)
        r2[cols]=1-ratios
    def resh(a,leading=()):
        return a.reshape((*leading,*shape))
    return Fit(resh(beta,(4,)),resh(cov,(4,4)),resh(pred,(len(obs),)),
        resh(se,(len(obs),)),resh(rmse),resh(r2),resh(df),
        resh(good.sum(axis=0)),resh(valid),resh(cond),resh(pre_noise))

def contrast_se(model,vector):
    return np.sqrt(np.maximum(0,np.einsum("i,ij...,j->...",vector,model.covariance,vector)))
