"""Exact area summaries; class denominators include NoData within owned study cells."""
import numpy as np
import pandas as pd
from .recovery import CLASSES, QUALITY, REASONS
from .rri import INTERPRETATION, STATUS
from .temporal import day_offsets


def itsa_tables(member,classes,reasons,quality,metrics,stack,obs,result,folder,name,area,c,fit_series=None):
    total=int(member.sum())
    rows=[dict(code=k,label=CLASSES[k],pixels=int(((classes==k)&member).sum()),
               area_m2=int(((classes==k)&member).sum())*area,
               percent_study=100*int(((classes==k)&member).sum())/total if total else np.nan,
               study_pixels=total) for k in range(6)]
    pd.DataFrame(rows).to_csv(folder/'Recovery_Class_Area_Statistics.csv',index=False)
    pd.DataFrame([dict(reason=k,label=label,pixels=int(((reasons==k)&member).sum()))
                  for k,label in REASONS.items()]).to_csv(folder/'Exclusion_Reasons.csv',index=False)
    pd.DataFrame([dict(bit=k,label=label,pixels=int((((quality&k)>0)&member).sum()))
                  for k,label in QUALITY.items()]).to_csv(folder/'Quality_Flags.csv',index=False)
    summaries=[]
    for key,a in metrics.items():
        x=a[member&np.isfinite(a)]
        summaries.append(dict(indicator=key,valid_pixels=len(x),mean=float(x.mean()) if len(x) else np.nan,
                              median=float(np.median(x)) if len(x) else np.nan,
                              minimum=float(x.min()) if len(x) else np.nan,maximum=float(x.max()) if len(x) else np.nan))
    pd.DataFrame(summaries).to_csv(folder/'ITSA_Indicator_Statistics.csv',index=False)
    days=day_offsets(obs,c.EVENT_ORIGIN)
    post=np.array([o.phase=='post' for o in obs])
    use=np.array([o.phase!='event' for o in obs])
    y=np.array([np.nanmedian(a[member]) if np.isfinite(a[member]).any() else np.nan for a in stack])
    if fit_series is not None:
        y=np.asarray(fit_series,dtype=float)
    x=np.column_stack((np.ones(len(obs)),days,post,days*post))
    residual=y-x@result.fit.beta[:,0]
    finite=use&np.isfinite(y)
    df=int(finite.sum())-4
    covariance=np.full((4,4),np.nan)
    if df>0 and np.linalg.matrix_rank(x[finite])==4:
        covariance=np.linalg.inv(x[finite].T@x[finite])*np.sum(residual[finite]**2)/df
    pd.DataFrame(dict(term=['beta0','beta1','beta2','beta3'],estimate=result.fit.beta[:,0],
                      standard_error=np.sqrt(np.maximum(0,np.diag(covariance))),residual_df=df,
                      scope='separate descriptive aggregate median fit; not map classification')).to_csv(
                          folder/'ITSA_Regression_Results.csv',index=False)
    pd.DataFrame(dict(observation=[o.name for o in obs],included=finite,residual=residual)).to_csv(folder/'ITSA_Residuals.csv',index=False)


def rri_tables(member,values,classes,status,raw_classes,folder,area,coverage,c):
    total=int(member.sum())
    valid=member&(status==1)&np.isfinite(values)
    a=values[valid]
    coverage_eligible=bool((coverage[:,member].mean(1)>=c.MIN_PER_DATE_VALID_COVERAGE).all()) if total else False
    row=dict(study_pixels=total,valid_rri_pixels=len(a),undefined_rri_pixels=total-len(a),
             coverage_eligible=coverage_eligible,mean=float(a.mean()) if len(a) and coverage_eligible else np.nan,
             median=float(np.median(a)) if len(a) and coverage_eligible else np.nan,
             minimum=float(a.min()) if len(a) and coverage_eligible else np.nan,
             maximum=float(a.max()) if len(a) and coverage_eligible else np.nan,
             formula='(modeled final POST - (beta0+beta2))/(-beta2)',
             eligibility='resolved raw pixel ITSA classes 1-3; unchanged modeled denominator gates')
    pd.DataFrame([row]).to_csv(folder/'RRI_Statistics.csv',index=False)
    pd.DataFrame([dict(code=k,label=STATUS[k],pixels=int(((status==k)&member).sum()),study_pixels=total)
                  for k in STATUS]).to_csv(folder/'RRI_Validity_Status.csv',index=False)
    pd.DataFrame([dict(code=k,label=INTERPRETATION[k],pixels=int(((classes==k)&valid).sum()),
                       area_m2=int(((classes==k)&valid).sum())*area) for k in range(1,6)]).to_csv(folder/'RRI_Class_Distribution.csv',index=False)
    pd.DataFrame([dict(raw_itsa_class=i,rri_class=j,pixels=int((valid&(raw_classes==i)&(classes==j)).sum()),
                       note='Independent definitions; not an accuracy confusion matrix')
                  for i in range(6) for j in range(1,6)]).to_csv(folder/'RRI_ITSA_Agreement.csv',index=False)
