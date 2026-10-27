"""Future-run object diagnostics, without scanning full rasters a second time."""
import numpy as np
import pandas as pd
from recovery import LABELS,EXCLUSION


def export(records,reports):
    table=pd.DataFrame(records.values())
    rows=[]
    for code,label in LABELS.items():
        subset=table[table.class_code==code]
        rows.append(dict(code=code,label=label,objects=len(subset),pixels=int(subset.pixels.sum()),
            pixel_denominator='PRE-segmented members only; see primary_map_class_counts for full analysis domain',
            object_percent=100*len(subset)/len(table),
            pixel_percent=100*subset.pixels.sum()/table.pixels.sum(),
            guard_supported_conditional=int((subset.evidence_status==1).sum()),
            provisional=int((subset.evidence_status==2).sum())))
    pd.DataFrame(rows).to_csv(reports/'object_class_evidence_counts.csv',index=False)
    pd.DataFrame([dict(code=code,reason=label,
        objects=int((table.exclusion_reason==code).sum()),
        pixels=int(table.loc[table.exclusion_reason==code,'pixels'].sum()))
        for code,label in EXCLUSION.items()]).to_csv(reports/'object_rejection_reasons.csv',index=False)
    distribution=[]
    for name in ('fitted_level_loss','post_slope','baseline_deficit_final','modeled_recovery_progress'):
        values=table[name].replace([np.inf,-np.inf],np.nan).dropna()
        if len(values):
            distribution.append(dict(indicator=name,objects=len(values),**{
                key:float(v) for key,v in zip(('p05','q25','median','q75','p95'),values.quantile([.05,.25,.5,.75,.95]))}))
    pd.DataFrame(distribution).to_csv(reports/'object_indicator_distributions.csv',index=False)
    candidates=[name for name in table if name.startswith('sensitivity_')]
    sensitivity=[]
    for name in candidates:
        changed=table[name]!=table.class_code
        for code in LABELS:
            sensitivity.append(dict(candidate=name,code=code,
                objects=int((table[name]==code).sum()),
                pixels=int(table.loc[table[name]==code,'pixels'].sum()),
                changed_objects=int(changed.sum()),changed_pixels=int(table.loc[changed,'pixels'].sum()),
                scope='fixed PRE objects; actual observations; primary pixel conflict reference held fixed'))
    if sensitivity:
        pd.DataFrame(sensitivity).to_csv(reports/'object_threshold_sensitivity.csv',index=False)
    guarded=table.guarded_reference_class
    pd.DataFrame([dict(operational=left,guarded=right,
        objects=int(((table.class_code==left)&(guarded==right)).sum()))
        for left in LABELS for right in LABELS]).to_csv(reports/'object_guarded_comparison.csv',index=False)
    flag_rows=[]
    for flag in ('below_mapping_unit','heterogeneous_temporal_evidence','conflicting_high_confidence_pixel_reference'):
        flagged=table.object_flags.fillna('').str.contains(flag,regex=False)
        flag_rows.append(dict(flag=flag,objects=int(flagged.sum()),
            indeterminate=int((flagged&(table.class_code==5)).sum())))
    pd.DataFrame(flag_rows).to_csv(reports/'object_spatial_quality_counts.csv',index=False)
