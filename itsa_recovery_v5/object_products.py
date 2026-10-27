"""Class/indicator projection retains a direct region-ID to model-table relationship."""
import json
import numpy as np
import pandas as pd
from recovery import QUALITY,EXCLUSION,EVIDENCE_METRICS,metric_unit,EVIDENCE_TAGS
from spatial_io import grid_windows,read_values,strata_for,writer,project_field


def export_products(resources,c,folder,labels,stack,reason,pixels,barrier,builtup,
                    records,regular,mapping,tags):
    fields={'Object_ITSA_Recovery_Classes':('class_code','uint8',0),
        'ITSA_Recovery_Classes':('class_code','uint8',0),
        'ITSA_Recovery_Confidence':('evidence_score','float32',c.NODATA),
        'object_quality_flags':('quality_flags','uint16',65535),
        'object_decision_flags':('object_decision_flags','uint8',255),
        'object_exclusion_reason':('exclusion_reason','uint8',255),
        'Object_vs_Pixel_Change':(None,'uint8',0)}
    # Names match the established maps/continuous metric interfaces.
    for name in ('pre_baseline','modeled_first_post','modeled_final_post','fitted_level_loss',
                 'event_loss','post_slope','baseline_deficit_final','modeled_recovery_progress',
                 'recovery_duration_days','final_guard_halfwidth'):
        fields[name]=(name,'float32',c.NODATA)
    for name in EVIDENCE_METRICS:
        fields[name]=(name,'float32',c.NODATA)
    fields['ITSA_Evidence_Status']=('evidence_status','uint8',255)
    if regular is not None:
        fields.update(Spatially_Regularized_ITSA_Recovery_Classes=('class_code','uint8',0),
                      Regularized_Object_ID=('regular_id','uint32',0),
                      Regularization_Change=(None,'uint8',0))
    writers={name:writer(resources,folder,name,stack,dtype,nodata,tags)
             for name,(_,dtype,nodata) in fields.items()}
    for name,dst in writers.items():
        unit=('days since event origin' if name=='recovery_duration_days' else
              'coherence/day' if name=='post_slope' else 'coherence' if fields[name][1]=='float32'
              else 'code')
        if name=='ITSA_Recovery_Confidence':
            unit='operational evidence score [0,1]; not probability'
        elif name in EVIDENCE_METRICS:
            unit=metric_unit(name)
        if name in EVIDENCE_TAGS:
            dst.update_tags(code_definitions=json.dumps(EVIDENCE_TAGS[name]))
        if name=='ITSA_Evidence_Status':
            dst.update_tags(code_definitions=json.dumps(EVIDENCE_TAGS['evidence_status']))
        dst.update_tags(1,units=unit)
        if name.endswith('Change'):
            dst.update_tags(code_0='outside/unavailable',code_1='unchanged',code_2='class changed')
        if name=='object_decision_flags':
            dst.update_tags(bit_1='below mapping unit',bit_2='heterogeneous temporal evidence',
                            bit_4='conflicting strong pixel evidence',bit_8='MMU quality flag or optional deferral',
                            bit_16='no PRE-supported object')
        if name=='object_quality_flags':
            dst.update_tags(bit_definitions=json.dumps(QUALITY))
        if name=='object_exclusion_reason':
            dst.update_tags(reason_definitions=json.dumps(EXCLUSION))
    maprows={key:dict(regular_id=value) for key,value in mapping.items()} if mapping else {}
    agreement=np.zeros((6,6),dtype='int64')
    primary_counts=np.zeros(6,dtype='int64')
    for window in grid_windows(labels,c.WINDOW_SIZE):
        ids=labels.read(1,window=window)
        strata=strata_for(window,ids.shape,barrier,builtup)
        domain=(reason.read(1,window=window)!=1)&(strata>0)
        if not domain.any():
            continue  # Unwritten GeoTIFF blocks retain declared NoData.
        available=domain&np.isfinite(read_values(stack,window)).any(axis=0)
        raw=project_field(ids,records,'class_code',0).astype('uint8')
        raw[(ids==0)&available]=5
        active=regular if regular is not None else records
        final=project_field(ids,active,'class_code',0).astype('uint8')
        final[(ids==0)&available]=5
        primary_counts+=np.bincount(final[domain],minlength=6)
        pixel=pixels.read(1,window=window)
        agreement+=np.bincount((raw[available]*6+pixel[available]),minlength=36).reshape(6,6)
        for name,(field,dtype,nodata) in fields.items():
            if name=='Object_ITSA_Recovery_Classes':
                values=raw
            elif name in ('ITSA_Recovery_Classes','Spatially_Regularized_ITSA_Recovery_Classes'):
                values=final
            elif name=='Regularization_Change':
                values=np.where(available,np.where(raw==final,1,2),0)
            elif name=='Object_vs_Pixel_Change':
                values=np.where(available,np.where(raw==pixel,1,2),0)
            elif name=='Regularized_Object_ID':
                values=project_field(ids,maprows,'regular_id',0)
            else:
                values=project_field(ids,active,field,nodata)
                if name=='ITSA_Recovery_Confidence':
                    values[(ids==0)&available]=0
                elif name=='ITSA_Evidence_Status':
                    values[(ids==0)&available]=0
                elif name=='object_quality_flags':
                    values[(ids==0)&available]=64
                elif name=='object_exclusion_reason':
                    values[(ids==0)&available]=2
                elif name=='object_decision_flags':
                    values[(ids==0)&available]=16
            values=np.where(np.isfinite(values),values,nodata).astype(dtype)
            writers[name].write(values,1,window=window)
    for dst in writers.values():
        dst.close()
    pd.DataFrame([dict(object_class=i,pixel_class=j,pixels=int(agreement[i,j]))
        for i in range(6) for j in range(6)]).to_csv(folder.parent/'Reports'/'pixel_object_agreement.csv',index=False)
    denominator=int(primary_counts.sum())
    pd.DataFrame([dict(code=i,pixels=int(n),percent_of_analysis_domain=100*n/denominator if denominator else np.nan,
        denominator_pixels=denominator,includes_structurally_unavailable=True)
        for i,n in enumerate(primary_counts)]).to_csv(folder.parent/'Reports'/'primary_map_class_counts.csv',index=False)


def export_tables(root,obs,records,regular,geometry,c):
    reports=root/'Reports'
    def table(rows):
        return pd.DataFrame([dict(row,feature_min=json.dumps(row['feature_min']),
            feature_max=json.dumps(row['feature_max'])) for row in rows])
    table(records.values()).to_csv(reports/'region_statistics.csv',index=False)
    from object_diagnostics import export as export_diagnostics
    export_diagnostics(records,reports)
    if regular is not None:
        unique={row['region_id']:row for row in regular.values()}
        table(unique.values()).to_csv(reports/'regularized_region_statistics.csv',index=False)
    pd.DataFrame([dict(observation=j+1,name=getattr(o,'name',str(j+1)),phase=o.phase,
        pair_start=o.start.isoformat(),pair_end=o.end.isoformat(),reference=o.reference.isoformat())
        for j,o in enumerate(obs)]).to_csv(reports/'region_observation_dates.csv',index=False)
    information=dict(mode=c.CLASSIFICATION_MODE,objects=len(records),geometry=geometry,
        builtup_mask_verified=bool(c.VERIFIED_BUILTUP_MASK),
        source='unaltered aligned unfilled temporal observations',
        confidence='operational guard-width score; not calibrated probability',
        segmentation_outcomes_used=False,regularization_indeterminate_promotion=False,
        pixel_reference_published=c.EXPORT_PIXEL_LEVEL_REFERENCE,
        parameter_status='provisional physical-scale and similarity defaults; independent validation required')
    (reports/'spatial_method.json').write_text(json.dumps(information,indent=2),encoding='utf-8')
