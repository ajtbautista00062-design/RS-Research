"""Conservative region adjacency regularization, never majority voting.

Only similar, same-class low-confidence neighbors can be pooled. A new ITSA fit
must reproduce their class. Indeterminate areas cannot acquire recovery labels.
Small isolated low-confidence regions are deferred to Class 5, not recolored.
"""
import numpy as np


def boundary_pairs(labels):
    pairs=set()
    for a,b in ((labels[:-1],labels[1:]),(labels[:,:-1],labels[:,1:])):
        valid=(a>0)&(b>0)&(a!=b)
        if valid.any():
            values=np.sort(np.column_stack((a[valid],b[valid])),axis=1)
            pairs.update(map(tuple,np.unique(values,axis=0).tolist()))
    return pairs


def compatible(a,b,c):
    if a['class_code'] not in (1,2,3,4) or a['class_code']!=b['class_code']:
        return False
    if max(a['evidence_score'],b['evidence_score'])>=c.CLASS_CONFIDENCE_THRESHOLD:
        return False
    if a.get('protected_high_confidence') or b.get('protected_high_confidence'):
        return False
    if a['stratum']!=b['stratum'] or a['area_m2']+b['area_m2']>c.MAX_REGION_AREA_M2:
        return False
    low=np.minimum(a['feature_min'],b['feature_min'])
    high=np.maximum(a['feature_max'],b['feature_max'])
    return bool(np.all(high-low<=np.asarray(c.PRE_FEATURE_TOLERANCES)*c.REGULARIZATION_STRENGTH))


def regularize(records,pairs,c,refit):
    """One disjoint-pair pass bounds memory/window size and avoids chaining.

    refit(ids) reaggregates original observations, fits ITSA, and returns a row.
    Returned mapping leaves the original object identifiers unchanged on disk.
    """
    output={key:dict(row) for key,row in records.items()}
    mapping={key:key for key in records}
    used=set(); accepted=set()
    for left,right in sorted(pairs):
        if left in used or right in used or left not in records or right not in records:
            continue
        a,b=records[left],records[right]
        if not compatible(a,b,c):
            continue
        merged=refit((left,right))
        if merged['class_code']!=a['class_code'] or merged['object_flags']:
            continue
        # A pooled region cannot rely on lower evidence than its stronger parent.
        if merged['evidence_score']<max(a['evidence_score'],b['evidence_score']):
            continue
        key=min(left,right)
        for member in (left,right):
            output[member]=dict(merged)
            mapping[member]=key
        used.update((left,right)); accepted.update((left,right))
    for key,row in records.items():
        if (key not in accepted and row['class_code'] in (1,2,3,4) and
                not row.get('protected_high_confidence') and
                row['area_m2']<c.MIN_REGION_AREA_M2 and
                row['evidence_score']<c.CLASS_CONFIDENCE_THRESHOLD):
            output[key]['object_decision_flags']=row.get('object_decision_flags',0)|8
            if getattr(c,'MMU_POLICY','defer_low_confidence')=='defer_low_confidence':
                output[key]['class_code']=5
                output[key]['evidence_score']=0.
                output[key]['exclusion_reason']=23
                output[key]['evidence_status']=0
                output[key]['object_flags']+=';isolated_below_MMU_deferred'
                output[key]['quality_flags']=row.get('quality_flags',0)|8
            else:
                output[key]['object_flags']+=';below_MMU_provisional'
                output[key]['quality_flags']=row.get('quality_flags',0)|16384
                output[key]['evidence_status']=2
    return output,mapping
