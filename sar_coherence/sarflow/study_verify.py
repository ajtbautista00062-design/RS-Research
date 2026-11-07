"""Read back model arrays, class accounting, and every named barangay clip."""
import numpy as np
import pandas as pd
import rasterio
from .study_layout import folders


def verify(root,data,paths,c):
    totals=np.zeros(6,dtype='int64')
    rri_totals=np.zeros(7,dtype='int64')
    with rasterio.open(paths['ITSA_RRI']) as ds:
        rri=ds.read(1,masked=True).filled(np.nan)
    with rasterio.open(paths['RRI_Validity_Status']) as ds:
        rri_status=ds.read(1)
    with rasterio.open(paths['RRI_Interpretation_Classes']) as ds:
        rri_classes=ds.read(1)
    with rasterio.open(paths['Raw_ITSA_Metrics']) as ds:
        indexes={label:i+1 for i,label in enumerate(ds.descriptions)}
        b0,b1,b2,b3=[ds.read(indexes[k],masked=True).filled(np.nan) for k in ('beta0','beta1','beta2','beta3')]
    from .temporal import day_offsets
    final_day=max(day_offsets(data.observations,c.EVENT_ORIGIN)[
        np.array([o.phase=='post' for o in data.observations])])
    valid=data.inside&(rri_status==1)
    formula=np.full(rri.shape,np.nan)
    np.divide((b1+b3)*final_day,-b2,out=formula,where=valid)
    if not np.allclose(rri[valid],formula[valid],rtol=2e-5,atol=2e-6):
        raise AssertionError('RRI formula differs from saved raw pixel coefficients.')
    if np.any(np.isfinite(rri[data.inside&~valid])) or np.any(rri_classes[data.inside&~valid]!=0):
        raise AssertionError('Undefined RRI was represented as numerical recovery.')
    with rasterio.open(paths['ITSA_Recovery_Classes']) as reference:
        a=reference.read(1)
        expected=np.bincount(a[data.inside],minlength=6)
        for path in paths.values():
            with rasterio.open(path) as ds:
                if ds.crs!=reference.crs or ds.transform!=reference.transform or ds.shape!=reference.shape:
                    raise AssertionError(f'Overall grid mismatch {path.name}')
        for i,b in enumerate(data.barangays,1):
            mask=data.ownership==i
            rr,cc=np.where(mask)
            if not len(rr):
                continue
            rd=folders(root,'ITSA',b.folder).rasters
            with rasterio.open(rd/'ITSA_Recovery_Classes.tif') as clip:
                cropped=clip.read(1)
                portion=mask[rr.min():rr.max()+1,cc.min():cc.max()+1]
                if not np.array_equal(cropped[portion],a[mask]):
                    raise AssertionError(f'Barangay map differs from municipal map: {b.name}')
                totals+=np.bincount(cropped[portion],minlength=6)
            table=pd.read_csv(folders(root,'ITSA',b.folder).csv/'Recovery_Class_Area_Statistics.csv')
            if not np.array_equal(table.pixels.to_numpy(),np.bincount(a[mask],minlength=6)):
                raise AssertionError(f'Barangay CSV counts differ: {b.name}')
            rd=folders(root,'RRI',b.folder).rasters
            for name,expected_array in [('ITSA_RRI',rri),('RRI_Validity_Status',rri_status),('RRI_Interpretation_Classes',rri_classes)]:
                with rasterio.open(rd/(name+'.tif')) as ds:
                    value=ds.read(1,masked=True).filled(np.nan) if name=='ITSA_RRI' else ds.read(1)
                    if not np.allclose(value[portion],expected_array[mask],equal_nan=True):
                        raise AssertionError(f'Barangay RRI clip mismatch {b.name}: {name}')
            rri_totals+=np.bincount(rri_status[mask],minlength=7)
    if not np.array_equal(totals,expected) or int(expected.sum())!=int(data.inside.sum()):
        raise AssertionError('Barangay ownership or fixed-denominator accounting mismatch.')
    target=pd.read_csv(root/'Validation'/'Target_Assessment.csv')
    overall=target[(target['mode']=='object')&(target.zone=='Overall')].iloc[0]
    if not np.array_equal(expected,np.array([overall[f'corrected_class_{k}'] for k in range(6)])):
        raise AssertionError('Target assessment does not match reread primary raster.')
    if not np.array_equal(rri_totals,np.bincount(rri_status[data.inside],minlength=7)):
        raise AssertionError('Barangay RRI status population does not equal overall.')
    for row in target.itertuples():
        if sum(getattr(row,f'original_class_{k}') for k in range(6))!=row.population:
            raise AssertionError('Before comparison denominator changed.')
        if sum(getattr(row,f'corrected_class_{k}') for k in range(6))!=row.population:
            raise AssertionError('After comparison denominator changed.')
    return dict(overall_grid='passed',barangay_class_values='passed',barangay_csv_counts='passed',
                rri_formula='passed',barangay_rri_values_status_classes='passed',target_counts='passed',
                sum_barangays_equals_overall='passed',fixed_population='passed',study_pixels=int(expected.sum()))
