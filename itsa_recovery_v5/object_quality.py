"""Future-run reread checks linking exported maps to region tables; never on import."""
from contextlib import ExitStack
from collections import Counter
import numpy as np
import rasterio
from spatial_io import project_field


def verify(folder,records,regular,c):
    active=regular if regular is not None else records
    with ExitStack() as resources:
        readers={name:resources.enter_context(rasterio.open(folder/f'{name}.tif','r'))
            for name in ('Spatial_Object_ID','Object_ITSA_Recovery_Classes','ITSA_Recovery_Classes',
                         'ITSA_Recovery_Confidence','modeled_recovery_progress')}
        ids=readers['Spatial_Object_ID'];counts=Counter()
        grid=(ids.crs,ids.transform,ids.shape)
        for src in readers.values():
            if (src.crs,src.transform,src.shape)!=grid:
                raise RuntimeError('Object output grid mismatch.')
        for _,window in ids.block_windows(1):
            labels=ids.read(1,window=window);valid=labels>0
            keys,n=np.unique(labels[valid],return_counts=True)
            counts.update({int(key):int(value) for key,value in zip(keys,n)})
            if any(int(key) not in records for key in keys):
                raise RuntimeError('Unknown object identifier in map.')
            for name,rows in (('Object_ITSA_Recovery_Classes',records),('ITSA_Recovery_Classes',active)):
                classes=readers[name].read(1,window=window)
                expected=project_field(labels,rows,'class_code',0)
                if not np.isin(classes,range(6)).all() or not np.array_equal(classes[valid],expected[valid]):
                    raise RuntimeError(f'Class map/table mismatch: {name}')
            for name,field in (('ITSA_Recovery_Confidence','evidence_score'),
                               ('modeled_recovery_progress','modeled_recovery_progress')):
                values=readers[name].read(1,window=window)
                expected=project_field(labels,active,field,c.NODATA)
                expected=np.where(np.isfinite(expected),expected,c.NODATA)
                if not np.allclose(values[valid],expected[valid],rtol=1e-6,atol=1e-6):
                    raise RuntimeError(f'Indicator map/table mismatch: {name}')
        for key,row in records.items():
            if counts[key]!=row['pixels']:
                raise RuntimeError('Object area/count mismatch.')
            final=active[key]
            if row['class_code']==5 and final['class_code']!=5:
                raise RuntimeError('Regularization promoted an indeterminate object.')
            if (row['evidence_score']>=c.CLASS_CONFIDENCE_THRESHOLD or
                    row.get('protected_high_confidence')) and row['class_code']!=final['class_code']:
                raise RuntimeError('Regularization overwrote protected evidence.')
    return dict(region_table_projection='passed',region_pixel_counts='passed',
                protected_evidence='passed',objects_verified=len(records))
