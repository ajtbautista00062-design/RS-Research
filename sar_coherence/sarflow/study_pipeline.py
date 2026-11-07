"""Seven-pair controlled comparison and requested RRI/ITSA product workflow."""
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
import hashlib
import json
import logging
import numpy as np
import pandas as pd
from .study_inputs import load, effective_settings
from .study_layout import folders
from .study_comparison import Comparison
from .study_rasters import RasterProducts
from .study_objects import evaluate, project
from .spatial_objects import cell_sizes
from .spatial import windows
from .recovery import assess
from .baseline_recovery import assess as original_assess
from .recovery_io import METRICS
from .rri import calculate
from .temporal import day_offsets
from .study_exports import export_all
from .study_verify import verify


def configuration(c, **changes):
    return SimpleNamespace(**({k:v for k,v in vars(c).items() if not k.startswith('_')} | changes))


def write_window(products, win, stack, mask, ownership, raw, primary, original_raw, original, ids, rri):
    arrays = dict(ITSA_Recovery_Classes=primary.classes,Original_ITSA_Recovery_Classes=original.classes,
                  Raw_Pixel_ITSA_Classes=raw.classes,Original_Raw_Pixel_ITSA_Classes=original_raw.classes,
                  ITSA_Exclusion_Reasons=primary.reasons,ITSA_Quality_Flags=primary.quality,
                  Raw_ITSA_Reasons=raw.reasons,Raw_ITSA_Quality_Flags=raw.quality,
                  Spatial_Object_ID=ids,Barangay_Owner=ownership,Study_Mask=mask,
                  Aligned_Coherence_Observations=stack,
                  Valid_Data_Coverage=np.where(mask[None],np.isfinite(stack),255),
                  ITSA_Metrics=np.stack([primary.metrics[k] for k in METRICS]),
                  Raw_ITSA_Metrics=np.stack([raw.metrics[k] for k in METRICS]),
                  ITSA_Gate_Outcomes=np.stack([np.where(mask,g,255) for g in primary.gates.values()]))
    for name,key in [('ITSA_Baseline','baseline'),('ITSA_Disturbance_Magnitude','modeled_level_loss'),
                     ('ITSA_Final_Deficit','baseline_deficit'),('ITSA_Recovery_Progress','modeled_progress'),
                     ('ITSA_Recovery_Fraction','recovery_fraction'),('ITSA_Crossing_Days','crossing_days')]:
        arrays[name]=primary.metrics[key]
    for name,value in zip(('ITSA_RRI','RRI_Interpretation_Classes','RRI_Validity_Status',
                           'RRI_Disturbance_Magnitude','RRI_Recovery_Magnitude'),rri):
        arrays[name]=value
    for name,a in arrays.items():
        products.write(name,a,win)


def sensitivity(comparison, series, ids, keys, stack, mask, owner, obs, c, area):
    scenarios = [('noise_2',dict(OBSERVABILITY_NOISE_MULTIPLIER=2.)),
                 ('noise_4',dict(OBSERVABILITY_NOISE_MULTIPLIER=4.)),
                 ('legacy_absolute',dict(BASELINE_ELIGIBILITY_MODE='legacy_absolute')),
                 ('loss_002',dict(EVENT_LOSS_MINIMUM=.02)),('loss_005',dict(EVENT_LOSS_MINIMUM=.05)),
                 ('return_003',dict(RECOVERY_TOLERANCE=.03)),('return_007',dict(RECOVERY_TOLERANCE=.07))]
    for name,changes in scenarios:
        settings=configuration(c,**changes)
        if keys.size:
            result=project(assess(series,obs,settings),ids,keys,stack,mask,area,settings)
        else:
            result=assess(stack,obs,settings)
        comparison.add_sensitivity(name,result.classes,mask,owner)


def process(data, root, c):
    obs=data.observations
    destinations={m:folders(root,m) for m in ('ITSA','RRI')}
    comparison=Comparison()
    area=abs(data.profile['transform'].a*data.profile['transform'].e)*data.profile['crs'].linear_units_factor[1]**2
    cell,_=cell_sizes(area,c)
    chunk=cell*max(1,c.WINDOW_SIZE//cell)
    pre=[o.phase=='pre' for o in obs]
    final_day=max(day_offsets(obs,c.EVENT_ORIGIN)[np.array([o.phase=='post' for o in obs])])
    with ExitStack() as context:
        products=RasterProducts(context,destinations,data.profile,c,obs)
        for win in windows(data.inside.shape[1],data.inside.shape[0],chunk):
            sl=np.s_[int(win.row_off):int(win.row_off+win.height),int(win.col_off):int(win.col_off+win.width)]
            stack,mask,owner=data.stack[:,sl[0],sl[1]],data.inside[sl],data.ownership[sl]
            raw=assess(stack,obs,c)
            original_raw=original_assess(stack,obs,c)
            primary,original,ids,keys,counts,series=evaluate(stack,mask,obs,c,area,
                int(data.source_window.row_off+win.row_off),int(data.source_window.col_off+win.col_off),
                int(data.profile['source_width']))
            for result in (raw,original_raw):
                result.reasons[~mask]=11
                result.quality[~mask]=65535
            products.ensure_gates(context,destinations['ITSA'].rasters,data.profile,list(primary.gates))
            comparison.add('object',original,primary,mask,owner)
            comparison.add('pixel',original_raw,raw,mask,owner)
            if mask.any():
                sensitivity(comparison,series,ids,keys,stack,mask,owner,obs,c,area)
            b=raw.metrics['beta0'];e=b+raw.metrics['beta2']
            final=e+(raw.metrics['beta1']+raw.metrics['beta3'])*final_day
            rri=calculate(b,e,final,np.isin(raw.classes,(1,2,3)),raw.gates['temporal_support'],
                          raw.gates['stable_model'],mask,raw.metrics['adaptive_loss_threshold'],
                          no_loss=raw.classes==4)
            rri[2][mask&(raw.classes==0)]=0
            old_rri=calculate(b,e,final,np.isin(original_raw.classes,(1,2,3)),original_raw.gates['temporal_support'],
                             original_raw.gates['stable_model'],mask,original_raw.metrics['adaptive_loss_threshold'],
                             no_loss=original_raw.classes==4)
            old_rri[2][mask&(original_raw.classes==0)]=0
            comparison.add_rri(old_rri[2],rri[2],mask,owner)
            write_window(products,win,stack,mask,owner,raw,primary,original_raw,original,ids,rri)
            logging.info('Processed row=%s col=%s study_cells=%s objects=%s',win.row_off,win.col_off,int(mask.sum()),len(keys))
        paths={name:Path(ds.name) for name,ds in products.outputs.items()}
    names={0:'Overall',**{i:b.name for i,b in enumerate(data.barangays,1)}}
    summaries=comparison.export(root/'Validation',names)
    return destinations,paths,summaries


def run(c, run_name=None):
    if c.CLASSIFICATION_MODE!='object_itsa':
        raise ValueError('This controlled primary implementation uses fixed PRE-only object ITSA; raw pixel products are separate.')
    if c.ENABLE_NEIGHBOR_CLASS_ASSIGNMENT or c.ENABLE_SPATIAL_REGULARIZATION or c.ENABLE_SPATIAL_REGION_GROWTH:
        raise ValueError('Spatial propagation and merging are disabled for this controlled experiment.')
    # Short default name: Windows MAX_PATH (260) limits the deepest per-barangay
    # export: <run>/B. ITSA/Per Barangay/<name_code>/ITSA Graph/<~40 char file>.
    root=c.VERSION_ROOT/(run_name or datetime.now(timezone.utc).strftime('%m%d_%H%M%S'))
    if not root.resolve().is_relative_to(c.RESULTS_ROOT.resolve()):
        raise ValueError('Run directory escapes the designated results workspace.')
    from .analysis_boundaries import load_barangays
    longest=max([len(b.folder) for b in load_barangays(c.BARANGAY_BOUNDARY,'EPSG:4326',c)[0]] or [0])
    deepest=len(str(root.resolve()))+73+longest
    if deepest>259:
        raise ValueError(f'Output paths would reach {deepest} characters (Windows limit 259); use a shorter '
                         f'--run-name or a shorter Output_Folder in config.py.')
    root.mkdir(parents=True,exist_ok=False)
    state=root/'run_state.json'
    state.write_text(json.dumps(dict(status='running')),encoding='utf-8')
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(message)s',force=True,
                        handlers=[logging.FileHandler(root/'processing.log',encoding='utf-8'),logging.StreamHandler()])
    try:
        frozen_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')}
        data=load(c)
        c=effective_settings(c,data.observations)
        counts=dict(pre=c.EXPECTED_PRE_COUNT,post=c.EXPECTED_POST_COUNT)
        logging.info('Inputs: %s PRE, 1 event, %s POST; ITSA_MIN_PRE=%s ITSA_MIN_POST=%s',
                     counts['pre'],counts['post'],c.ITSA_MIN_PRE,c.ITSA_MIN_POST)
        for o in data.observations:
            logging.info('  %-6s %s  %s -> %s  %s',o.phase,o.name,o.start.date(),o.end.date(),o.path)
        # This key is carried separately from the output crop dimensions.
        data.profile['source_width']=data.inventory[0]['width']
        destinations,paths,summaries=process(data,root,c)
        export_all(data,root,paths,c)
        verification=verify(root,data,paths,c)
        if frozen_hashes!={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')}:
            raise RuntimeError('Source changed during processing; rerun from frozen reviewed code.')
        manifest=dict(settings={k:v for k,v in vars(c).items() if k.isupper() or k=='Study_Area'},
            dates=[vars(o) for o in data.observations],study_pixels=int(data.inside.sum()),
            boundary=data.boundary_audit,source_window=str(data.source_window),
            reference_grid='first chronological PRE; same CRS pixel size origin; cropped to study bounding window',
            source_inventory=data.inventory,targets=[r for r in summaries if r['zone']=='Overall'],
            source_sha256=frozen_hashes,
            scientific_status=f'Exploratory {len(data.observations)}-pair operational amendment (signal observability gate, multi-scale zonal product); no field-accuracy claim',
            verification=verification)
        if getattr(c,'ENABLE_ZONAL_RECOVERY',False):
            from .zonal_recovery import run_zonal
            zonal_root,zonal=run_zonal(root,c)
            manifest['zonal_recovery']=dict(folder=str(zonal_root),overall=zonal)
            logging.info('Zonal: indeterminate reduction %.2f%%, classes 1-3 increase %.1f%%',
                         zonal['indeterminate_reduction_percent'],zonal['recovery_class_increase_percent'])
        (root/'run_manifest.json').write_text(json.dumps(manifest,indent=2,default=str),encoding='utf-8')
        state.write_text(json.dumps(dict(status='complete',verification=verification)),encoding='utf-8')
        return root
    except BaseException as exc:
        state.write_text(json.dumps(dict(status='interrupted' if isinstance(exc,KeyboardInterrupt) else 'failed',error=str(exc))),encoding='utf-8')
        raise
    finally:
        logging.shutdown()
