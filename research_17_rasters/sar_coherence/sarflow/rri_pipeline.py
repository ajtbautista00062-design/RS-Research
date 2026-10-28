"""Saved raw pixel ITSA -> independent RRI rasters and exact zonal statistics."""
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window, transform as window_transform
from rasterio.features import geometry_window
from rasterio.errors import WindowError
from .rri import calculate, INTERPRETATION, STATUS
from .spatial import windows
from .temporal import local_path
from .analysis_boundaries import load_barangays, owners
from .barangay_outputs import DiskValues

NAMES=('ITSA_RRI','RRI_Interpretation_Classes','RRI_Validity_Status',
       'RRI_Disturbance_Magnitude','RRI_Recovery_Magnitude')


class Summary:
    def __init__(self,scratch,name):
        self.name=name; self.values=DiskValues(scratch)
        self.status=np.zeros(7,dtype='int64');self.classes=np.zeros(6,dtype='int64')
        self.squares={};self.extrema={};self.at_least_one=0;self.disagreements=np.zeros((6,6),dtype='int64')

    def add(self,products,recovery_classes,member):
        rri,category,status=products[:3]
        self.status+=np.bincount(status[member],minlength=7)
        valid=member&(status==1)
        self.classes+=np.bincount(category[valid],minlength=6)
        self.at_least_one+=int((valid&(rri>=1)).sum())
        self.disagreements+=np.bincount((recovery_classes[valid]*6+category[valid]),minlength=36).reshape(6,6)
        for key,mask in [('rri',valid)]+[(f'class_{k}',valid&(recovery_classes==k)) for k in range(6)]:
            a=rri[mask]
            self.values.add(key,a)
            if a.size:
                self.squares[key]=self.squares.get(key,0.)+float(np.square(a.astype(float)).sum())
                lo,hi=self.extrema.get(key,(np.inf,-np.inf))
                self.extrema[key]=(min(lo,float(a.min())),max(hi,float(a.max())))

    def row(self,area,key='rri'):
        stats=self.values.stats(key);n=stats['count'];lo,hi=self.extrema.get(key,(np.nan,np.nan))
        variance=max(0.,self.squares.get(key,0.)/n-stats['mean']**2) if n else np.nan
        return dict(zone=self.name,**stats,minimum=lo,maximum=hi,std_population=np.sqrt(variance),
                    valid_area_m2=n*area,valid_area_ha=n*area/10000)

    def export(self,folder,area,barangay=False):
        row=self.row(area)
        row.update(percent_rri_ge_1=100*self.at_least_one/row['count'] if row['count'] else np.nan)
        pd.DataFrame([row]).to_csv(folder/('RRI_Statistics.csv' if barangay else 'RRI_Study_Area_Statistics.csv'),index=False)
        pd.DataFrame([dict(category=k,interpretation=INTERPRETATION[k],pixels=int(n),area_m2=n*area,
                           area_ha=n*area/10000,percent_valid=100*n/max(1,self.status[1]))
                      for k,n in enumerate(self.classes) if k]).to_csv(folder/'RRI_Class_Distribution.csv',index=False)
        pd.DataFrame([dict(recovery_class=k,**self.row(area,f'class_{k}')) for k in range(6)]).to_csv(
            folder/'RRI_By_Recovery_Class.csv',index=False)
        total=int(self.status.sum());inside=total-int(self.status[0])
        pd.DataFrame([dict(status=k,meaning=STATUS[k],pixels=int(n),study_pixels=inside,
                           percent_study=100*n/inside if k and inside else np.nan,
                           total_owned_pixels=total) for k,n in enumerate(self.status)]).to_csv(
            folder/'RRI_Processing_Diagnostics.csv',index=False)
        pd.DataFrame([dict(original_recovery_class=i,rri_interpretation=j,pixels=int(self.disagreements[i,j]),
                           disagreement=(i==1 and j!=4) or (i==2 and j!=3) or (i==3 and j not in (1,2)),
                           note='Different reference and endpoint definitions; disagreement is not an accuracy error')
                      for i in (1,2,3) for j in range(1,6)]).to_csv(folder/'RRI_Classification_Agreement.csv',index=False)
        return row


def create_rasters(context,folder,ref,settings):
    result={}
    for name in NAMES:
        categorical=name in NAMES[1:3]
        profile=ref.profile.copy()
        profile.update(count=1,dtype='uint8' if categorical else 'float32',nodata=0 if categorical else settings.NODATA,
                       tiled=True,blockxsize=256,blockysize=256,compress='deflate',predictor=2 if categorical else 3,BIGTIFF='IF_SAFER')
        result[name]=context.enter_context(rasterio.open(folder/(name+'.tif'),'w',**profile))
        result[name].set_band_description(1,name)
        result[name].update_tags(method='ITSA final-observation RRI',baseline='beta0 at intervention',
                                event='beta0+beta2 modeled POST intercept at intervention; not event-pair instantaneous coherence',
                                zero_tolerance=settings.ZERO_TOLERANCE,one_tolerance=settings.ONE_TOLERANCE,
                                units='code' if categorical else 'dimensionless' if name==NAMES[0] else 'coherence',
                                codebook=json.dumps(STATUS if name==NAMES[2] else INTERPRETATION) if categorical else '')
    return result


def run(source,manifest_path,s):
    root=Path(__file__).resolve().parents[2]
    source=local_path(source,root)
    manifest_path=local_path(manifest_path or root/'CSV files'/source.name/'run_manifest.json',root)
    manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
    dates=manifest['dates'];settings=manifest['settings']
    phases=[d['phase'] for d in dates]
    if [phases.count(k) for k in ('pre','event','post')]!=[8,1,8]:
        raise ValueError('RRI requires a completed 8 PRE / 1 EVENT / 8 POST run.')
    origin=pd.Timestamp(settings['EVENT_ORIGIN'])
    final_day=max((pd.Timestamp(d['reference'])-origin).total_seconds()/86400 for d in dates if d['phase']=='post')
    raster_dir=root/'RII'/'Raster Files';csv_dir=root/'RII'/'CSV'
    if any(p.exists() and any(p.iterdir()) for p in (raster_dir,csv_dir)):
        raise FileExistsError('RII output folders must be empty; archive the previous RII run before rerunning.')
    required=('Raw_ITSA_Metrics','Raw_Pixel_ITSA_Classes','Raw_ITSA_Reasons','Original_ITSA_Recovery_Classification',
              'Valid_Data_Coverage','Barangay_Owner')
    with ExitStack() as context:
        readers={k:context.enter_context(rasterio.open(source/(k+'.tif'))) for k in required}
        ref=readers['Raw_ITSA_Metrics']
        if not ref.crs or not ref.crs.is_projected: raise ValueError('Projected analysis grid required.')
        if any((ds.crs,ds.transform,ds.shape)!=(ref.crs,ref.transform,ref.shape) for ds in readers.values()):
            raise ValueError('Saved ITSA inputs do not share an exact grid.')
        bands={name:i+1 for i,name in enumerate(ref.descriptions)}
        needed=('beta0','beta1','beta2','beta3','pre_count','post_count','adaptive_loss_threshold')
        if any(k not in bands for k in needed): raise ValueError('Raw metrics lack required named ITSA bands.')
        if readers['Valid_Data_Coverage'].count!=17: raise ValueError('Expected 17 acquisition coverage bands.')
        if pd.Timestamp(ref.tags()['time_origin']) != origin:
            raise ValueError('Saved coefficient time origin differs from the manifest.')
        for band,date in enumerate(dates,1):
            tags=readers['Valid_Data_Coverage'].tags(band)
            if (pd.Timestamp(tags['master'])!=pd.Timestamp(date['master']) or
                    pd.Timestamp(tags['slave'])!=pd.Timestamp(date['slave'])):
                raise ValueError('Coverage acquisition metadata differs from the manifest.')
        c=SimpleNamespace(**settings);c.PROJECT_ROOT=root
        records,boundary_report=load_barangays(root/'Boundary'/'INFANTA_BOUNDARY.kml',ref.crs,c)
        raster_dir.mkdir(parents=True,exist_ok=True);csv_dir.mkdir(parents=True,exist_ok=True)
        outputs=create_rasters(context,raster_dir,ref,s)
        area=abs(ref.transform.a*ref.transform.e)*ref.crs.linear_units_factor[1]**2
        with TemporaryDirectory(prefix='.rri_stats_',dir=csv_dir) as scratch:
            summary=Summary(scratch,'Study area')
            for win in windows(ref.width,ref.height,s.WINDOW_SIZE):
                m={k:ref.read(bands[k],window=win,masked=True).filled(np.nan) for k in needed}
                raw=readers['Raw_Pixel_ITSA_Classes'].read(1,window=win)
                reasons=readers['Raw_ITSA_Reasons'].read(1,window=win)
                inside=reasons!=11
                inside&=raw!=0
                support=(m['pre_count']>=c.ITSA_MIN_PRE)&(m['post_count']>=c.ITSA_MIN_POST)&(reasons!=7)
                baseline=m['beta0'];event=baseline+m['beta2']
                final=event+(m['beta1']+m['beta3'])*final_day
                products=calculate(baseline,event,final,np.isin(raw,(1,2,3)),support,
                                   reasons!=3,inside,m['adaptive_loss_threshold'],s.MIN_DENOMINATOR,
                                   s.ZERO_TOLERANCE,s.ONE_TOLERANCE,no_loss=raw==4)
                products=(products[0].astype('float32'),*products[1:])
                for name,a in zip(NAMES,products):
                    ds=outputs[name];ds.write(np.where(np.isfinite(a),a,ds.nodata).astype(ds.dtypes[0]),1,window=win)
                original=readers['Original_ITSA_Recovery_Classification'].read(1,window=win)
                summary.add(products,original,inside)
            summary.export(csv_dir,area)
        paths={k:Path(ds.name) for k,ds in outputs.items()}
    for path in paths.values():
        with rasterio.open(path) as ds:
            if (ds.crs,ds.transform,ds.shape)!=(ref.crs,ref.transform,ref.shape): raise RuntimeError('RRI output grid mismatch.')
    barangay_rows=export_barangays(paths,source,records,raster_dir,csv_dir,area,s)
    pd.DataFrame(barangay_rows).to_csv(csv_dir/'RRI_Barangay_Summary.csv',index=False)
    (csv_dir/'RRI_Provenance.json').write_text(json.dumps(dict(source=str(source),manifest=str(manifest_path),
       dates=dates,boundary=boundary_report,formula='(Cfinal-(beta0+beta2))/(-beta2)',
       interpretation='final-observation adaptation, not maximum-post-recovery RRI',settings=
       {k:getattr(s,k) for k in ('MIN_DENOMINATOR','ZERO_TOLERANCE','ONE_TOLERANCE')}),indent=2),encoding='utf-8')
    print('RRI complete:',raster_dir)


def export_barangays(paths,source,records,raster_dir,csv_dir,area,s):
    rows=[]
    with ExitStack() as context:
        readers={k:context.enter_context(rasterio.open(p)) for k,p in paths.items()}
        coverage=context.enter_context(rasterio.open(source/'Valid_Data_Coverage.tif'))
        original=context.enter_context(rasterio.open(source/'Original_ITSA_Recovery_Classification.tif'))
        ref=readers[NAMES[0]]
        for record in records:
            try: crop=geometry_window(ref,[record.geometry]).intersection(Window(0,0,ref.width,ref.height))
            except WindowError: continue
            rd=raster_dir/'Barangay'/record.folder;cd=csv_dir/'Barangay'/record.folder
            rd.mkdir(parents=True);cd.mkdir(parents=True)
            with ExitStack() as local, TemporaryDirectory(prefix='.stats_',dir=cd) as scratch:
                writers={}
                for name,src in readers.items():
                    profile=src.profile.copy();profile.update(width=int(crop.width),height=int(crop.height),transform=window_transform(crop,ref.transform))
                    ds=local.enter_context(rasterio.open(rd/(name+'.tif'),'w',**profile));ds.update_tags(**src.tags())
                    ds.update_tags(barangay=record.name,barangay_code=record.code);writers[name]=ds
                summary=Summary(scratch,record.name);counts=np.zeros(17,dtype='int64');total=0
                # Same code-sorted ownership rule as the existing barangay workflow.
                index=records.index(record)+1
                for win in windows(int(crop.width),int(crop.height),s.WINDOW_SIZE):
                    sw=Window(crop.col_off+win.col_off,crop.row_off+win.row_off,win.width,win.height)
                    member=owners(records,(int(win.height),int(win.width)),window_transform(sw,ref.transform))==index
                    total+=int(member.sum())
                    products=[ds.read(1,window=sw,masked=True).filled(np.nan) if ds.dtypes[0]=='float32'
                              else ds.read(1,window=sw) for ds in readers.values()]
                    for name,a in zip(NAMES,products):
                        ds=writers[name];ds.write(np.where(member,a,ds.nodata).astype(ds.dtypes[0]),1,window=win)
                    counts+=(coverage.read(window=sw)==1)[:,member].sum(axis=1)
                    summary.add(products,original.read(1,window=sw),member)
                eligible=total>0 and bool((counts/max(total,1)>=.8).all())
                row=summary.export(cd,area,barangay=True)
                row.update(barangay_code=record.code,coverage_eligible=eligible,
                           minimum_date_valid_fraction=float(counts.min()/total) if total else 0,
                           interpretation='supported region summary' if eligible else 'descriptive available-pixel summary; below 80% date coverage')
                if not eligible:
                    for key in ('mean','median','std_population','minimum','maximum','q25','q75','percent_rri_ge_1'):
                        row[key]=np.nan
                    by_class=pd.read_csv(cd/'RRI_By_Recovery_Class.csv')
                    for key in ('mean','median','std_population','minimum','maximum','q25','q75'):
                        by_class[key]=np.nan
                    by_class['coverage_eligible']=False
                    by_class.to_csv(cd/'RRI_By_Recovery_Class.csv',index=False)
                pd.DataFrame([row]).to_csv(cd/'RRI_Statistics.csv',index=False)
                pd.DataFrame(dict(observation_index=np.arange(17),valid_pixels=counts,
                                  owned_pixels=total,valid_fraction=counts/max(total,1))).to_csv(cd/'RRI_Date_Coverage.csv',index=False)
                rows.append(row)
    return rows
