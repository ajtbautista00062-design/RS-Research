"""Named barangay clips retain municipal model values; descriptive fits stay separate."""
from contextlib import ExitStack
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_window
from rasterio.windows import Window, transform
from .study_layout import folders
from .study_figures import itsa_graph, map_classes, rri_graph, save
from .study_tables import itsa_tables, rri_tables
from .recovery import assess, CLASSES
from .recovery_io import METRICS


def clip(sources, member, crop, folder, record):
    for name,src in sources.items():
        p=src.profile.copy()
        p.update(height=int(crop.height),width=int(crop.width),transform=transform(crop,src.transform))
        a=src.read(window=crop)
        a[:,~member]=src.nodata
        with rasterio.open(folder/(name+'.tif'),'w',**p) as dst:
            try:
                dst.write_colormap(1,src.colormap(1))
            except ValueError:
                pass
            dst.write(a)
            dst.update_tags(**src.tags())
            dst.update_tags(barangay=record.name,barangay_code=record.code,model_scope='clipped municipal fit; no barangay refit')
            for i,label in enumerate(src.descriptions,1):
                dst.set_band_description(i,label)
                dst.update_tags(i,**src.tags(i))


def time_rows(stack,member,obs):
    rows=[]
    total=int(member.sum())
    for a,o in zip(stack,obs):
        x=a[member&np.isfinite(a)]
        q=np.quantile(x,[.25,.5,.75]) if len(x) else [np.nan]*3
        rows.append(dict(observation=o.name,phase=o.phase,pair_start=o.start,pair_end=o.end,
                         reference=o.reference,pair_days=(o.end-o.start).days,valid_pixels=len(x),
                         study_pixels=total,valid_fraction=len(x)/total if total else 0,
                         q25=q[0],median=q[1],q75=q[2],support='per-date available cohort'))
    return rows


def comparison_graph(old,new,member,directory,title,c):
    import matplotlib.pyplot as plt
    a=np.bincount(old[member],minlength=6)
    b=np.bincount(new[member],minlength=6)
    fig,ax=plt.subplots(figsize=(11,5),layout='constrained')
    ix=np.arange(6)
    ax.bar(ix-.2,a,.4,label='Original rule on selected inputs',color='.65')
    ax.bar(ix+.2,b,.4,label='Corrected observability rule',color='#245a8d')
    ax.set_xticks(ix,['NoData','Recovered','Partial','Not recovered','No detectable\nloss','Indeterminate'])
    ax.set(title=title+' fixed population comparison',ylabel='Study pixels')
    ax.legend(frameon=False,fontsize=8)
    save(fig,directory,'Original_to_Corrected_Class_Counts',c)


def zone_outputs(root,name,record,crop,member,sources,data,c,do_clip=True):
    itsa=folders(root,'ITSA',record.folder if record else None)
    rri=folders(root,'RRI',record.folder if record else None)
    grouped={m:{k:s for k,s in sources.items() if ('RRI' in k)==(m=='RRI')} for m in ('ITSA','RRI')}
    if do_clip:
        clip(grouped['ITSA'],member,crop,itsa.rasters,record)
        clip(grouped['RRI'],member,crop,rri.rasters,record)
    sl=np.s_[int(crop.row_off):int(crop.row_off+crop.height),int(crop.col_off):int(crop.col_off+crop.width)]
    stack=data.stack[:,sl[0],sl[1]]
    rows=time_rows(stack,member,data.observations)
    pd.DataFrame(rows).to_csv(itsa.csv/'Coherence_Time_Series.csv',index=False)
    series=np.array([r['median'] if r['valid_fraction']>=c.MIN_PER_DATE_VALID_COVERAGE else np.nan for r in rows])
    result=assess(series[:,None],data.observations,c)
    itsa_graph(rows,result,data.observations,itsa.graphs,name,c)
    read=lambda k:sources[k].read(1,window=crop)
    classes=read('ITSA_Recovery_Classes')
    original=read('Original_ITSA_Recovery_Classes')
    raw=read('Raw_Pixel_ITSA_Classes')
    metrics_array=sources['ITSA_Metrics'].read(window=crop,masked=True).filled(np.nan)
    metrics=dict(zip(METRICS,metrics_array))
    area=abs(data.profile['transform'].a*data.profile['transform'].e)*data.profile['crs'].linear_units_factor[1]**2
    itsa_tables(member,classes,read('ITSA_Exclusion_Reasons'),read('ITSA_Quality_Flags'),metrics,
                stack,data.observations,result,itsa.csv,name,area,c,series)
    map_classes(classes,member,itsa.graphs,name+' object ITSA recovery',c)
    comparison_graph(original,classes,member,itsa.graphs,name,c)
    value=sources['ITSA_RRI'].read(1,window=crop,masked=True).filled(np.nan)
    category,status=read('RRI_Interpretation_Classes'),read('RRI_Validity_Status')
    rri_tables(member,value,category,status,raw,rri.csv,area,np.isfinite(stack),c)
    map_classes(category,member,rri.graphs,name+' independent RRI interpretation',c,'RRI')
    rri_graph(value,category,status,member,rri.graphs,name,c)
    pd.DataFrame(rows).to_csv(rri.csv/'Source_Date_Coverage.csv',index=False)
    if record:
        targets=pd.read_csv(root/'Validation'/'Target_Assessment.csv')
        targets[targets.zone==name].to_csv(itsa.csv/'Before_After_Target_Assessment.csv',index=False)


def export_all(data,root,paths,c):
    with ExitStack() as context:
        sources={k:context.enter_context(rasterio.open(p)) for k,p in paths.items()}
        shape=data.inside.shape
        zone_outputs(root,'Overall Infanta',None,Window(0,0,shape[1],shape[0]),data.inside,sources,data,c,False)
        report=[]
        for i,record in enumerate(data.barangays,1):
            rr,cc=np.where(data.ownership==i)
            if not len(rr):
                report.append(dict(code=record.code,name=record.name,status='no owned reference-grid cells'))
                continue
            crop=Window(int(cc.min()),int(rr.min()),int(np.ptp(cc)+1),int(np.ptp(rr)+1))
            member=data.ownership[int(crop.row_off):int(crop.row_off+crop.height),int(crop.col_off):int(crop.col_off+crop.width)]==i
            zone_outputs(root,record.name,record,crop,member,sources,data,c)
            report.append(dict(code=record.code,name=record.name,status='exported',pixels=int(member.sum()),folder=record.folder))
            print('Exported barangay',record.name,flush=True)
        pd.DataFrame(report).to_csv(root/'Validation'/'Barangay_Export_Status.csv',index=False)
