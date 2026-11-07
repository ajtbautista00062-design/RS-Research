"""Windowed GeoTIFF writer with explicit units, codebooks, and immutable sources."""
import json
import numpy as np
import rasterio
from .recovery import CLASSES, QUALITY, REASONS
from .recovery_io import METRICS, COLORS, unit
from .rri import INTERPRETATION, STATUS


class RasterProducts:
    def __init__(self, context, folders, profile, c, obs):
        self.outputs = {}
        self.c = c
        definitions = {
            'ITSA_Recovery_Classes':('uint8',0,1), 'Original_ITSA_Recovery_Classes':('uint8',0,1),
            'Raw_Pixel_ITSA_Classes':('uint8',0,1), 'Original_Raw_Pixel_ITSA_Classes':('uint8',0,1),
            'ITSA_Exclusion_Reasons':('uint8',255,1),'ITSA_Quality_Flags':('uint16',65535,1),
            'Raw_ITSA_Reasons':('uint8',255,1),'Raw_ITSA_Quality_Flags':('uint16',65535,1),
            'ITSA_Metrics':('float32',c.NODATA,len(METRICS)),
            'Raw_ITSA_Metrics':('float32',c.NODATA,len(METRICS)),
            'Spatial_Object_ID':('uint32',0,1),'Barangay_Owner':('uint16',0,1),
            'Study_Mask':('uint8',0,1), 'Aligned_Coherence_Observations':('float32',c.NODATA,len(obs)),
            'Valid_Data_Coverage':('uint8',255,len(obs)),
            'ITSA_Baseline':('float32',c.NODATA,1),'ITSA_Disturbance_Magnitude':('float32',c.NODATA,1),
            'ITSA_Final_Deficit':('float32',c.NODATA,1),'ITSA_Recovery_Progress':('float32',c.NODATA,1),
            'ITSA_Recovery_Fraction':('float32',c.NODATA,1),'ITSA_Crossing_Days':('float32',c.NODATA,1),
        }
        rri = {'ITSA_RRI':('float32',c.NODATA,1),'RRI_Interpretation_Classes':('uint8',0,1),
               'RRI_Validity_Status':('uint8',255,1),'RRI_Disturbance_Magnitude':('float32',c.NODATA,1),
               'RRI_Recovery_Magnitude':('float32',c.NODATA,1)}
        self.create(context,folders['ITSA'].rasters,profile,definitions,'ITSA',obs)
        self.create(context,folders['RRI'].rasters,profile,rri,'RRI',obs)

    def create(self,context,folder,profile,definitions,method,obs):
        for name,(dtype,nodata,count) in definitions.items():
            p = profile.copy()
            p.pop('source_width',None)
            p.update(driver='GTiff',count=count,dtype=dtype,nodata=nodata,tiled=True,
                     blockxsize=256,blockysize=256,compress='deflate',predictor=3 if dtype=='float32' else 2,
                     BIGTIFF='IF_SAFER')
            ds = context.enter_context(rasterio.open(folder/(name+'.tif'),'w',**p))
            self.outputs[name] = ds
            ds.update_tags(method=method,classification_mode=self.c.CLASSIFICATION_MODE,
                           scientific_status='experimental; unknown coherence estimator bias floor; no field accuracy validation',
                           source_semantics='unconfirmed; metadata masks honored; undeclared zeros valid',
                           spatial_filling='none',event_observation='excluded from pure-phase OLS',
                           original='unmodified local Codes classifier on identical selected inputs and regions')
            if name.endswith('Classes') and method=='ITSA':
                ds.write_colormap(1,COLORS)
                ds.update_tags(codebook=json.dumps(CLASSES))
            if name.endswith('Reasons'):
                ds.update_tags(codebook=json.dumps(REASONS))
            if name.endswith('Flags'):
                ds.update_tags(codebook=json.dumps(QUALITY))
            if name=='RRI_Interpretation_Classes':
                ds.update_tags(codebook=json.dumps(INTERPRETATION))
            if name=='RRI_Validity_Status':
                ds.update_tags(codebook=json.dumps(STATUS),note='0 is valid outside/unavailable code; NoData is 255')
            for i in range(count):
                label = METRICS[i] if name.endswith('ITSA_Metrics') else obs[i].name if count==len(obs) else name
                ds.set_band_description(i+1,label)
                units = unit(label) if name.endswith('ITSA_Metrics') else 'code' if dtype!='float32' else 'coherence'
                if name in ('ITSA_RRI','ITSA_Recovery_Fraction'):
                    units='dimensionless; unclipped'
                if name=='ITSA_Crossing_Days':
                    units='days since intervention origin'
                ds.update_tags(i+1,units=units)

    def write(self,name,values,window):
        ds = self.outputs[name]
        a = np.asarray(values)
        if a.ndim==2:
            a=a[None]
        ds.write(np.where(np.isfinite(a),a,ds.nodata).astype(ds.dtypes[0]),window=window)

    def ensure_gates(self,context,folder,profile,names):
        if 'ITSA_Gate_Outcomes' in self.outputs:
            return
        definition={'ITSA_Gate_Outcomes':('uint8',255,len(names))}
        self.create(context,folder,profile,definition,'ITSA',[])
        ds=self.outputs['ITSA_Gate_Outcomes']
        for i,name in enumerate(names,1):
            ds.set_band_description(i,name)
            ds.update_tags(i,units='0 false / 1 true / 255 outside')
