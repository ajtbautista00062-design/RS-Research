"""Tiled GeoTIFF outputs with exact grid, units, quality codes and reread QA."""
from pathlib import Path
import json
import numpy as np
import rasterio
from recovery import METRICS,LABELS,COLORS,QUALITY,EXCLUSION,metric_unit,EVIDENCE_TAGS

def qml(path):
    from xml.etree.ElementTree import Element,SubElement,ElementTree
    root=Element("qgis",version="3.34")
    pipe=SubElement(root,"pipe")
    render=SubElement(pipe,"rasterrenderer",type="paletted",band="1",opacity="1",alphaBand="-1")
    palette=SubElement(render,"colorPalette")
    for code,label in LABELS.items():
        if code:
            r,g,b,_=COLORS[code]
            SubElement(palette,"paletteEntry",value=str(code),label=label,
                       color=f"#{r:02x}{g:02x}{b:02x}",alpha="255")
    SubElement(pipe,"rasterresampler",zoomedInResampler="nearestNeighbour",
               zoomedOutResampler="nearestNeighbour",maxOversampling="2")
    ElementTree(root).write(path.with_suffix(".qml"),encoding="utf-8",xml_declaration=True)

class Outputs:
    def __init__(self,resources,folder,ref,obs,c,metadata,full=True,provenance=False):
        folder.mkdir(parents=True,exist_ok=False)
        self.paths={}; self.writers={}; self.c=c
        self.grid=(ref.crs,ref.transform,ref.shape)
        profile=dict(driver="GTiff",crs=ref.crs,transform=ref.transform,width=ref.width,height=ref.height,
            tiled=True,blockxsize=256,blockysize=256,compress="deflate",BIGTIFF="IF_SAFER")
        specs={"ITSA_Recovery_Classes":("uint8",0,1)}
        if full:
            specs|={"Guarded_ITSA_Recovery_Classes":("uint8",0,1),"ITSA_Evidence_Status":("uint8",255,1)}
            specs|={name:("float32",c.NODATA,1) for name in METRICS}
            specs|={"quality_flags":("uint16",65535,1),"exclusion_reason":("uint8",255,1),
                    "input_coherence":("float32",c.NODATA,len(obs)),
                    "input_provenance":("uint8",0,len(obs)),
                    "model_coefficients_diagnostics":("float32",c.NODATA,9)}
        if provenance:
            specs|={"input_provenance":("uint8",0,len(obs)),"input_coherence":("float32",c.NODATA,len(obs))}
        self.specs=specs
        for name,(dtype,nodata,count) in specs.items():
            path=folder/f"{name}.partial.tif"
            dst=resources.enter_context(rasterio.open(path,"w",**profile,count=count,dtype=dtype,
                    nodata=nodata,predictor=3 if dtype=="float32" else 2))
            dst.update_tags(**metadata,product=name)
            unit=metric_unit(name) if name in METRICS else 'code'
            dst.update_tags(1,units=unit)
            for j in range(count):
                description=(obs[j].name if name in ("input_coherence","input_provenance") else
                    ("b0","b1_per_day","b2","b3_per_day","post_slope_per_day","RMSE","R_squared","n","residual_df")[j]
                    if name=="model_coefficients_diagnostics" else name)
                dst.set_band_description(j+1,description)
                band_unit=("coherence/day" if "per_day" in description else
                    "count" if description in ("n","residual_df") else
                    "dimensionless" if description=="R_squared" else unit)
                dst.update_tags(j+1,units=band_unit)
            if name.endswith("Recovery_Classes"):
                dst.write_colormap(1,COLORS)
                dst.update_tags(classes=json.dumps(LABELS))
            if name in EVIDENCE_TAGS:
                dst.update_tags(code_definitions=json.dumps(EVIDENCE_TAGS[name]))
            if name=='ITSA_Evidence_Status':
                dst.update_tags(code_definitions=json.dumps(EVIDENCE_TAGS['evidence_status']))
            if name=="input_provenance":
                dst.update_tags(code_0="missing/outside",code_1="valid aligned input",code_2="nearest-neighbor spatial estimate")
            if name=="quality_flags":
                dst.update_tags(bit_definitions=json.dumps(QUALITY),outside="65535")
            if name=="exclusion_reason":
                dst.update_tags(reason_definitions=json.dumps(EXCLUSION),outside="1; 255 reserved")
            self.paths[name]=path; self.writers[name]=dst
    def write(self,result,stack,provenance,window):
        m=result.model
        values={"ITSA_Recovery_Classes":result.classes,"quality_flags":result.quality,
            "Guarded_ITSA_Recovery_Classes":np.nan_to_num(result.metrics['guarded_reference_class'],nan=0),
            "ITSA_Evidence_Status":result.metrics['evidence_status'],
            "exclusion_reason":result.exclusion,**result.metrics,"input_coherence":stack,
            "input_provenance":provenance,
            "model_coefficients_diagnostics":np.concatenate((m.beta,
                (m.beta[1]+m.beta[3])[None],m.rmse[None],m.r_squared[None],m.n[None],m.df[None]))}
        values["model_coefficients_diagnostics"][:,result.exclusion==1]=np.nan
        for name,dst in self.writers.items():
            a=values[name]
            if a.ndim==2:
                a=a[None]
            dtype,nodata,_=self.specs[name]
            a=np.where(np.isfinite(a),a,nodata).astype(dtype)
            dst.write(np.ascontiguousarray(a),window=window)
    def empty(self,shape,window):
        for name,dst in self.writers.items():
            dtype,nodata,count=self.specs[name]
            fill=1 if name=="exclusion_reason" else nodata
            dst.write(np.full((count,*shape),fill,dtype=dtype),window=window)
    def close(self):
        for dst in self.writers.values():
            dst.close()
    def verify(self):
        for name,path in self.paths.items():
            with rasterio.open(path) as src:
                dtype,nodata,count=self.specs[name]
                if (src.crs,src.transform,src.shape)!=self.grid or src.count!=count or src.nodata!=nodata or src.dtypes!=(dtype,)*count:
                    raise RuntimeError(f"Output metadata verification failed: {name}")
                # Reread all tiles without modifying statistics.
                for _,window in src.block_windows(1):
                    a=src.read(window=window)
                    if not np.isfinite(a).all():
                        raise RuntimeError(f"Nonfinite output values: {name}")
                    if name.endswith("Recovery_Classes") and not np.isin(a,range(6)).all():
                        raise RuntimeError("Unknown output class code.")
                    if name=="exclusion_reason" and not np.isin(a,range(max(EXCLUSION)+1)).all():
                        raise RuntimeError("Unknown exclusion reason.")
    def finalize(self):
        for name,path in self.paths.items():
            target=path.with_name(path.name.replace(".partial.tif",".tif"))
            path.rename(target)
            if name.endswith("Recovery_Classes"):
                qml(target)
