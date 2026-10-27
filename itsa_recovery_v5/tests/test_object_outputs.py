"""Tiny synthetic round-trip test for a future authorized test run; not executed."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
sys.path.insert(1,str(Path(__file__).resolve().parents[1]/'preserved_v1'))
import tempfile
import unittest
from types import SimpleNamespace
import numpy as np
import rasterio
from rasterio.transform import from_origin
import config
from test_model_recovery import observations


class ObjectOutputTests(unittest.TestCase):
    def test_small_object_pipeline_preserves_grid_nodata_and_model_link(self):
        from object_pipeline import execute
        values={key:getattr(config,key) for key in dir(config) if key.isupper()}
        c=SimpleNamespace(**values,analysis_settings=config.analysis_settings)
        c.MIN_REGION_AREA_M2=16.;c.TARGET_REGION_AREA_M2=64.;c.MAX_REGION_AREA_M2=256.
        c.CLASSIFICATION_MODE='object_itsa';c.FIGURE_FORMATS=()
        c.SPATIAL_BARRIER_RASTER='';c.VERIFIED_BUILTUP_MASK=''
        obs=observations()
        for j,o in enumerate(obs):
            o.name=f'synthetic_{j}'
            o.start=o.reference;o.end=o.reference
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);pixel=root/'Pixel ITSA Reference';pixel.mkdir()
            (root/'Reports').mkdir();(root/'Time Series Graphs').mkdir()
            transform=from_origin(100,100,2,2)
            profile=dict(driver='GTiff',height=8,width=8,crs='EPSG:32651',transform=transform)
            stack=np.broadcast_to(np.array([.7,.7,.7,.3,.5,.72,.94])[:,None,None],(7,8,8)).copy()
            stack[:,:,7]=c.NODATA
            specifications={'input_coherence':(stack,'float32',c.NODATA),
                'exclusion_reason':(np.zeros((1,8,8)),'uint8',255),
                'ITSA_Recovery_Classes':(np.ones((1,8,8)),'uint8',0),
                'final_guard_halfwidth':(np.full((1,8,8),.001),'float32',c.NODATA)}
            specifications['exclusion_reason'][0][:,:,7]=1
            for name,(a,dtype,nodata) in specifications.items():
                with rasterio.open(pixel/f'{name}.tif','w',**profile,count=len(a),dtype=dtype,nodata=nodata) as dst:
                    dst.write(a.astype(dtype))
            execute(c,obs,root,pixel,{})
            folder=root/'Recovery Map Raster File'
            with rasterio.open(folder/'Spatial_Object_ID.tif') as ids,rasterio.open(folder/'ITSA_Recovery_Classes.tif') as classes:
                self.assertEqual(ids.transform,transform)
                self.assertEqual(classes.crs.to_epsg(),32651)
                self.assertEqual(classes.nodata,0)
                self.assertTrue(np.all(ids.read(1)[:,7]==0))
                self.assertTrue(np.all(classes.read(1)[:,7]==0))
                self.assertTrue(np.isin(classes.read(1),range(6)).all())
            import pandas as pd
            table=pd.read_csv(root/'Reports'/'region_statistics.csv')
            self.assertIn('b2',table)
            self.assertIn('date_1_fraction',table)
            self.assertTrue(table.region_id.is_unique)
