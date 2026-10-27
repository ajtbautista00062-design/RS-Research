"""Synthetic GeoTIFF roundtrip only; never reads the pilot rasters."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest,tempfile,shutil
from contextlib import ExitStack
from types import SimpleNamespace
import numpy as np
import rasterio
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from rasterio.windows import Window
import config
from recovery import assess
from raster_outputs import Outputs
from test_model_recovery import observations

class RasterTests(unittest.TestCase):
    def test_roundtrip_grid_nodata_quality_and_classes(self):
        settings=config.analysis_settings()
        obs=observations()
        data=np.tile(np.array([.7,.7,.7,.3,.5,.72,.94],dtype="float32")[:,None,None],(1,8,8))
        inside=np.ones((8,8),bool); inside[0]=False
        data[:,~inside]=np.nan
        result=assess(data,obs,config.EVENT_ORIGIN,settings,inside)
        workspace=Path(__file__).resolve().parents[1]
        folder=Path(tempfile.mkdtemp(prefix="synthetic_test_",dir=workspace)).resolve()
        self.assertTrue(folder.is_relative_to(workspace))
        try:
            with MemoryFile() as mem, ExitStack() as resources:
                ref=resources.enter_context(mem.open(driver="GTiff",height=8,width=8,count=1,dtype="float32",
                    crs="EPSG:3123",transform=from_origin(500000,1600000,14,14)))
                named=[SimpleNamespace(name=f"SYNTHETIC_{j}") for j in range(7)]
                output=Outputs(resources,folder/"rasters",ref,named,config,dict(test="synthetic"),full=True)
                output.write(result,data,np.isfinite(data).astype("uint8"),Window(0,0,8,8))
                output.close(); output.verify(); output.finalize()
                with rasterio.open(folder/"rasters"/"ITSA_Recovery_Classes.tif") as src:
                    np.testing.assert_array_equal(src.read(1),result.classes)
                    self.assertEqual(src.nodata,0)
                with rasterio.open(folder/"rasters"/"model_coefficients_diagnostics.tif") as src:
                    self.assertTrue((src.read()[:,0,:]==config.NODATA).all())
        finally:
            if not folder.is_relative_to(workspace):
                raise RuntimeError("Refusing cleanup outside versioned workspace.")
            shutil.rmtree(folder)

if __name__=="__main__":
    unittest.main()
