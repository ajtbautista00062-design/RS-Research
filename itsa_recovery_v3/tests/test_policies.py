"""Policy tests use synthetic arrays only; no research input rasters are accessed."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
sys.path.insert(1,str(HERE/"preserved_v1"))
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from rasterio.transform import from_origin
from rasterio.windows import Window
import config
from preparation import read_stack,read_window
from recovery import assess
from validation import variants
from test_model_recovery import observations

class Fake:
    def __init__(self,a):
        self.a=np.ma.array(a,dtype="float32",mask=False)
    def read(self,band,window,masked):
        r,c=int(window.row_off),int(window.col_off)
        return self.a[r:r+int(window.height),c:c+int(window.width)]

class PolicyTests(unittest.TestCase):
    def source(self):
        return SimpleNamespace(scales=(1.,),offsets=(0.,),name="SYNTHETIC")
    def test_valid_zero_retained_and_infinity_invalid(self):
        reader=config.reader_settings()
        a=np.array([[0,.5,np.inf]],dtype="float32")
        stack=read_stack([Fake(a)],[Fake(np.zeros_like(a))],[self.source()],Window(0,0,3,1),reader)
        self.assertEqual(stack[0,0,0],0)
        self.assertTrue(np.isnan(stack[0,0,2]))
    def test_zero_exclusion_is_separate_policy(self):
        a=np.zeros((1,1),dtype="float32")
        stack=read_stack([Fake(a)],[Fake(a)],[self.source()],Window(0,0,1,1),
                         config.reader_settings(exclude_zero=True))
        self.assertTrue(np.isnan(stack).all())
    def test_primary_unfilled_and_sensitivity_filled(self):
        a=np.full((9,9),.6,dtype="float32"); a[4,4]=np.nan
        ref=SimpleNamespace(width=9,height=9,transform=from_origin(0,9,1,1))
        readers=([Fake(a)],[Fake(a)])
        raw,flags=read_window(readers,[self.source()],ref,Window(0,0,9,9),config.reader_settings(),[])
        filled,provenance=read_window(readers,[self.source()],ref,Window(0,0,9,9),
                                     config.reader_settings(filled=True),[])
        self.assertTrue(np.isnan(raw[0,4,4]))
        self.assertEqual(provenance[0,4,4],2)
        self.assertAlmostEqual(filled[0,4,4],.6,places=6)
    def test_threshold_sensitivity_changes_model_class(self):
        obs=observations()
        stack=np.array([.7,.7,.7,.3,.62,.64,.66])[:,None,None]
        settings=config.analysis_settings(); settings.USE_UNCERTAINTY_GUARD=False
        settings.RETURN_TOLERANCE=.03
        strict=assess(stack,obs,config.EVENT_ORIGIN,settings)
        settings.RETURN_TOLERANCE=.07
        loose=assess(stack,obs,config.EVENT_ORIGIN,settings)
        self.assertEqual(strict.classes.item(),2)
        self.assertEqual(loose.classes.item(),1)
    def test_runtime_modules_import_without_running_pipeline(self):
        import workflow,reporting,raster_outputs
        self.assertTrue(callable(workflow.execute))
        self.assertTrue(callable(reporting.Summary))
    def test_full_stage_never_runs_implicitly(self):
        import run_analysis
        with self.assertRaises(ValueError):
            run_analysis.stage_gate("full",observations(8,8),None)
        with self.assertRaises(ValueError):
            run_analysis.stage_gate("pilot",observations(8,8),None)
    def test_resolved_zero_policy_needs_documented_evidence(self):
        import run_analysis
        with patch.object(config,"ZERO_SEMANTICS","valid"),patch.object(config,"ZERO_SEMANTICS_JUSTIFICATION",""):
            with self.assertRaises(ValueError):
                run_analysis.stage_gate("pilot",observations(),None)
    def test_sensitive_settings_are_not_mutated(self):
        from dataclasses import dataclass
        import pandas as pd
        @dataclass
        class Obs:
            phase: str
            reference: object
            start: object
            end: object
        simple=observations()
        obs=[Obs(o.phase,o.reference,o.reference-pd.Timedelta(days=12),o.reference) for o in simple]
        settings=config.analysis_settings()
        cases=variants(settings,obs)
        self.assertEqual(settings.LOSS_MIN,config.CALIBRATED_LOSS_MIN)
        self.assertEqual(settings.RETURN_TOLERANCE,.05)
        self.assertTrue(any(name=="midpoint_indexing" for name,_,_ in cases))

if __name__=="__main__":
    unittest.main()
