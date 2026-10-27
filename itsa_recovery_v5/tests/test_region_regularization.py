"""Small synthetic adjacency tests. These tests have not been executed."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
from types import SimpleNamespace
from unittest.mock import Mock
import numpy as np
import config
from spatial_regularization import boundary_pairs,regularize,compatible
from spatial_settings import validate


def record(code=1,score=.3,area=1000,group=1,protected=False):
    return dict(class_code=code,evidence_score=score,area_m2=area,stratum=group,
        feature_min=[.6,.01,.01],feature_max=[.62,.02,.02],object_flags='',
        protected_high_confidence=protected,quality_flags=0,object_decision_flags=0)


class RegionRegularizationTests(unittest.TestCase):
    def setUp(self):
        from unittest.mock import patch
        historical=patch.object(config,'MMU_POLICY','defer_low_confidence')
        historical.start();self.addCleanup(historical.stop)
    def test_isolated_low_confidence_artifact_is_deferred(self):
        rows,mapping=regularize({1:record()},set(),config,Mock())
        self.assertEqual(rows[1]['class_code'],5)
        self.assertEqual(mapping[1],1)

    def test_small_high_confidence_degraded_region_survives(self):
        rows,_=regularize({1:record(3,.95,200),2:record(1,.95,3000)},
                           {(1,2)},config,Mock())
        self.assertEqual(rows[1]['class_code'],3)

    def test_strong_pixel_evidence_protects_small_region(self):
        rows,_=regularize({1:record(3,.3,200,protected=True)},set(),config,Mock())
        self.assertEqual(rows[1]['class_code'],3)

    def test_indeterminate_is_never_promoted(self):
        refit=Mock(return_value=record(1,.9,2000))
        rows,_=regularize({1:record(5),2:record()}, {(1,2)},config,refit)
        self.assertEqual(rows[1]['class_code'],5)
        refit.assert_not_called()

    def test_same_class_merge_requires_supporting_refit(self):
        merged=record(2,.4,2000)
        rows,mapping=regularize({1:record(2),2:record(2)}, {(1,2)},config,
                                Mock(return_value=merged))
        self.assertEqual(mapping,{1:1,2:1})
        self.assertEqual(rows[1]['class_code'],2)
        self.assertEqual(rows[2]['area_m2'],2000)

    def test_unsupported_refit_cannot_assign_recovery(self):
        rows,mapping=regularize({1:record(2),2:record(2)}, {(1,2)},config,
                                Mock(return_value=record(1,.9,2000)))
        self.assertEqual(mapping,{1:1,2:2})
        self.assertNotIn(1,[row['class_code'] for row in rows.values()])

    def test_nodata_and_diagonal_contact_do_not_connect(self):
        self.assertEqual(boundary_pairs(np.array([[1,0],[0,2]])),set())
        self.assertEqual(boundary_pairs(np.array([[1,2],[0,0]])),{(1,2)})

    def test_physical_boundaries_and_maximum_area(self):
        self.assertFalse(compatible(record(group=1),record(group=2),config))
        self.assertFalse(compatible(record(area=30000),record(area=30000),config))
        other=record();other['feature_max']=[.9,.02,.02]
        self.assertFalse(compatible(record(),other,config))

    def test_none_parameters_fail_before_processing(self):
        values={key:getattr(config,key) for key in dir(config) if key.isupper()}
        values['MIN_REGION_AREA_M2']=None
        with self.assertRaises(ValueError):
            validate(SimpleNamespace(**values))

    def test_degree_grid_area_is_rejected(self):
        from affine import Affine
        grid=SimpleNamespace(crs='EPSG:4326',transform=Affine.scale(.001,-.001),
                             width=10,height=10)
        with self.assertRaises(ValueError):
            validate(config,grid)
