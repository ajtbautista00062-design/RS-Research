"""Object-series/model dependence and conservative aggregation tests; not executed."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
import numpy as np
import config
from object_statistics import aggregate,classify_region
from test_model_recovery import observations


class ObjectITSATests(unittest.TestCase):
    def setUp(self):
        self.settings=config.analysis_settings()
        self.settings.ZERO_UNRESOLVED=False
        self.settings.COMPARABILITY_UNVERIFIED=False
        self.obs=observations()

    def region(self,post,event=.3):
        y=np.array([.7,.7,.7,event,*post])
        return classify_region(np.tile(y[:,None],(1,16)),self.obs,config.EVENT_ORIGIN,
                               self.settings,config,3200.)

    def test_exact_medians_counts_fractions_and_missing_dates(self):
        values=np.array([[.1,.2,.3,.4,.5],[.5,np.nan,np.nan,np.nan,.7]])
        median,count,fraction,mad=aggregate(values,3,.8)
        self.assertAlmostEqual(median[0],.3)
        self.assertTrue(np.isnan(median[1]))
        np.testing.assert_array_equal(count,[5,2])
        np.testing.assert_allclose(fraction,[1,.4])
        self.assertGreater(mad[0],0)

    def test_contiguous_recovered_neighborhood_uses_itsa(self):
        row=self.region([.5,.72,.94])
        self.assertEqual(row['class_code'],1)
        self.assertTrue(np.isfinite(row['b2']))
        self.assertGreater(row['modeled_recovery_progress'],0)

    def test_partial_neighborhood(self):
        self.assertEqual(self.region([.4,.45,.5])['class_code'],2)

    def test_degraded_neighborhood(self):
        self.assertEqual(self.region([.5,.45,.4])['class_code'],3)

    def test_missing_temporal_evidence_is_indeterminate(self):
        self.assertEqual(self.region([.4,np.nan,.6])['class_code'],5)

    def test_same_endpoint_different_itsa_trajectory(self):
        self.settings.USE_UNCERTAINTY_GUARD=False
        self.assertEqual(self.region([.5,.5,.5])['class_code'],3)
        self.assertEqual(self.region([.3,.5,.5])['class_code'],2)

    def test_strong_degraded_pixel_evidence_is_not_concealed(self):
        y=np.array([.7,.7,.7,.3,.5,.72,.94])
        row=classify_region(np.tile(y[:,None],(1,16)),self.obs,config.EVENT_ORIGIN,
            self.settings,config,3200.,np.array([3]+[1]*15),np.ones(16))
        self.assertEqual(row['class_code'],5)
        self.assertIn('conflicting_high_confidence',row['object_flags'])

    def test_spatial_heterogeneity_does_not_become_precision(self):
        from unittest.mock import patch
        values=np.tile(np.array([.7,.7,.7,.3,.4,.45,.5])[:,None],(1,16))
        values[:,8:]+= .25
        with patch.object(config,'OBJECT_HETEROGENEITY_POLICY','strict'):
            row=classify_region(values,self.obs,config.EVENT_ORIGIN,self.settings,config,3200.)
        self.assertEqual(row['class_code'],5)
