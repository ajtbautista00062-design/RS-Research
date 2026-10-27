"""Object support and spatial uncertainty policies: tests are not executed here."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
from unittest.mock import patch
import numpy as np
import config
from object_statistics import classify_region
from spatial_regularization import regularize
from test_model_recovery import observations
from test_region_regularization import record


class OperationalObjectTests(unittest.TestCase):
    def setUp(self):
        self.c=config.analysis_settings();self.c.EVIDENCE_MODE='operational'
        self.c.ZERO_UNRESOLVED=False;self.obs=observations()

    def test_one_pixel_object_has_real_series_not_impossible_three_pixel_rule(self):
        values=np.array([.7,.7,.7,.3,.4,.45,.5])[:,None]
        row=classify_region(values,self.obs,config.EVENT_ORIGIN,self.c,config,196.)
        self.assertEqual(row['class_code'],2)
        self.assertEqual(row['evidence_status'],2)
        self.assertEqual(row['n_fit'],6)

    def test_pixels_are_not_additional_temporal_observations(self):
        values=np.tile(np.array([.7,.7,.7,.3,.4,.45,.5])[:,None],(1,20))
        row=classify_region(values,self.obs,config.EVENT_ORIGIN,self.c,config,4000.)
        self.assertEqual(row['n_fit'],6)
        self.assertEqual(row['residual_df'],2)

    def test_some_missing_members_do_not_invalidate_an_adequate_date(self):
        values=np.tile(np.array([.7,.7,.7,.3,.4,.45,.5])[:,None],(1,20))
        values[4,:3]=np.nan
        row=classify_region(values,self.obs,config.EVENT_ORIGIN,self.c,config,4000.)
        self.assertEqual(row['class_code'],2)
        self.assertEqual(row['date_5_count'],17)
        self.assertAlmostEqual(row['date_5_fraction'],.85)

    def test_heterogeneity_is_flagged_without_automatic_rejection(self):
        values=np.tile(np.array([.7,.7,.7,.3,.4,.45,.5])[:,None],(1,16))
        values[:,8:]+=.25
        row=classify_region(values,self.obs,config.EVENT_ORIGIN,self.c,config,3200.)
        self.assertEqual(row['class_code'],2)
        self.assertEqual(row['evidence_status'],2)
        with patch.object(config,'OBJECT_HETEROGENEITY_POLICY','strict'):
            strict=classify_region(values,self.obs,config.EVENT_ORIGIN,self.c,config,3200.)
        self.assertEqual(strict['class_code'],5)
        self.assertEqual(strict['exclusion_reason'],21)

    def test_mapping_unit_flag_does_not_erase_a_supported_class(self):
        rows,_=regularize({1:record(3,.3,196)},set(),config,lambda ids:None)
        self.assertEqual(rows[1]['class_code'],3)
        self.assertEqual(rows[1]['evidence_status'],2)

    def test_spatial_regularization_never_promotes_indeterminate(self):
        rows,_=regularize({1:record(5,.0,196),2:record(1,.3,2000)},
            {(1,2)},config,lambda ids:record(1,.9,2196))
        self.assertEqual(rows[1]['class_code'],5)
