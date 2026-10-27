"""Controlled synthetic tests ONLY. These are never reported as research observations."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
from types import SimpleNamespace
from dataclasses import replace
from unittest.mock import patch
import numpy as np
import pandas as pd
import config
from model import design,fit
from recovery import assess
from validation import independent_check

def observations(npre=3,npost=3,reference_shift=0):
    times=list(np.linspace(-60,-12,npre))+[6]+list(np.linspace(12,36,npost))
    origin=pd.Timestamp(config.EVENT_ORIGIN)
    return [SimpleNamespace(phase="pre" if i<npre else "event" if i==npre else "post",
        reference=origin+pd.Timedelta(days=float(t)+reference_shift)) for i,t in enumerate(times)]

class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.c=config.analysis_settings()
        self.c.EVIDENCE_MODE='guarded'  # preserve historical-rule regression tests
        self.c.REQUIRE_COMPLETE=True
        self.c.ZERO_UNRESOLVED=False
        self.c.COMPARABILITY_UNVERIFIED=False
        self.obs=observations()
    def run_case(self,post,event=.3):
        stack=np.array([.7,.7,.7,event,*post],float)[:,None,None]
        return assess(stack,self.obs,config.EVENT_ORIGIN,self.c),stack

    def test_complete_recovery_and_supported_crossing(self):
        result,_=self.run_case([.5,.72,.94])
        self.assertEqual(result.classes.item(),1)
        self.assertTrue(12<result.metrics["recovery_duration_days"].item()<24)

    def test_partial(self):
        result,_=self.run_case([.4,.45,.5])
        self.assertEqual(result.classes.item(),2)
        self.assertTrue(np.isnan(result.metrics["recovery_duration_days"].item()))

    def test_no_recovery_and_degradation(self):
        for post in ([.4,.4,.4],[.5,.45,.4]):
            result,_=self.run_case(post)
            self.assertEqual(result.classes.item(),3)

    def test_no_loss(self):
        result,_=self.run_case([.7,.7,.7],event=.7)
        self.assertEqual(result.classes.item(),4)

    def test_insufficient_post(self):
        result,_=self.run_case([.4,np.nan,.6])
        self.assertEqual(result.classes.item(),5)
        self.assertEqual(result.exclusion.item(),4)

    def test_missing_pre_event_and_outside_distinct(self):
        _,stack=self.run_case([.4,.45,.5])
        for index,reason in ((0,2),(3,3)):
            data=stack.copy(); data[index]=np.nan
            result=assess(data,self.obs,config.EVENT_ORIGIN,self.c)
            self.assertEqual(result.classes.item(),5)
            self.assertEqual(result.exclusion.item(),reason)
        result=assess(stack,self.obs,config.EVENT_ORIGIN,self.c,inside=np.array([[False]]))
        self.assertEqual(result.classes.item(),0)
        self.assertEqual(result.exclusion.item(),1)
        self.assertEqual(result.quality.item(),65535)

    def test_conflicting_event_and_model(self):
        result,_=self.run_case([.7,.7,.7],event=.3)
        self.assertEqual(result.classes.item(),5)
        self.assertTrue(result.quality.item()&16)

    def test_conflicting_model_and_late_observations(self):
        self.c.USE_UNCERTAINTY_GUARD=False
        result,_=self.run_case([.4,.6,1.])
        self.assertEqual(result.classes.item(),5)
        self.assertTrue(result.quality.item()&32)

    def test_model_dependency_same_legacy_endpoint(self):
        # Both sequences have the same PRE, event, POST median and final two values.
        constant,_=self.run_case([.5,.5,.5])
        rising,_=self.run_case([.3,.5,.5])
        self.assertEqual(constant.classes.item(),3)
        self.assertEqual(rising.classes.item(),5)
        self.c.USE_UNCERTAINTY_GUARD=False
        rising,_=self.run_case([.3,.5,.5])
        self.assertEqual(rising.classes.item(),2)

    def test_independent_coefficients_and_predictions(self):
        result,stack=self.run_case([.4,.45,.5])
        self.assertEqual(independent_check(stack,self.obs,config.EVENT_ORIGIN,self.c,result,3),1)

    def test_event_excluded_from_regression_but_used_for_evidence(self):
        a,stack=self.run_case([.4,.45,.5])
        stack[3]=.95
        b=assess(stack,self.obs,config.EVENT_ORIGIN,self.c)
        np.testing.assert_allclose(a.model.beta,b.model.beta)
        self.assertNotEqual(a.classes.item(),b.classes.item())

    def test_candidate_beyond_period_not_duration(self):
        result,_=self.run_case([.4,.45,.5])
        self.assertGreater(result.metrics["candidate_crossing_days"].item(),36)
        self.assertTrue(result.quality.item()&256)
        self.assertTrue(np.isnan(result.metrics["recovery_duration_days"].item()))

    def test_recovered_before_first_post_has_no_observed_crossing(self):
        # Event loss is strong, model step is not: uncertainty class, no invented duration.
        result,_=self.run_case([.7,.72,.74])
        self.assertTrue(np.isnan(result.metrics["recovery_duration_days"].item()))

    def test_large_uncertainty_is_indeterminate(self):
        self.c.ERROR_FLOOR=.2
        result,_=self.run_case([.4,.45,.5])
        self.assertEqual(result.classes.item(),5)

    def test_seventeen_observations(self):
        obs=observations(8,8)
        times=design(obs,config.EVENT_ORIGIN)[2]
        y=np.array([.7 if o.phase=="pre" else .3 if o.phase=="event" else .28+.018*t
                    for o,t in zip(obs,times)])[:,None,None]
        result=assess(y,obs,config.EVENT_ORIGIN,self.c)
        self.assertEqual(result.model.n.item(),16)
        self.assertEqual(result.model.df.item(),12)
        self.assertEqual(result.classes.item(),1)

    def test_partial_support_is_explicit(self):
        obs=observations(8,8)
        stack=np.full((17,1,1),.7); stack[8]=.7; stack[0]=np.nan
        strict=assess(stack,obs,config.EVENT_ORIGIN,self.c)
        self.assertEqual(strict.classes.item(),5)
        self.c.REQUIRE_COMPLETE=False
        relaxed=assess(stack,obs,config.EVENT_ORIGIN,self.c)
        self.assertEqual(relaxed.classes.item(),4)

    def test_time_reference_sensitivity_is_real(self):
        result,stack=self.run_case([.4,.45,.5])
        shifted=assess(stack,observations(reference_shift=-6),config.EVENT_ORIGIN,self.c)
        self.assertNotAlmostEqual(result.metrics["candidate_crossing_days"].item(),
                                  shifted.metrics["candidate_crossing_days"].item())

    def test_neutral_event_does_not_veto_model_loss(self):
        strong,_=self.run_case([.4,.45,.5],event=.3)
        neutral,_=self.run_case([.4,.45,.5],event=.7)
        self.assertEqual(strong.classes.item(),2)
        self.assertEqual(neutral.classes.item(),2)
        self.assertTrue(neutral.quality.item()&2048)

    def test_uncertain_no_loss_is_not_class_four(self):
        self.c.ERROR_FLOOR=.2
        result,_=self.run_case([.7,.7,.7],event=.7)
        self.assertEqual(result.classes.item(),5)

    def test_unverified_all_zero_is_indeterminate(self):
        self.c.ZERO_UNRESOLVED=True
        stack=np.zeros((7,1,1))
        result=assess(stack,self.obs,config.EVENT_ORIGIN,self.c)
        self.assertEqual(result.classes.item(),5)
        self.assertEqual(result.exclusion.item(),10)

    def test_structurally_missing_is_nodata(self):
        result=assess(np.full((7,1,1),np.nan),self.obs,config.EVENT_ORIGIN,self.c)
        self.assertEqual(result.classes.item(),0)
        self.assertEqual(result.exclusion.item(),9)

    def test_absolute_sweep_not_clamped_by_old_spread_floor(self):
        self.c.ERROR_FLOOR=.001
        self.c.BASELINE_GUARD_FLOOR=.001
        self.c.LOSS_MIN=.03
        detected,_=self.run_case([.62,.62,.62],event=.62)
        self.c.LOSS_MIN=.10
        below,_=self.run_case([.62,.62,.62],event=.62)
        self.assertEqual(detected.classes.item(),3)
        self.assertEqual(below.classes.item(),4)

    def test_event_opposition_uses_fixed_calibrated_tolerance(self):
        self.c.USE_UNCERTAINTY_GUARD=False
        self.c.LOSS_MIN=.01
        small,_=self.run_case([.4,.45,.5],event=.72)
        self.c.LOSS_MIN=.10
        large,_=self.run_case([.4,.45,.5],event=.72)
        self.assertEqual(small.classes.item(),2)
        self.assertEqual(large.classes.item(),2)
        self.assertFalse(small.quality.item()&16)

    def test_quality_provenance(self):
        _,stack=self.run_case([.4,.45,.5])
        self.c.ZERO_UNRESOLVED=True
        result=assess(stack,self.obs,config.EVENT_ORIGIN,self.c,filled=np.array([[True]]))
        self.assertTrue(result.quality.item()&1)
        self.assertTrue(result.quality.item()&128)

if __name__=="__main__":
    unittest.main()
