"""V5 operational rules: synthetic specifications, written but not executed."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
import numpy as np
from scipy.stats import t as student_t
import config
from recovery import assess
from test_model_recovery import observations


class OperationalRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.c=config.analysis_settings();self.c.EVIDENCE_MODE='operational'
        self.c.ZERO_UNRESOLVED=False;self.c.COMPARABILITY_UNVERIFIED=False
        self.obs=observations()

    def case(self,post,pre=(.7,.7,.7),event=.3):
        stack=np.array([*pre,event,*post],float)[:,None,None]
        return assess(stack,self.obs,config.EVENT_ORIGIN,self.c)

    def test_recovered(self):
        result=self.case([.5,.72,.94])
        self.assertEqual(result.classes.item(),1)
        self.assertGreater(result.metrics['modeled_recovery_progress'].item(),0)

    def test_partial(self):
        self.assertEqual(self.case([.4,.45,.5]).classes.item(),2)

    def test_negligible_recovery(self):
        self.assertEqual(self.case([.4,.4,.4]).classes.item(),3)

    def test_further_degradation(self):
        self.assertEqual(self.case([.5,.45,.4]).classes.item(),3)

    def test_no_loss_supported_by_observations(self):
        self.assertEqual(self.case([.7,.7,.7],event=.7).classes.item(),4)

    def test_nonsignificant_level_change_can_have_provisional_recovery(self):
        result=self.case([.5,.72,.94],pre=(.72,.65,.72))
        statistic=abs(result.model.beta[2].item())/result.metrics['level_change_se'].item()
        self.assertLess(statistic,student_t.ppf(.975,result.model.df.item()))
        self.assertEqual(result.classes.item(),1)
        self.assertEqual(result.metrics['evidence_status'].item(),2)
        self.assertEqual(result.metrics['guarded_reference_class'].item(),5)

    def test_nonsignificant_but_persistently_degraded(self):
        result=self.case([.4,.4,.4],pre=(.75,.65,.75))
        self.assertEqual(result.classes.item(),3)
        self.assertEqual(result.metrics['evidence_status'].item(),2)

    def test_confidence_alpha_does_not_change_operational_class(self):
        self.c.DETECTION_ALPHA=.01
        strict=self.case([.4,.45,.5],pre=(.72,.65,.72))
        self.c.DETECTION_ALPHA=.20
        loose=self.case([.4,.45,.5],pre=(.72,.65,.72))
        self.assertEqual(strict.classes.item(),loose.classes.item())

    def test_effect_near_detection_boundary_is_genuinely_ambiguous(self):
        level=.7-self.c.LOSS_MIN
        result=self.case([level,level,level])
        self.assertEqual(result.classes.item(),5)
        self.assertEqual(result.exclusion.item(),12)

    def test_contradictory_event_stays_indeterminate(self):
        result=self.case([.4,.45,.5],event=.95)
        self.assertEqual(result.classes.item(),5)
        self.assertEqual(result.exclusion.item(),14)

    def test_baseline_deficit_is_not_the_same_as_level_change(self):
        result=self.case([.60,.55,.50],pre=(.9,.8,.7),event=.8)
        self.assertEqual(result.classes.item(),5)
        self.assertEqual(result.exclusion.item(),24)

    def test_missing_event_can_use_modeled_and_observed_early_loss(self):
        result=self.case([.4,.45,.5],event=np.nan)
        self.assertEqual(result.classes.item(),2)
        self.assertEqual(result.metrics['event_agreement_code'].item(),0)
        self.c.REQUIRE_EVENT=True
        self.assertEqual(self.case([.4,.45,.5],event=np.nan).exclusion.item(),3)

    def test_return_by_first_post_does_not_invent_crossing_time(self):
        result=self.case([.68,.72,.76])
        self.assertEqual(result.classes.item(),1)
        self.assertEqual(result.metrics['crossing_status'].item(),2)
        self.assertTrue(np.isnan(result.metrics['recovery_duration_days'].item()))

    def test_positive_slope_without_initial_disturbance_is_not_recovery(self):
        self.assertEqual(self.case([.71,.72,.73],event=.7).classes.item(),4)

    def test_hard_model_failure_has_explicit_reason(self):
        result=self.case([.1,.9,.1])
        self.assertEqual(result.classes.item(),5)
        self.assertEqual(result.exclusion.item(),17)

    def test_incomplete_longer_series_uses_actual_supported_dates(self):
        self.obs=observations(3,6)
        result=self.case([.30,np.nan,.36,.39,.42,.45])
        self.assertEqual(result.classes.item(),2)
        self.assertEqual(result.metrics['valid_post_count'].item(),5)
        self.assertEqual(result.metrics['fit_observation_count'].item(),8)

    def test_missing_late_confirmation_remains_indeterminate(self):
        result=self.case([.4,np.nan,.5])
        self.assertEqual(result.classes.item(),5)
        self.assertEqual(result.exclusion.item(),4)  # only two valid POST dates

    def test_structural_missingness_and_outside_remain_nodata(self):
        data=np.full((7,1,1),np.nan)
        result=assess(data,self.obs,config.EVENT_ORIGIN,self.c)
        self.assertEqual(result.classes.item(),0)
        self.assertEqual(result.exclusion.item(),9)

    def test_missing_last_confirmation_has_separate_reason(self):
        self.obs=observations(3,4)
        result=self.case([.3,.35,.4,np.nan])
        self.assertEqual(result.classes.item(),5)
        self.assertEqual(result.exclusion.item(),19)

    def test_late_model_observation_contradiction_is_not_recovered(self):
        result=self.case([.4,.6,1.])
        self.assertEqual(result.classes.item(),5)
        self.assertEqual(result.exclusion.item(),15)
