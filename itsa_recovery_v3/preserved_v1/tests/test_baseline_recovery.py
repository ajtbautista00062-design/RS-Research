"""Small unit fixtures only; never execute the simulation or read source rasters."""
import unittest
from types import SimpleNamespace
import numpy as np
from sarflow.recovery import baseline_recovery, validate
from sarflow.zonal import ZoneStats

class BaselineRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.c = SimpleNamespace(CLASS_MIN_PRE=3, CLASS_MIN_POST=3, END_POST_COUNT=3,
            CONFIRM_POST_COUNT=2, LOSS_MIN=.10, PRE_SPREAD_FLOOR=.02,
            LOSS_SPREAD_MULTIPLIER=2., RETURN_TOLERANCE=.05)
        self.obs = [SimpleNamespace(phase=p) for p in ["pre"]*3+["event"]+["post"]*3]
        self.stack = np.array([.7,.7,.7,.4,.45,.66,.68],dtype="float32")[:,None,None]

    def test_recovered(self):
        codes, fraction = baseline_recovery(self.stack,self.obs,self.c)
        self.assertEqual(codes.item(),1)
        self.assertAlmostEqual(fraction.item(),(.66-.4)/.3,places=5)

    def test_partial_and_temporary_peak(self):
        self.stack[4:,0,0] = [.5,.7,.5]
        self.assertEqual(baseline_recovery(self.stack,self.obs,self.c)[0].item(),2)

    def test_further_decline(self):
        self.stack[4:,0,0] = [.3,.35,.38]
        codes,fraction=baseline_recovery(self.stack,self.obs,self.c)
        self.assertEqual(codes.item(),3)
        self.assertLess(fraction.item(),0)

    def test_no_loss_does_not_require_post(self):
        self.stack[3]=.68
        self.stack[4:]=np.nan
        codes,fraction=baseline_recovery(self.stack,self.obs,self.c)
        self.assertEqual(codes.item(),4)
        self.assertTrue(np.isnan(fraction.item()))

    def test_endpoint_missing_indeterminate(self):
        self.stack[-1]=np.nan
        codes,fraction=baseline_recovery(self.stack,self.obs,self.c)
        self.assertEqual(codes.item(),5)
        self.assertTrue(np.isnan(fraction.item()))

    def test_missing_pre_or_event_nodata(self):
        for i in (0,3):
            stack=self.stack.copy()
            stack[i]=np.nan
            self.assertEqual(baseline_recovery(stack,self.obs,self.c)[0].item(),0)

    def test_overshoot_retained(self):
        self.stack[4:]=.8
        codes,fraction=baseline_recovery(self.stack,self.obs,self.c)
        self.assertEqual(codes.item(),1)
        self.assertGreater(fraction.item(),1)

    def test_threshold_validation(self):
        validate(self.c)
        self.c.RETURN_TOLERANCE=.2
        with self.assertRaises(ValueError):
            validate(self.c)

    def test_zonal_denominators(self):
        zone=ZoneStats("test",("value",))
        codes=np.array([[0,1,2,3,4,5]],dtype="uint8")
        zone.add(np.zeros((1,1,6)),np.ones((1,6),dtype=bool),codes)
        row=zone.classes_row()
        self.assertEqual(row["assessed_count"],4)
        self.assertEqual(row["recovery_eligible_count"],3)
        self.assertAlmostEqual(row["recovered_share_eligible"],1/3)
        self.assertAlmostEqual(row["class_1_share_assessed"],.25)
        self.assertEqual(row["indeterminate_count"],1)
        self.assertTrue(np.isnan(ZoneStats("empty").classes_row()["recovered_share_eligible"]))

if __name__ == "__main__":
    unittest.main()
