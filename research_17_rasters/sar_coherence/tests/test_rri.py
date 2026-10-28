import unittest
import numpy as np
from sarflow.rri import calculate


class RRITests(unittest.TestCase):
    def run_value(self,final,**changes):
        kw=dict(baseline=np.array([.6]),event=np.array([.4]),final=np.array([final]),
                detected=True,supported=True,model_valid=True,inside=True,loss_threshold=.03)
        kw.update(changes)
        return calculate(**kw)

    def test_recovery_values(self):
        for final,expected in ((.6,1),(.5,.5),(.4,0),(.3,-.5),(.8,2)):
            value,_,status,_,_=self.run_value(final)
            self.assertAlmostEqual(value[0],expected)
            self.assertEqual(status[0],1)

    def test_negligible_loss(self):
        out=self.run_value(.6,event=np.array([.599]))
        self.assertTrue(np.isnan(out[0][0]));self.assertEqual(out[2][0],4)

    def test_zero_denominator(self):
        self.assertTrue(np.isnan(self.run_value(.6,event=np.array([.6]))[0][0]))

    def test_invalid_model(self):
        for final in (np.nan,1.1):
            out=self.run_value(final)
            self.assertTrue(np.isnan(out[0][0]));self.assertEqual(out[2][0],5)

    def test_missing_support(self):
        out=self.run_value(.6,supported=False)
        self.assertEqual(out[2][0],3);self.assertTrue(np.isnan(out[0][0]))

    def test_no_loss_and_outside(self):
        self.assertEqual(self.run_value(.6,detected=False,no_loss=True)[2][0],2)
        self.assertEqual(self.run_value(.6,inside=False)[2][0],0)

    def test_ambiguous_not_undisturbed(self):
        self.assertEqual(self.run_value(.6,detected=False,no_loss=False)[2][0],6)

    def test_tolerances_and_no_clamping(self):
        for final,code in ((.39,2),(.595,4),(.8,5),(.3,1),(.5,3)):
            self.assertEqual(self.run_value(final)[1][0],code)

    def test_invalid_tolerance(self):
        with self.assertRaises(ValueError): self.run_value(.6,zero_tolerance=.6)


if __name__=='__main__': unittest.main()
