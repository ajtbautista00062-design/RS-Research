"""Small software-test fixtures; not study observations or simulation results."""
import unittest
from types import SimpleNamespace
import numpy as np
from sarflow.itsa import raster_bands, fit_complete

class PilotModelTests(unittest.TestCase):
    def setUp(self):
        t=np.array([-60.,-36.,-12.,6.,18.,30.,42.])
        post=np.array([0,0,0,0,1,1,1])
        self.x=np.column_stack((np.ones(7),t,post,t*post))
        self.x[3]=0
        self.beta=np.array([.65,.001,-.2,.002])
        self.y=self.x@self.beta
        self.y[3]=.25

    def test_event_has_no_influence(self):
        expected=fit_complete(self.y,self.x)
        self.y[3]=np.nan
        np.testing.assert_allclose(fit_complete(self.y,self.x),expected,atol=1e-12)

    def test_independent_fit_and_post_slope(self):
        bands=raster_bands(self.y[:,None,None],self.x)[:,0,0]
        use=self.x[:,0]==1
        independent=np.linalg.lstsq(self.x[use],self.y[use],rcond=None)[0]
        np.testing.assert_allclose(bands[7:11],independent,atol=1e-6)
        self.assertAlmostEqual(bands[11],.003,places=6)
        self.assertEqual(bands[13],6)

    def test_missing_event_retains_model(self):
        self.y[3]=np.nan
        bands=raster_bands(self.y[:,None,None],self.x)[:,0,0]
        self.assertTrue(np.isfinite(bands[7:11]).all())

    def test_missing_post_blocks_fit_not_count(self):
        self.y[-1]=np.nan
        bands=raster_bands(self.y[:,None,None],self.x)[:,0,0]
        self.assertTrue(np.isnan(bands[7:11]).all())
        self.assertEqual(bands[13],5)

    def test_constant_r_squared_undefined(self):
        y=np.ones((7,1,1))*.5
        bands=raster_bands(y,self.x)[:,0,0]
        self.assertTrue(np.isnan(bands[-1]))
        self.assertLess(bands[-2],1e-10)

if __name__=="__main__":
    unittest.main()
