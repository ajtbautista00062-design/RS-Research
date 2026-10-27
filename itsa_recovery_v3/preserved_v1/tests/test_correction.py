"""Unit fixtures for gap correction; never substitute for study observations."""
import unittest
from types import SimpleNamespace
import numpy as np
from sarflow.correction import fill_holes

class CorrectionTests(unittest.TestCase):
    def setUp(self):
        self.c=SimpleNamespace(FILL_SMALL_HOLES=True,FILL_MAX_HOLE_PIXELS=9,FILL_RADIUS_PIXELS=2)

    def test_enclosed_gap_and_original_preservation(self):
        a=np.full((15,15),.6,dtype="float32")
        a[7,7]=np.nan
        result,flag=fill_holes(a,self.c)
        self.assertEqual(result[7,7],a[7,6])
        self.assertEqual(flag[7,7],2)
        np.testing.assert_array_equal(result[np.isfinite(a)],a[np.isfinite(a)])

    def test_edge_connected_gaps_untouched(self):
        a=np.full((15,15),.6)
        a[:8,7]=np.nan
        result,flag=fill_holes(a,self.c)
        self.assertTrue(np.isnan(result[:8,7]).all())
        self.assertEqual((flag==2).sum(),0)

    def test_large_gap_untouched(self):
        a=np.full((15,15),.6)
        a[5:9,5:9]=np.nan
        self.assertTrue(np.isnan(fill_holes(a,self.c)[0][5:9,5:9]).all())

    def test_distance_limit_preserves_whole_hole(self):
        self.c.FILL_RADIUS_PIXELS=1
        a=np.full((15,15),.6)
        a[5:8,5:8]=np.nan
        self.assertTrue(np.isnan(fill_holes(a,self.c)[0][5:8,5:8]).all())

    def test_disabled_and_empty(self):
        self.c.FILL_SMALL_HOLES=False
        a=np.full((15,15),np.nan)
        result,flag=fill_holes(a,self.c)
        self.assertTrue(np.isnan(result).all())
        self.assertFalse(flag.any())

    def test_tile_halo_matches_full_array(self):
        a=np.arange(80*80,dtype="float32").reshape(80,80)/6400
        a[31:34,31:34]=np.nan  # straddles a tile boundary
        full,_=fill_holes(a,self.c)
        tiled=a.copy()
        halo=self.c.FILL_MAX_HOLE_PIXELS+self.c.FILL_RADIUS_PIXELS+1
        for row in range(0,80,32):
            for col in range(0,80,32):
                top,left=max(0,row-halo),max(0,col-halo)
                bottom,right=min(80,row+32+halo),min(80,col+32+halo)
                block,_=fill_holes(a[top:bottom,left:right],self.c)
                h,w=min(32,80-row),min(32,80-col)
                tiled[row:row+h,col:col+w]=block[row-top:row-top+h,col-left:col-left+w]
        np.testing.assert_array_equal(full,tiled)

if __name__=="__main__":
    unittest.main()
