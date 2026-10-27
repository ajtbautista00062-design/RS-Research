"""Synthetic specification tests. Written for V4; NOT executed during implementation."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
sys.path.insert(1,str(Path(__file__).resolve().parents[1]/'preserved_v1'))
import unittest
import numpy as np
from scipy.ndimage import label as connected
import config
from object_segmentation import pre_features,segment_cell,stable_ids
from spatial_io import active_cells


class FakeBoundary:
    def __init__(self,a):
        self.a=a; self.height,self.width=a.shape
    def read(self,band,window):
        r,c=int(window.row_off),int(window.col_off)
        return self.a[r:r+int(window.height),c:c+int(window.width)]


class ObjectSegmentationTests(unittest.TestCase):
    def segment(self,features,valid=None,strata=None,target=16):
        shape=features.shape[1:]
        valid=np.ones(shape,bool) if valid is None else valid
        strata=np.ones(shape,'uint32') if strata is None else strata
        return segment_cell(features,valid,strata,target,64,4,np.array([.12,.06,.06]))

    def test_connected_compact_regions_and_maximum(self):
        features=np.zeros((3,8,8)); features[0]=.7
        ids=self.segment(features)
        self.assertTrue(np.all(ids>0))
        for key in np.unique(ids):
            self.assertEqual(connected(ids==key)[1],1)
            self.assertLessEqual(int((ids==key).sum()),64)

    def test_heterogeneous_pre_regions_remain_separate(self):
        features=np.zeros((3,8,8)); features[0,:,:4]=.2; features[0,:,4:]=.8
        ids=self.segment(features)
        self.assertFalse(set(ids[:,:4].ravel())&set(ids[:,4:].ravel()))

    def test_nodata_and_hard_boundaries_are_not_crossed(self):
        f=np.zeros((3,8,8)); strata=np.ones((8,8),'uint32'); strata[:,5:]=2
        valid=np.ones((8,8),bool); valid[:,3]=False
        ids=self.segment(f,valid,strata)
        self.assertTrue(np.all(ids[:,3]==0))
        self.assertFalse(set(ids[:,4][ids[:,4]>0])&set(ids[:,5][ids[:,5]>0]))

    def test_pre_features_do_not_mix_independent_categories(self):
        pre=np.full((3,5,6),.2); pre[:,:,3:]=.8
        strata=np.ones((5,6),'uint32'); strata[:,3:]=2
        features=pre_features(pre,strata,1)
        np.testing.assert_allclose(features[0,:,:3],.2)
        np.testing.assert_allclose(features[0,:,3:],.8)

    def test_minimum_mapping_unit_merge_requires_pre_similarity(self):
        f=np.zeros((3,3,3)); valid=np.ones((3,3),bool)
        ids=segment_cell(f,valid,np.ones((3,3)),1,9,4,np.ones(3)*.1)
        self.assertLessEqual(np.unique(ids).size,2)
        self.assertTrue(all((ids==key).sum()>=4 for key in np.unique(ids)))
        f[0,0,0]=.9
        ids=segment_cell(f,valid,np.ones((3,3)),1,9,4,np.ones(3)*.1)
        self.assertEqual(int((ids==ids[0,0]).sum()),1)

    def test_io_tile_boundaries_do_not_change_global_region_ids(self):
        source=FakeBoundary(np.zeros((16,16),dtype='uint8'))
        def build(io_size):
            out=np.zeros((16,16),'uint32')
            for window in active_cells(source,8,io_size):
                f=np.zeros((3,8,8))
                local=self.segment(f,target=64)
                ids=stable_ids(local,int(window.row_off),int(window.col_off),16)
                r,c=int(window.row_off),int(window.col_off)
                out[r:r+8,c:c+8]=ids
            return out
        np.testing.assert_array_equal(build(3),build(11))
        self.assertEqual(build(3)[1,2],build(3)[1,4])

    def test_ids_are_globally_unique(self):
        labels=np.ones((2,2),'uint32')
        self.assertNotEqual(stable_ids(labels,0,0,20)[0,0],stable_ids(labels,0,2,20)[0,0])
