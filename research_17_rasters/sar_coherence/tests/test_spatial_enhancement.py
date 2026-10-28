"""Small synthetic software checks; no research raster or simulation is opened."""
import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
import numpy as np
import pandas as pd
import config
from sarflow.recovery import Assessment, assess
from sarflow.temporal import Observation
from sarflow.spatial_enhancement import evidence_score, refine, SIGNATURES
from sarflow.output_layout import create_layout
try:
    from sarflow.spatial_objects import aggregate, cell_sizes, grow_pre_regions
    SPATIAL_AVAILABLE = True
except ImportError:
    SPATIAL_AVAILABLE = False


def settings(**changes):
    c = SimpleNamespace(**{k:v for k,v in vars(config).items() if not k.startswith('__')})
    for k,v in changes.items():
        setattr(c,k,v)
    return c


def result(shape=(5,7), code=1):
    metrics = {k:np.full(shape,.1) for k in SIGNATURES}
    metrics.update({k:np.zeros(shape) for k in ('level_se','progress_se','final_prediction_se')})
    metrics.update(pre_count=np.full(shape,8.),post_count=np.full(shape,8.),
                   adaptive_loss_threshold=np.full(shape,.03))
    return Assessment(np.full(shape,code,'uint8'),np.zeros(shape,'uint16'),
                      np.zeros(shape,'uint8'),metrics,None)


class EnhancementTests(unittest.TestCase):
    def setUp(self):
        self.c = settings(NEIGHBOR_RADIUS_M=40.)
        self.raw = result()
        self.original = np.full((5,7),5,'uint8')
        self.original[2,2] = 1
        self.score = np.zeros((5,7),'float32')
        self.score[2,2] = 1
        self.strata = np.ones((5,7),'uint32')

    def run_refine(self):
        return refine(self.original,self.score,self.raw,self.strata,self.c,(10.,10.))

    def test_score_range_and_uncertainty(self):
        self.assertTrue(np.all(evidence_score(self.raw,self.c)==1))
        self.raw.metrics['level_se'][:] = .10
        self.assertTrue(np.all(evidence_score(self.raw,self.c)<.8))

    def test_own_supported_candidate_can_change(self):
        out,score,reason = self.run_refine()
        self.assertEqual(out[2,3],1)
        self.assertGreaterEqual(score[2,3],.8)
        self.assertEqual(reason[2,3],1)

    def test_indeterminate_pixel_cannot_be_filled(self):
        self.raw.classes[:]=5
        out,_,_=self.run_refine()
        np.testing.assert_array_equal(out,self.original)

    def test_raw_no_loss_cannot_be_promoted(self):
        self.raw.classes[:]=4
        out,_,_=self.run_refine()
        np.testing.assert_array_equal(out,self.original)

    def test_high_confidence_original_protected(self):
        self.original[2,3]=3
        self.score[2,3]=.9
        self.assertEqual(self.run_refine()[0][2,3],3)

    def test_nodata_and_landcover_barrier(self):
        self.original[:,3]=0
        self.strata[:,3]=0
        self.assertTrue(np.all(self.run_refine()[0][:,4:]==5))
        self.original[:,3]=5
        self.strata[:,3:]=2
        self.assertTrue(np.all(self.run_refine()[0][:,3:]==5))

    def test_temporal_signature_disagreement_blocks_growth(self):
        self.raw.metrics['modeled_final'][:,3:]=.8
        self.assertTrue(np.all(self.run_refine()[0][:,3:]==5))

    def test_insufficient_own_support_blocks_growth(self):
        self.raw.metrics['post_count'][:]=5
        np.testing.assert_array_equal(self.run_refine()[0],self.original)

    def test_halo_reproduces_interior(self):
        # Candidate distance is bounded: enough identical context yields identical core.
        raw=result((11,17)); original=np.full((11,17),5,'uint8');original[5,7]=1
        score=np.zeros((11,17),'float32');score[5,7]=1
        full=refine(original,score,raw,np.ones_like(original),self.c,(10.,10.))[0]
        cropped=Assessment(raw.classes[:,3:14],raw.quality[:,3:14],raw.reasons[:,3:14],
                           {k:v[:,3:14] for k,v in raw.metrics.items()},None)
        sub=refine(original[:,3:14],score[:,3:14],cropped,np.ones((11,11),'uint8'),self.c,(10.,10.))[0]
        np.testing.assert_array_equal(full[:,7:10],sub[:,4:7])

    @unittest.skipUnless(SPATIAL_AVAILABLE, 'SciPy runtime unavailable')
    def test_coverage_threshold_no_zero_fill(self):
        stack=np.full((17,2,5),.5)
        stack[0,0,:3]=np.nan
        series,counts,fractions,_=aggregate(stack,np.ones((2,5),bool),self.c)
        self.assertTrue(np.isnan(series[0]))
        self.assertEqual(fractions[0],.7)

    @unittest.skipUnless(SPATIAL_AVAILABLE, 'SciPy runtime unavailable')
    def test_physical_area_limits(self):
        cell,domain=cell_sizes(200.,self.c)
        self.assertLessEqual(domain*domain*200.,90000.)
        self.assertLessEqual(cell*cell*200.,30000.)

    @unittest.skipUnless(SPATIAL_AVAILABLE, 'SciPy runtime unavailable')
    def test_pre_only_growth_preserves_barrier(self):
        ids=np.array([[1,1,3,3],[1,1,3,3]],dtype='uint32')
        features=np.zeros((3,2,4))
        barriers=np.ones((2,4),'uint32');barriers[:,2:]=2
        np.testing.assert_array_equal(grow_pre_regions(ids,features,barriers,20,self.c),ids)
        merged=grow_pre_regions(ids,features,np.ones_like(barriers),20,self.c)
        self.assertEqual(len(np.unique(merged)),1)

    def test_output_layout_no_overwrite(self):
        # Temporary fixtures remain under the designated project, never source inputs.
        with TemporaryDirectory(dir=config.PROJECT_ROOT) as temp:
            c=settings(PROJECT_ROOT=Path(temp),Output_Folder=temp)
            layout=create_layout(c,'test_run')
            self.assertEqual(layout.graphs.parent.name,'ITSA Graphs')
            (layout.tables/'existing.txt').write_text('fixture')
            with self.assertRaises(FileExistsError):
                create_layout(c,'test_run')

    def test_actual_17_date_model_classes(self):
        origin = pd.Timestamp('2022-09-25')
        days = [-90,-78,-66,-54,-42,-30,-18,-6,6,18,30,42,54,66,78,90,102]
        phases = ['pre']*8+['event']+['post']*8
        obs = []
        for i,(day,phase) in enumerate(zip(days,phases)):
            end=origin+pd.Timedelta(days=day);start=end-pd.Timedelta(days=12)
            obs.append(Observation(str(i),Path('synthetic.tif'),phase,start,end,start,end,end))
        c=settings(EVENT_ORIGIN=str(origin),SOURCE_VALIDITY_CONFIRMED=True)
        for slope,loss,expected in ((.003,.2,1),(.0008,.2,2),(0,.2,3),(0,0,4)):
            y=np.full(17,.6);y[8]=.6-loss;y[9:]=.6-loss+slope*np.array(days[9:])
            self.assertEqual(int(assess(y[:,None],obs,c).classes[0]),expected)
        y[:5]=np.nan
        self.assertEqual(int(assess(y[:,None],obs,c).classes[0]),5)


@unittest.skipUnless(all(importlib.util.find_spec(n) for n in ('rasterio','shapely','pyproj','defusedxml')),
                     'Geospatial runtime dependencies unavailable')
class GeographicTests(unittest.TestCase):
    def test_actual_barangay_names_and_geometry(self):
        from sarflow.analysis_boundaries import load_barangays, safe_name
        records,report=load_barangays(config.BARANGAY_BOUNDARY,'EPSG:32651',settings())
        self.assertEqual(report['valid_features'],36)
        self.assertEqual(len({b.code for b in records}),36)
        self.assertIn('Abiawin',[b.name for b in records])
        self.assertNotIn('/',safe_name('A/B'))

    def test_owner_no_double_count(self):
        from rasterio.transform import from_origin
        from sarflow.analysis_boundaries import Barangay, owners
        geom={'type':'Polygon','coordinates':[[[0,0],[2,0],[2,2],[0,2],[0,0]]]}
        records=[Barangay('A','1','A',geom),Barangay('B','2','B',geom)]
        self.assertTrue(np.all(owners(records,(2,2),from_origin(0,2,1,1))==1))

    def test_exact_disk_quantiles(self):
        from sarflow.barangay_outputs import DiskValues
        with TemporaryDirectory(dir=config.PROJECT_ROOT) as temp:
            values=DiskValues(temp)
            values.add('date',np.array([.2,np.nan,.4,.6,.8]))
            stats=values.stats('date')
            self.assertEqual(stats['count'],4)
            self.assertAlmostEqual(stats['median'],.5,places=6)
            self.assertEqual(values.stats('missing')['count'],0)

    def test_global_components_cross_row_and_strip(self):
        import rasterio
        from rasterio.transform import from_origin
        from sarflow.enhancement_reports import connected_region_report
        with TemporaryDirectory(dir=config.PROJECT_ROOT) as temp:
            folder=Path(temp); path=folder/'tiny.tif'
            data=np.ones((130,3),dtype='uint8');data[:,2]=3
            with rasterio.open(path,'w',driver='GTiff',width=3,height=130,count=1,dtype='uint8',
                               crs='EPSG:32651',transform=from_origin(0,1300,10,10),nodata=0) as ds:
                ds.write(data,1)
            connected_region_report(path,folder,settings(),100.)
            table=pd.read_csv(folder/'Spatial_Region_Statistics.csv')
            self.assertEqual(len(table),2)
            self.assertEqual(table.pixels.sum(),390)


if __name__=='__main__':
    unittest.main()
