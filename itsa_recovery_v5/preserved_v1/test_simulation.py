"""Focused numerical, alignment, masking, validation, and export regressions."""
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np
import rasterio
from rasterio.transform import from_origin

import config
from sarflow.metrics import calculate_fraction
from sarflow.sampling import Sample
from sarflow.spatial import aligned, read_stack, windows, verify_source_values, read_zones, mask_for
from sarflow.temporal import validate
from sarflow.zonal import ZoneStats


class SimulationTests(unittest.TestCase):
    def setUp(self):
        self.obs = validate(config)

    def test_classes_and_unclipped_fraction(self):
        data = np.array([[.8]*7, [.4]*5 + [.75, .75],
                         [.8, .6, .4, .2, 1., .75, np.nan]], dtype='float32')
        result = calculate_fraction(data, self.obs, config)
        np.testing.assert_array_equal(result.classes, [1, 2, 3, 3, 1, 4, -9999])
        np.testing.assert_allclose(result.fraction[:5], [1, .5, 0, -.5, 1.5], atol=1e-6)
        self.assertTrue(np.isnan(result.fraction[5:]).all())
        self.assertTrue(np.isnan(result.spread).all())

    def test_invalid_and_low_pre_coherence(self):
        data = np.array([[np.nan, .8, .8, .3], [.4, np.nan, .4, .1],
                         [.8, .8, np.nan, .3]], dtype='float32')
        result = calculate_fraction(data, self.obs, config)
        self.assertTrue((result.classes == -9999).all())

    def test_threshold_boundaries(self):
        # Exactly representable inputs and config avoid ambiguity at boundaries.
        c = SimpleNamespace(**{key: getattr(config, key) for key in dir(config) if key.isupper()})
        c.LOSS_ABSOLUTE, c.RECOVERY_TOLERANCE = .125, .0625
        data = np.array([[.75, .75, .75], [.625, .625, .6875], [.6875, .6874, .2]])
        result = calculate_fraction(data, self.obs, c)
        np.testing.assert_array_equal(result.classes, [1, 2, 4])

    def test_reject_persistence(self):
        with self.assertRaisesRegex(ValueError, 'Persistent'):
            calculate_fraction(np.ones((3, 2)), self.obs, config, confirmations=2)

    def test_interactive_html_without_plotly_theme_loading(self):
        from unittest.mock import patch
        from sarflow.figures import interactive
        import plotly.graph_objects as go
        sample = Sample(2, 1)
        sample.add(np.array([[.8, .7], [.4, .3], [.75, .65]], dtype='float32'),
                   np.arange(2), np.arange(2))
        with TemporaryDirectory(dir=config.PROJECT_ROOT) as folder:
            # Export must not instantiate graph objects or copy a bundled theme.
            with patch.object(go, 'Figure', side_effect=AssertionError('Unexpected graph-object loading')):
                interactive(sample, self.obs, Path(folder), config)
            html = (Path(folder) / 'interactive_timeseries.html').read_text(encoding='utf-8')
            self.assertIn('Plotly.newPlot', html)
            self.assertIn('no recovery trend', html)
            self.assertIn('2022-10-13T00:00:00', html)

    def test_untagged_zero_footprint_is_masked_before_warp(self):
        with TemporaryDirectory(dir=config.PROJECT_ROOT) as folder:
            path = Path(folder) / 'zero_footprint.tif'
            with rasterio.open(path, 'w', driver='GTiff', width=4, height=4, count=1,
                               crs='EPSG:32651', transform=from_origin(500000, 1600000, 10, 10),
                               dtype='float32') as dst:
                data = np.full((4, 4), .8, dtype='float32')
                data[:, :2] = 0
                dst.write(data, 1)
            with rasterio.open(path) as src, aligned(src, src, config) as vrt:
                stack = read_stack([vrt], [src], next(windows(4, 4, 256)), config)
                self.assertTrue(np.isnan(stack[0, :, :2]).all())
                np.testing.assert_allclose(stack[0, :, 2:], .8)

    def test_temporal_order_and_comparability_gate(self):
        self.assertEqual([(o.end-o.start).days for o in self.obs], [12, 12, 12])
        self.assertEqual([o.name for o in self.obs], ['Pre1', 'During', 'Post1'])
        values = {key: getattr(config, key) for key in dir(config) if not key.startswith('__')}
        c = SimpleNamespace(**values)
        c.PAIR_DATES = dict(config.PAIR_DATES, Post1=('2022-09-20', '2022-10-02'))
        with self.assertRaisesRegex(ValueError, 'post pair'):
            validate(c)
        c.PAIR_DATES = config.PAIR_DATES
        c.ALLOW_UNVERIFIED_COMPARABILITY_FOR_SIMULATION = False
        with self.assertRaisesRegex(ValueError, 'comparability'):
            validate(c)

    def test_sampling_and_zonal_counts(self):
        data = np.array([[.8]*4, [.4]*3+[.75], [.8, .6, .4, .75]], dtype='float32')
        a, b = Sample(2, 1), Sample(2, 1)
        for sample in (a, b):
            sample.add(data, np.arange(4), np.arange(4))
        np.testing.assert_array_equal(a.values, b.values)
        np.testing.assert_allclose(a.table(list(config.PATHS)).mean_exact, data.mean(axis=1))
        result = calculate_fraction(data, self.obs, config)
        zone = ZoneStats('test')
        zone.add(result.fraction, result.classes, np.ones(4, dtype=bool))
        self.assertEqual(zone.counts.tolist(), [0, 1, 1, 1, 1])
        self.assertEqual(zone.summary()['recovered_fraction_among_eligible'], 1/3)

    def test_alignment_source_range_and_kml(self):
        with TemporaryDirectory(dir=config.PROJECT_ROOT) as folder:
            root = Path(folder)
            profile = dict(driver='GTiff', width=4, height=4, count=1,
                           crs='EPSG:32651', transform=from_origin(500000, 1600000, 10, 10),
                           dtype='float32', nodata=-9999)
            paths = []
            for i, value in enumerate((.8, .4, .75)):
                path = root / f'{i}.tif'
                p = dict(profile)
                if i == 2:
                    p.update(width=2, height=2, transform=from_origin(500000, 1600000, 20, 20))
                with rasterio.open(path, 'w', **p) as dst:
                    dst.write(np.full((p['height'], p['width']), value, dtype='float32'), 1)
                paths.append(path)
            from contextlib import ExitStack
            with ExitStack() as context:
                srcs = [context.enter_context(rasterio.open(p)) for p in paths]
                vrts = [context.enter_context(aligned(src, srcs[0], config)) for src in srcs]
                stack = read_stack(vrts, srcs, next(windows(4, 4, 256)), config)
                self.assertEqual(stack.shape, (3, 4, 4))
                np.testing.assert_allclose(stack[2], .75)
            with rasterio.open(paths[0], 'r+') as dst:
                dst.write(np.full((4, 4), 1.5, dtype='float32'), 1)
            with rasterio.open(paths[0]) as src:
                with self.assertRaisesRegex(ValueError, 'raw coherence range'):
                    verify_source_values(src, config)
            kml = root / 'test.kml'
            kml.write_text('<kml><Placemark><name>test</name><Polygon><outerBoundaryIs>'
                           '<LinearRing><coordinates>121,14 121.1,14 121.1,14.1 121,14.1 121,14</coordinates>'
                           '</LinearRing></outerBoundaryIs></Polygon></Placemark></kml>')
            zones = read_zones(str(kml), rasterio.crs.CRS.from_epsg(32651), config.PROJECT_ROOT)
            self.assertEqual(len(zones), 1)
            self.assertEqual(mask_for([zones[0][1]], (4, 4), profile['transform'], False).shape, (4, 4))


if __name__ == '__main__':
    unittest.main(verbosity=2)
