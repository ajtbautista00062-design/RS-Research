"""End-to-end run on small real GeoTIFFs: 8 PRE + 1 event + 8 POST through config.PATHS."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from pyproj import Transformer
import config
from sarflow.study_pipeline import configuration, run
from sarflow.study_inputs import SLOTS, observations
from sarflow.study_layout import folders


def pair_names():
    """Consecutive 12-day pairs: 8 before, 1 straddling 2022-09-25, 8 after."""
    first = pd.Timestamp('2022-06-15')
    dates = [first + pd.Timedelta(days=12 * i) for i in range(18)]
    return [f'coh_{a:%Y%m%d}_{b:%Y%m%d}.tif' for a, b in zip(dates, dates[1:])]


def write_inputs(base, values):
    x, y = Transformer.from_crs(4326, 3123, always_xy=True).transform(121.6, 14.7)
    affine = from_origin(x, y + 40, 10, 10)
    paths = {}
    for slot, name, value in zip(SLOTS, pair_names(), values):
        path = base / 'Coherence' / name
        with rasterio.open(path, 'w', driver='GTiff', height=4, width=4, count=1, dtype='float32',
                           crs='EPSG:3123', transform=affine, nodata=-9999) as ds:
            a = np.full((4, 4), value, dtype='float32')
            a[0, 0] = -9999
            ds.write(a, 1)
        paths[slot] = str(path)
    reverse = Transformer.from_crs(3123, 4326, always_xy=True)
    pieces = []
    for index, (left, right) in enumerate([(x, x + 20), (x + 20, x + 40)], 1):
        coords = ' '.join(','.join(map(str, reverse.transform(a, b))) for a, b in
                          [(left, y), (right, y), (right, y + 40), (left, y + 40), (left, y)])
        pieces.append(f'<Placemark><ExtendedData><SchemaData><SimpleData name="ADM4_EN">Zone {index}</SimpleData>'
                      f'<SimpleData name="ADM4_PCODE">Z{index}</SimpleData></SchemaData></ExtendedData>'
                      f'<Polygon><outerBoundaryIs><LinearRing><coordinates>{coords}</coordinates></LinearRing>'
                      '</outerBoundaryIs></Polygon></Placemark>')
    boundary = base / 'Boundary' / 'boundary.kml'
    boundary.write_text('<kml xmlns="http://www.opengis.net/kml/2.2"><Document>' + ''.join(pieces) + '</Document></kml>')
    return paths, boundary


def settings(base, paths, boundary, **changes):
    values = dict(PROJECT_ROOT=base, RESULTS_ROOT=base / 'Results', VERSION_ROOT=base / 'Results',
                  PATHS=paths, PAIR_DATES={s: ('', '') for s in SLOTS},
                  Study_Area=str(boundary), BARANGAY_BOUNDARY=str(boundary), FIGURE_FORMATS=('png',),
                  ACQUISITION_COMPARABILITY_CONFIRMED=True, COMPARABILITY_NOTES='synthetic fixture',
                  ZONAL_MIN_VALID_PIXELS=1, ZONAL_BASE_CELLS=2, ZONAL_LEVELS=2,
                  ZONAL_SENSITIVITY_BASE_CELLS=(1,))
    values.update(changes)
    return configuration(config, **values)


class ExportIntegrationTests(unittest.TestCase):
    def test_seventeen_inputs_run_end_to_end(self):
        with TemporaryDirectory(dir=config.PROJECT_ROOT) as tmp:
            base = Path(tmp)
            (base / 'Coherence').mkdir()
            (base / 'Boundary').mkdir()
            values = [.3] * 8 + [.1] + [.12, .16, .20, .24, .27, .29, .30, .30]
            paths, boundary = write_inputs(base, values)
            c = settings(base, paths, boundary)
            obs = observations(c)
            self.assertEqual([o.phase for o in obs].count('pre'), 8)
            self.assertEqual([o.phase for o in obs].count('post'), 8)
            destination = run(c, 'synthetic')
            state = json.loads((destination / 'run_state.json').read_text())
            self.assertEqual(state['status'], 'complete')
            self.assertEqual(state['verification']['study_pixels'], 16)
            targets = pd.read_csv(destination / 'Validation' / 'Target_Assessment.csv')
            self.assertTrue((targets.population == targets[[f'corrected_class_{k}' for k in range(6)]].sum(axis=1)).all())
            with rasterio.open(folders(destination, 'ITSA').rasters / 'ITSA_Recovery_Classes.tif') as ds:
                self.assertEqual(int((ds.read(1) == 0).sum()), 1)
            for zone in ('Zone 1_Z1', 'Zone 2_Z2'):
                self.assertTrue((folders(destination, 'ITSA', zone).graphs / 'ITSA_Time_Series.png').is_file())
            # Zonal products use the standard A. RRI / B. ITSA layout, never a separate folder.
            self.assertEqual(sorted(p.name for p in destination.iterdir() if p.is_dir()),
                             ['A. RRI', 'B. ITSA', 'Validation'])
            with rasterio.open(folders(destination, 'ITSA').rasters / 'Zonal_Recovery_Classes.tif') as ds:
                classes = ds.read(1)
                self.assertEqual(ds.colormap(1)[1][:3], (26, 150, 65))
            self.assertEqual(int((classes == 1).sum()), 15)  # loss then full return; 1 NoData pixel
            for zone in (None, 'Zone 1_Z1', 'Zone 2_Z2'):
                itsa, rri = folders(destination, 'ITSA', zone), folders(destination, 'RRI', zone)
                for path in (itsa.rasters / 'Zonal_Recovery_Classes.tif', itsa.rasters / 'Zonal_Recovery_Classes.qml',
                             itsa.graphs / 'Zonal_Recovery_Classes_Map.png', itsa.csv / 'Zonal_Class_Area_Statistics.csv',
                             rri.rasters / 'Zonal_RRI_Interpretation_Classes.tif', rri.graphs / 'Zonal_RRI_Interpretation_Map.png',
                             rri.csv / 'Zonal_RRI_Class_Statistics.csv'):
                    self.assertTrue(path.is_file(), path)
            self.assertTrue((folders(destination, 'ITSA').rasters / 'Zonal_Recovery_Zones.geojson').is_file())
            from sarflow.zonal_recovery import run_zonal
            with self.assertRaises(FileExistsError):
                run_zonal(destination, c)
            run_zonal(destination, c, replace=True)
            manifest = json.loads((destination / 'run_manifest.json').read_text())
            self.assertIn('zonal_recovery', manifest)

    def test_empty_slots_skipped_and_minimums_enforced(self):
        with TemporaryDirectory(dir=config.PROJECT_ROOT) as tmp:
            base = Path(tmp)
            (base / 'Coherence').mkdir()
            (base / 'Boundary').mkdir()
            paths, boundary = write_inputs(base, [.3] * 17)
            for slot in ('Pre1', 'Pre2', 'Post7', 'Post8'):
                paths[slot] = ''
            self.assertEqual(len(observations(settings(base, paths, boundary))), 13)
            relative = {k: (str(Path(v).relative_to(base)) if v else '') for k, v in paths.items()}
            self.assertEqual(len(observations(settings(base, relative, boundary))), 13)
            with TemporaryDirectory() as outside:  # rasters on another folder/drive are accepted
                moved = {}
                for k, v in paths.items():
                    if v:
                        target = Path(outside) / Path(v).name
                        target.write_bytes(Path(v).read_bytes())
                        moved[k] = str(target)
                    else:
                        moved[k] = ''
                self.assertEqual(len(observations(settings(base, moved, boundary))), 13)
            for slot in ('Pre3', 'Pre4', 'Pre5', 'Pre6'):
                paths[slot] = ''
            with self.assertRaisesRegex(ValueError, 'at least 3 PRE'):
                observations(settings(base, paths, boundary))
            paths, boundary = write_inputs(base, [.3] * 17)
            with self.assertRaisesRegex(ValueError, 'ACQUISITION_COMPARABILITY_CONFIRMED'):
                observations(settings(base, paths, boundary, ACQUISITION_COMPARABILITY_CONFIRMED=False))
            paths['During'] = ''
            with self.assertRaisesRegex(ValueError, 'event-straddling raster is required'):
                observations(settings(base, paths, boundary))


if __name__ == '__main__':
    unittest.main()
