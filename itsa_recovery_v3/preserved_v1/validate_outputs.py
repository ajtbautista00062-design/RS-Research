"""Independent window-by-window arithmetic and GeoTIFF validation on actual inputs."""
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT

import config as c
from sarflow.spatial import windows, read_zones, mask_for
from sarflow.temporal import validate


def main():
    observations = validate(c)
    out = Path(c.Output_Folder)
    output = out / 'Infanta_Coherence_Recovery_Fraction_and_Class.tif'
    counts = {str(code): 0 for code in (-9999, 1, 2, 3, 4)}
    max_error, checked, warps = 0., 0, []
    with ExitStack() as context:
        srcs = [context.enter_context(rasterio.open(o.path)) for o in observations]
        ref = srcs[0]
        dst = context.enter_context(rasterio.open(output))
        assert dst.driver == 'GTiff' and dst.count == 2
        assert dst.crs == ref.crs and dst.transform == ref.transform and dst.shape == ref.shape
        assert dst.dtypes == ('float32', 'float32') and dst.nodata == -9999
        assert dst.is_tiled and dst.compression.name == 'deflate'
        for source in srcs:
            warps.append(context.enter_context(WarpedVRT(
                source, crs=ref.crs, transform=ref.transform, width=ref.width,
                height=ref.height, dtype='float32', nodata=np.nan,
                resampling=Resampling[c.RESAMPLING],
                **({'src_nodata': 0} if c.ZERO_IS_INVALID and source.nodata is None else {}))))
        zones = read_zones(c.Study_Area, ref.crs, c.PROJECT_ROOT) if c.Study_Area else []
        for win in windows(ref.width, ref.height, c.WINDOW_SIZE):
            data = np.stack([(v.read(c.BAND, window=win, masked=True).filled(np.nan)
                              * s.scales[c.BAND-1] + s.offsets[c.BAND-1])
                             * c.COHERENCE_SCALE + c.COHERENCE_OFFSET
                             for v, s in zip(warps, srcs)])
            data = np.clip(data, 0, 1)
            if c.ZERO_IS_INVALID:
                data[data == 0] = np.nan
            mask = mask_for([g for _, g in zones], data.shape[1:],
                            rasterio.windows.transform(win, ref.transform), c.ALL_TOUCHED)
            pre, event, post = data
            loss = pre-event
            valid = np.isfinite(data).all(axis=0) & (pre >= c.BASELINE_MIN) & mask
            affected = valid & (loss >= max(c.LOSS_ABSOLUTE, c.LOSS_MAD_MULTIPLIER*c.NOISE_FLOOR))
            expected = np.full(pre.shape, -9999, dtype='float32')
            expected[valid & ~affected] = 4
            returned = affected & (post >= np.maximum(0, pre-c.RECOVERY_TOLERANCE))
            expected[returned] = 1
            expected[affected & ~returned & (post > event)] = 2
            expected[affected & ~returned & (post <= event)] = 3
            actual = dst.read(2, window=win)
            np.testing.assert_array_equal(actual, expected)
            fraction = dst.read(1, window=win)
            np.testing.assert_array_equal(fraction == -9999, ~affected)
            if affected.any():
                error = np.abs(fraction[affected].astype('float64') -
                               (post[affected].astype('float64')-event[affected]) /
                               (pre[affected].astype('float64')-event[affected]))
                max_error = max(max_error, float(error.max()))
                assert error.max() < 2e-6
            for code in counts:
                counts[code] += int((actual == int(code)).sum())
            checked += actual.size
        assert checked == ref.width * ref.height
        metadata = {'shape': dst.shape, 'crs': str(dst.crs), 'transform': tuple(dst.transform),
                    'bands': dst.descriptions, 'nodata': dst.nodata, 'compression': dst.compression.name}
    original = json.loads((c.PROJECT_ROOT / 'original_source_hashes.json').read_text(encoding='utf-8-sig'))
    unchanged, differences = 0, []
    for entry in original:
        path = c.PROJECT_ROOT.parent / entry['path']
        if not path.is_file():
            differences.append({'path': entry['path'], 'status': 'missing'})
        elif hashlib.sha256(path.read_bytes()).hexdigest().upper() != entry['sha256']:
            differences.append({'path': entry['path'], 'status': 'changed_since_initial_snapshot'})
        else:
            unchanged += 1
    assert pd.read_csv(out / 'observation_dates.csv').shape[0] == 3
    assert pd.read_csv(out / 'timeseries_summary.csv').shape[0] == 3
    for number in ('fig_coherence_timeseries_20220925', '02_coherence_anomaly', '03_phase_distributions'):
        for ext in c.FIGURE_FORMATS:
            assert (out / f'{number}.{ext}').stat().st_size > 0
    assert (out / 'interactive_timeseries.html').stat().st_size > 0
    report = {'status': 'PASS', 'pixels_independently_checked': checked,
              'class_counts': counts, 'maximum_fraction_absolute_error_float64_reference': max_error,
              'geotiff': metadata, 'original_files_sha256_unchanged': unchanged,
              'original_snapshot_status': 'DIFFERS' if differences else 'UNCHANGED',
              'original_snapshot_differences': differences,
              'qgis_validation': 'GDAL/Rasterio grid, CRS, band, NoData and class checks passed; QGIS GUI not tested'}
    (out / 'validation_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
