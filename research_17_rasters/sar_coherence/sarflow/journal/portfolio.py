"""Adaptive 3/17-observation journal figure suite; no raster processing."""
import json
import numpy as np
import matplotlib as mpl
import pandas as pd
from .style import apply_publication_style, DESIGN
from .data import prepare
from .export import export
from .timeseries import main_figure
from .panels import (diagnostic_figure, anomaly_panel, distribution_panel, recovery_panel,
                     fraction_panel, overview_figure)
from .zones import heatmap, zone_series
from .interactive import interactive


def render_portfolio(sample, zones, obs, out, c, fraction, metrics_module,
                     only_main=False, single_column=False):
    if not sample.count:
        raise ValueError('No common valid pixels: publication statistics unavailable')
    data = prepare(sample, zones, obs, c, fraction, metrics_module)
    generated = []
    with mpl.rc_context():
        font = apply_publication_style()
        def save(figure, stem, kind):
            generated.append(export(figure, out, stem, data, kind, len(generated)+1, font))
        save(main_figure(data), DESIGN['stem'], 'main')
        if single_column:
            save(main_figure(data, single=True), DESIGN['stem']+'_single_column', 'main')
        if not only_main:
            save(diagnostic_figure(data, anomaly_panel), '02_coherence_anomaly', 'anomaly')
            save(diagnostic_figure(data, distribution_panel), '03_phase_distributions', 'distribution')
            if data.eligible.any():
                stem = '04_single_post_response' if len(data.post_indices) == 1 else '04_recovery_trajectory'
                save(diagnostic_figure(data, recovery_panel), stem, 'recovery')
            else:
                generated.append({'kind':'recovery','status':'unavailable','reason':'No eligible sampled pixels'})
            if np.isfinite(fraction).any():
                save(diagnostic_figure(data, fraction_panel), '06_recovery_fraction_distribution', 'fraction')
            else:
                generated.append({'kind':'fraction','status':'unavailable','reason':'No finite supplied return fractions'})
            real_zones = {name:s for name,s in zones.items() if s.count}
            if real_zones:
                names = list(real_zones)
                for start in range(0, len(names), 12):
                    selected = {name:real_zones[name] for name in names[start:start+12]}
                    suffix = '' if len(names) <= 12 else f'_page_{start//12+1:02d}'
                    save(heatmap(data, selected), '05_zone_heatmap'+suffix, 'heatmap')
                for start in range(0, len(names), 6):
                    selected = {name:real_zones[name] for name in names[start:start+6]}
                    suffix = '' if len(names) <= 6 else f'_page_{start//6+1:02d}'
                    save(zone_series(data, selected), '07_zone_trajectories'+suffix, 'zones')
            else:
                generated.append({'kind':'heatmap/zones','status':'unavailable',
                                  'reason':'No real nonempty polygon-zone samples supplied'})
            save(overview_figure(data), 'fig_coherence_synthesis_20220925', 'overview')
            interactive(sample, obs, out, c)
            # Preserve the original descriptive trend table for the actual 8-post workflow.
            # No trend line is added to the graphics, and no slope is estimated for one post.
            if len(data.post_indices) == 8 and data.eligible.any():
                from scipy.stats import theilslopes
                import pandas as pd
                q = np.quantile(data.normalized, [.25,.5,.75], axis=1)
                days = np.array([(obs[i].reference-pd.Timestamp(c.EVENT_ORIGIN)).total_seconds()/86400
                                 for i in data.post_indices])
                trend = theilslopes(q[1], days)
                pd.DataFrame([{'series':'spatial_median_normalized_recovery',
                               'theil_sen_slope_per_day':trend.slope,
                               'interpretation':'descriptive; no independent-observation p-value or CI claimed'}]
                             ).to_csv(out/'recovery_trend.csv', index=False)
    manifest = {'observation_count':len(obs),'common_pixels':sample.count,
                'sample_pixels':sample.values.shape[1],'figures':generated,
                'single_post_limit':len(data.post_indices)==1,
                'synthetic_component_test':getattr(c,'SYNTHETIC_VISUALIZATION_TEST',False)}
    filename = 'publication_main_manifest.json' if only_main else 'publication_manifest.json'
    (out/filename).write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return manifest
