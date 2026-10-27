"""Fixed-size publication files plus the exact statistical inputs used to draw them."""
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from .style import DESIGN
from .captions import make_caption


def export(fig, out, stem, data, kind, number, font):
    boundary = Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                         facecolor='none', edgecolor='none', linewidth=0)
    fig.add_artist(boundary)
    metadata = {
        'figure_inches': fig.get_size_inches().tolist(), 'font': font, 'dpi': DESIGN['dpi'],
        'axis_limits': [{'x': list(ax.get_xlim()), 'y': list(ax.get_ylim())} for ax in fig.axes],
        'observation_dates': [o.end.isoformat() for o in data.obs],
        'pair_intervals': [[o.start.isoformat(), o.end.isoformat()] for o in data.obs],
        'event_date': DESIGN['event_date'],
        'statistics': data.table.to_dict(orient='records'),
        'anomaly_quantiles': np.quantile(data.anomaly, [.25, .5, .75], axis=1).tolist(),
        'response_quantiles': (np.quantile(data.normalized, [.25, .5, .75], axis=1).tolist()
                               if data.eligible.any() else []),
        'phase_quantiles': [np.quantile(group, [.25, .5, .75]).tolist() for group in data.groups],
        'eligible_sample_pixels': int(data.eligible.sum()),
        'titles_empty': all(not ax.get_title() for ax in fig.axes),
        'synthetic_component_test': getattr(data.c, 'SYNTHETIC_VISUALIZATION_TEST', False),
        'zone_support': {name:{'common_pixels':sample.count,'sample_pixels':sample.values.shape[1],
                                'mean':(sample.sums/sample.count).tolist()}
                         for name,sample in data.zones.items() if sample.count},
    }
    out.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(out/f'{stem}.{ext}', dpi=DESIGN['dpi'], bbox_inches='tight', pad_inches=0,
                    bbox_extra_artists=[boundary])
    (out/f'{stem}_caption.md').write_text(make_caption(data, kind, number, font, metadata), encoding='utf-8')
    (out/f'{stem}_data.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    plt.close(fig)
    return {'stem': stem, 'kind': kind, 'status': 'generated'}
