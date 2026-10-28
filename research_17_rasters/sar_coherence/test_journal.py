"""Exact numeric and layout checks; --synthetic17 is component validation only."""
import argparse
import csv
import json
import os
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import config as c
from sarflow import metrics
from sarflow.sampling import Sample
from sarflow.temporal import validate, Observation
from sarflow.journal.data import prepare
from sarflow.journal.cache import read_cache, write_cache
from sarflow.journal.style import apply_publication_style
from sarflow.journal.timeseries import main_figure
from sarflow.journal.panels import (overview_figure, diagnostic_figure, anomaly_panel,
                                    distribution_panel, recovery_panel, fraction_panel)
from sarflow.journal.zones import heatmap, zone_series
from sarflow.journal.portfolio import render_portfolio


def synthetic_fixture():
    if len(c.PATHS)!=17:
        raise ValueError('Run the synthetic17 check from the actual Codes package')
    config = SimpleNamespace(**{k:getattr(c,k) for k in dir(c) if not k.startswith('__')})
    config.SYNTHETIC_VISUALIZATION_TEST = True
    config.INTERACTIVE_HTML = False
    origin = pd.Timestamp('2022-09-19')
    names = list(c.PATHS)
    obs = []
    for i,name in enumerate(names):
        phase = 'pre' if i<8 else 'event' if i==8 else 'post'
        offset = (i-8)*12 + (max(0,i-9)//2)*3
        start = origin+pd.Timedelta(days=offset)
        end = start+pd.Timedelta(days=12)
        obs.append(Observation(name,ROOT/'figure_qa'/f'SYNTHETIC_{name}.tif',phase,start,end,start,end,end))
    rng = np.random.default_rng(13)
    pixels = rng.uniform(0,1,256).astype('float32')
    values = np.array([.65+.08*pixels+.004*np.sin(i) if i<8 else .28+.05*pixels
                       if i==8 else .38+.055*(i-9)+.04*pixels for i in range(17)],dtype='float32')
    sample = Sample(256,13)
    sample.add(values,np.arange(256),np.zeros(256,dtype=int))
    zones = {}
    for i in range(2):
        zone = Sample(128,13+i)
        zone.add(values[:,i*128:(i+1)*128],np.arange(128),np.zeros(128,dtype=int))
        zones[f'SYNTHETIC zone {i+1}'] = zone
    fraction = metrics.calculate_fraction(values,obs,config).fraction
    return config,obs,sample,zones,fraction


def check_layout(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    page = fig.bbox
    for ax in fig.axes:
        assert not ax.get_title()
        xlabels = [a for a,x in zip(ax.get_xticklabels(),ax.get_xticks())
                   if ax.get_xlim()[0]<=x<=ax.get_xlim()[1]]
        ylabels = [a for a,y in zip(ax.get_yticklabels(),ax.get_yticks())
                   if min(ax.get_ylim())<=y<=max(ax.get_ylim())]
        artists = [ax.xaxis.label,ax.yaxis.label,*ax.texts,*xlabels,*ylabels]
        if ax.get_legend():
            artists.append(ax.get_legend())
        for artist in artists:
            if not artist.get_visible() or (hasattr(artist,'get_text') and not artist.get_text()):
                continue
            box = artist.get_window_extent(renderer)
            assert box.x0>=page.x0-1 and box.x1<=page.x1+1, f'Horizontal clipping: {artist}'
            assert box.y0>=page.y0-1 and box.y1<=page.y1+1, f'Vertical clipping: {artist}'
        labels = [a for a in xlabels if a.get_visible()]
        boxes = [a.get_window_extent(renderer) for a in labels]
        for left,right in zip(boxes,boxes[1:]):
            assert left.x1+1<right.x0, 'Overlapping x tick labels'
    for artist in [*fig.texts,*fig.legends]:
        box = artist.get_window_extent(renderer)
        assert box.x0>=page.x0-1 and box.x1<=page.x1+1, 'Figure note/legend clipping'
        assert box.y0>=page.y0-1 and box.y1<=page.y1+1, 'Figure note/legend clipping'
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--synthetic17',action='store_true')
    args = parser.parse_args()
    if args.synthetic17:
        config,obs,sample,zones,fraction = synthetic_fixture()
        out = ROOT/'figure_qa/synthetic17_component'
        out.mkdir(parents=True,exist_ok=True)
    else:
        config,obs = c,validate(c)
        out = Path(c.Output_Folder)
        sample,zones,fraction = read_cache(out,c,obs,metrics,Sample)
    apply_publication_style()
    data = prepare(sample,zones,obs,config,fraction,metrics)
    np.testing.assert_array_equal(data.anomaly,sample.values-data.metrics.baseline)
    for phase,group in zip(('pre','event','post'),data.groups):
        indices = [i for i,o in enumerate(obs) if o.phase==phase]
        np.testing.assert_array_equal(group,np.median(sample.values[indices],axis=0))
    event = next(i for i,o in enumerate(obs) if o.phase=='event')
    expected = (sample.values[data.post_indices][:,data.eligible]-sample.values[event,data.eligible])/data.metrics.loss[data.eligible]
    np.testing.assert_array_equal(data.normalized,expected)
    fig = main_figure(data)
    np.testing.assert_array_equal(fig.axes[0].lines[0].get_ydata(),data.table.mean_exact)
    np.testing.assert_array_equal(fig.axes[0].lines[0].get_xdata(),mdates.date2num([o.end for o in obs]))
    for line,o in zip(fig.axes[1].lines[::3],obs):
        np.testing.assert_array_equal(line.get_xdata(),mdates.date2num([o.start,o.end]))
    check_layout(fig)
    check_layout(main_figure(data,single=True))
    check_layout(overview_figure(data))
    for panel in (anomaly_panel,distribution_panel,recovery_panel,fraction_panel):
        check_layout(diagnostic_figure(data,panel))
    if zones:
        check_layout(heatmap(data,zones))
        check_layout(zone_series(data,zones))
    if args.synthetic17:
        render_portfolio(sample,zones,obs,out,config,fraction,metrics)
        write_cache(out,sample,zones,obs,config,fraction,metrics)
    report = {'status':'PASS','observations':len(obs),'synthetic_component_test':args.synthetic17,
              'exact_statistics':True,'phase_aggregation':True,'normalization_unchanged':True,
              'chronological_positions_and_intervals':True,'titles_empty':True,
              'layout':'All figures checked for page clipping and tick overlap'}
    (out/'journal_test_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
