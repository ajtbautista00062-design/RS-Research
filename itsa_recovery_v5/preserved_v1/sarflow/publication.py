"""Three-pair publication figure consuming the existing statistical table unchanged."""
from pathlib import Path
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator, FormatStrFormatter

from .publication_export import (PUBLICATION, apply_publication_style,
                                 export_publication_figure, caption_text)
from .temporal import timestamp


def plot_three_pair_coherence(table, obs, out: Path, c, single_column: bool = False) -> dict:
    required = ['mean_exact', 'q25_sample', 'median_sample', 'q75_sample']
    if len(obs) != 3 or table.observation.tolist() != [o.name for o in obs]:
        raise ValueError('Statistics must match the three configured observations in order')
    values = table[required].to_numpy(dtype='float64')
    if not np.isfinite(values).all() or not ((values >= 0) & (values <= 1)).all():
        raise ValueError('Publication figure requires finite coherence statistics in [0,1]')
    if not (table.q25_sample.le(table.median_sample) & table.median_sample.le(table.q75_sample)).all():
        raise ValueError('Invalid percentile ordering')
    if table.complete_case_pixels.nunique() != 1 or table.sample_pixels.nunique() != 1:
        raise ValueError('Fixed common-pixel support required')
    out.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context():
        font = apply_publication_style()
        cfg, colors = PUBLICATION, PUBLICATION['colors']
        size = cfg['single_inches'] if single_column else cfg['figure_inches']
        fig = plt.figure(figsize=size)
        ax = fig.add_axes([.205, .43, .76, .43] if single_column else [.105, .38, .87, .48])
        dates = [o.end for o in obs]  # explicit pair end, independent of sensitivity reference setting
        x = mdates.date2num(dates)
        mean, q25, median, q75 = values.T
        ax.fill_between(x, q25, q75, facecolor=colors['iqr'], alpha=cfg['iqr_alpha'],
                        linewidth=0, zorder=1)
        ax.plot(x, mean, color=colors['mean'], lw=cfg['mean_width'], zorder=4)
        ax.plot(x, median, color=colors['median'], lw=cfg['median_width'],
                linestyle=(0, (4, 2.5)), zorder=3)
        span = max(float(values.max()-values.min()), .01)
        interval_y = q75 + span*.14
        for j, o in enumerate(obs):
            marker = cfg['markers'][o.phase]
            ax.plot(x[j], mean[j], marker=marker, ms=cfg['marker_size'],
                    mfc=colors['mean'], mec=colors['mean'], mew=.7, linestyle='none', zorder=6)
            ax.plot(x[j], median[j], marker=marker, ms=cfg['marker_size'],
                    mfc='white', mec=colors['median'], mew=.9, linestyle='none', zorder=5)
            start = mdates.date2num(o.start)
            ax.errorbar(x[j], interval_y[j], xerr=np.array([[x[j]-start], [0.]]),
                        fmt='none', ecolor=colors['interval'], elinewidth=cfg['interval_width'],
                        capsize=2, capthick=cfg['interval_width'], zorder=2)
            label = cfg['phase_labels'][o.phase]
            ax.text((start+x[j])/2, interval_y[j]+span*.065, label,
                    ha='center', va='bottom', fontsize=7.5, color=colors['interval'])
        ylim = (max(0., float(values.min())-span*.16),
                min(1., float(interval_y.max())+span*.28))
        ax.set_ylim(*ylim)
        ax.set_xlim(mdates.date2num(obs[0].start)-.9, x[-1]+.9)
        step = cfg['single_date_tick_days'] if single_column else cfg['date_tick_days']
        ticks = np.arange(mdates.date2num(obs[0].start), x[-1]+.01, step)
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        if single_column:
            labels = ax.get_xticklabels()
            labels[0].set_ha('left')
            labels[-1].set_ha('right')
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.set_xlabel('Pair end date (UTC)', labelpad=7)
        ax.set_ylabel('Interferometric coherence (unitless)', labelpad=8)
        ax.grid(axis='y', which='major', color=colors['grid'], lw=cfg['grid_width'])
        ax.set_axisbelow(True)
        event = mdates.date2num(timestamp(cfg['event_date']))
        interval_clearance = (float(interval_y.min())-span*.04-ylim[0])/(ylim[1]-ylim[0])
        ax.axvline(event, ymax=interval_clearance, color=colors['event'], lw=.7,
                   ls=(0, (3, 3)), zorder=0)
        ax.annotate('Typhoon Noru · 25 Sep 2022', xy=(event, 1), xycoords=('data', 'axes fraction'),
                    xytext=(0, 16), textcoords='offset points', ha='center', va='bottom',
                    fontsize=8, color=colors['event'], annotation_clip=False)
        handles = [Line2D([], [], color=colors['mean'], lw=cfg['mean_width'], label='Spatial mean (all common pixels)'),
                   Line2D([], [], color=colors['median'], lw=cfg['median_width'], ls='--', label='Spatial median (sample)'),
                   Patch(facecolor=colors['iqr'], alpha=cfg['iqr_alpha'], edgecolor='none', label='Spatial IQR (P25–P75)')]
        ax.legend(handles=handles, loc=cfg['legend_location'], bbox_to_anchor=cfg['legend_anchor'],
                  ncol=1, frameon=False, borderaxespad=0, handlelength=2.4, labelspacing=.5)
        fig.text(.105 if not single_column else .205, .045,
                 'Expanded vertical scale; theoretical coherence range: 0–1' if not single_column
                 else 'Expanded vertical scale;\ntheoretical coherence range: 0–1',
                 fontsize=8, color=colors['event'], va='bottom')
        stem = cfg['stem'] + ('_single_column' if single_column else '')
        export_publication_figure(fig, out, stem)
        details = {
            'observation_dates': [o.end.isoformat() for o in obs],
            'pair_intervals': [[o.start.isoformat(), o.end.isoformat()] for o in obs],
            'event_date': cfg['event_date'], 'ylim': list(ylim), 'font_family': font,
            'font_path': mpl.font_manager.findfont(font), 'figure_inches': list(size),
            'png_dpi': cfg['dpi'], 'statistics': table.to_dict(orient='records'),
            'date_tick_days': step,
            'no_title': ax.get_title() == '', 'iqr_is_confidence_interval': False,
        }
        (out/f'{stem}_caption.md').write_text(caption_text(table, obs, c, ylim, font), encoding='utf-8')
        (out/f'{stem}_data.json').write_text(json.dumps(details, indent=2), encoding='utf-8')
        plt.close(fig)
    return details
