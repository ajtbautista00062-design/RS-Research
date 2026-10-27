"""Clear chronological statistics with acquisition support on a separate strip."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from .style import DESIGN as D, clean_axis, date_axis, frame_values, event_line


def statistical_series(ax, data, compact=False):
    table = data.table
    dates = mdates.date2num([o.end for o in data.obs])
    q25, median, q75 = [table[key].to_numpy() for key in ('q25_sample', 'median_sample', 'q75_sample')]
    ax.fill_between(dates, q25, q75, facecolor=D['iqr'], alpha=D['iqr_alpha'], lw=0)
    ax.plot(dates, table.mean_exact, lw=D['line_width'], color=D['mean'])
    ax.plot(dates, median, lw=1., ls=(0, (4, 2.5)), color=D['median'])
    for i, obs in enumerate(data.obs):
        marker = D['markers'][obs.phase]
        ax.plot(dates[i], table.mean_exact.iloc[i], marker=marker, ms=D['marker_size'],
                mfc=D['mean'], mec=D['mean'], mew=.7, ls='none', zorder=5)
        ax.plot(dates[i], median[i], marker=marker, ms=D['marker_size'],
                mfc='white', mec=D['median'], mew=.9, ls='none', zorder=5)
    clean_axis(ax)
    date_axis(ax, data.obs, coverage=not compact, compact=compact)
    frame_values(ax, table[['mean_exact', 'q25_sample', 'q75_sample']], theoretical=True)
    ax.set_ylabel('Interferometric coherence\n(unitless)' if compact else 'Interferometric coherence (unitless)')
    return [Line2D([], [], color=D['mean'], lw=1.4, label='Spatial mean'),
            Line2D([], [], color=D['median'], lw=1., ls='--', label='Spatial median'),
            Patch(facecolor=D['iqr'], alpha=D['iqr_alpha'], edgecolor='none', label='Spatial IQR (P25–P75)')]


def coverage_strip(ax, obs, compact=False):
    yvalues = {'pre': 2, 'event': 1, 'post': 0}
    for o in obs:
        start, end = mdates.date2num([o.start, o.end])
        y = yvalues[o.phase]
        ax.plot([start, end], [y, y], color=D[o.phase], lw=.8, solid_capstyle='butt')
        ax.plot(end, y, marker=D['markers'][o.phase], ms=3.4, color=D[o.phase], ls='none')
        ax.plot(start, y, marker='|', ms=3.8, color=D[o.phase], ls='none', mew=.7)
    ax.set_ylim(-.6, 2.6)
    ax.set_yticks([2, 1, 0], ['Pre', 'Co-event', 'Post'])
    ax.tick_params(axis='y', length=0, labelsize=7.5, pad=5)
    ax.spines[['left', 'top', 'right']].set_visible(False)
    date_axis(ax, obs, coverage=True, compact=compact)


def main_figure(data, single=False):
    fig = plt.figure(figsize=(3.5, 4.4) if single else (7, 4.1))
    left, right = (.215, .95) if single else (.12, .975)
    ax = fig.add_axes([left, .43, right-left, .43])
    handles = statistical_series(ax, data, compact=single)
    # Both panels use the same chronological domain, including the first pair start.
    date_axis(ax, data.obs, coverage=True, compact=single)
    event_line(ax, label=True)
    ax.set_xlabel('')
    ax.tick_params(axis='x', labelbottom=False, bottom=False)
    support = fig.add_axes([left, .265, right-left, .105], sharex=ax)
    coverage_strip(support, data.obs, compact=single)
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(left, .065),
               ncol=1 if single else 3, frameon=False, handlelength=2.1,
               columnspacing=1.4, labelspacing=.35, borderaxespad=0)
    note = 'Expanded vertical scale; theoretical coherence range: 0–1'
    if single:
        note = 'Expanded vertical scale;\ntheoretical coherence range: 0–1'
    fig.text(left, .018, note, fontsize=8, color=D['reference'], va='bottom')
    return fig
