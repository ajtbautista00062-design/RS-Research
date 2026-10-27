"""One design system used identically by both standalone research projects."""
import numpy as np
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib import font_manager
from matplotlib.ticker import MaxNLocator

DESIGN = {
    'mean': '#262b30', 'median': '#496b79', 'iqr': '#8c969e',
    'pre': '#687882', 'event': '#98684b', 'post': '#3e6c79',
    'grid': '#e5e8ea', 'reference': '#7c8186',
    'markers': {'pre': 'o', 'event': 's', 'post': '^'},
    'labels': {'pre': 'Pre-disaster', 'event': 'Co-event', 'post': 'Post-disaster'},
    'line_width': 1.4, 'marker_size': 4.7, 'iqr_alpha': .18,
    'width': 7., 'dpi': 600, 'event_date': '2022-09-25',
    'stem': 'fig_coherence_timeseries_20220925',
}


def apply_publication_style():
    import logging
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    installed = {font.name for font in font_manager.fontManager.ttflist}
    family = next((n for n in ('Arial', 'Helvetica', 'DejaVu Sans') if n in installed), 'DejaVu Sans')
    mpl.rcParams.update({
        'font.family': family, 'font.size': 8.5, 'axes.labelsize': 9.5,
        'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'axes.edgecolor': '#495058', 'axes.linewidth': .65,
        'axes.spines.top': False, 'axes.spines.right': False,
        'xtick.major.width': .65, 'ytick.major.width': .65,
        'xtick.major.size': 3, 'ytick.major.size': 3,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    })
    return family


def clean_axis(ax):
    ax.grid(axis='y', color=DESIGN['grid'], lw=.5)
    ax.set_axisbelow(True)
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.spines[['top', 'right']].set_visible(False)


def frame_values(ax, values, include=(), theoretical=False):
    finite = np.asarray(values).ravel()
    finite = finite[np.isfinite(finite)]
    finite = np.r_[finite, include]
    if not finite.size:
        return
    low, high = float(finite.min()), float(finite.max())
    padding = max(high-low, .025)*.14
    ax.set_ylim(max(0., low-padding) if theoretical else low-padding,
                min(1., high+padding) if theoretical else high+padding)


def date_axis(ax, obs, reference='end', coverage=False, compact=False):
    dates = [o.end if reference == 'end' else o.reference for o in obs]
    start = mdates.date2num(obs[0].start if coverage else dates[0])
    end = mdates.date2num(dates[-1])
    span = max(end-start, 1.)
    target = 2 if compact else 6
    options = (6, 7, 12, 14, 18, 21, 28, 42, 56, 84, 112)
    step = next((s for s in options if span/s <= target), 112)
    ticks = np.arange(start, end+.01, step)
    if span < 2:
        ticks = [start]
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    pad = max(.8, span*.02)
    ax.set_xlim(start-pad, end+pad)
    ax.set_xlabel('Pair end date (UTC)' if reference == 'end' else 'Pair reference date (UTC)', labelpad=7)
    ax.tick_params(axis='x', labelrotation=0)
    if compact:
        labels = ax.get_xticklabels()
        if labels:
            labels[0].set_ha('left')
            labels[-1].set_ha('right')


def event_line(ax, label=False):
    event = mdates.datestr2num(DESIGN['event_date'])
    if ax.get_xlim()[0] <= event <= ax.get_xlim()[1]:
        ax.axvline(event, color=DESIGN['reference'], lw=.7, ls=(0, (3, 3)), zorder=0)
        if label:
            ax.annotate('Typhoon Noru · 25 Sep 2022', (event, 1),
                        xycoords=('data', 'axes fraction'), xytext=(0, 10),
                        textcoords='offset points', ha='center', va='bottom', fontsize=8,
                        color=DESIGN['reference'], annotation_clip=False)
