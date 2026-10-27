"""Real-zone plots with actual chronological spacing and bounded page density."""
import textwrap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from .style import DESIGN as D, date_axis, clean_axis, frame_values, event_line


def heatmap(data, selected):
    names = list(selected)
    values = np.array([sample.sums/sample.count for sample in selected.values()])
    dates = mdates.date2num([o.reference for o in data.obs])
    boundaries = np.r_[dates[0]-(dates[1]-dates[0])/2,
                       (dates[:-1]+dates[1:])/2, dates[-1]+(dates[-1]-dates[-2])/2]
    fig = plt.figure(figsize=(7, max(3.6, len(names)*.35+1.8)))
    ax = fig.add_axes([.28, .24, .59, .59])
    image = ax.pcolormesh(boundaries, np.arange(len(names)+1), values,
                          cmap='cividis', vmin=0, vmax=1, shading='flat', rasterized=False)
    ax.set_yticks(np.arange(len(names))+.5, ['\n'.join(textwrap.wrap(name, 21)) for name in names])
    ax.invert_yaxis()
    ax.tick_params(axis='y', length=0)
    date_axis(ax, data.obs, reference=data.c.REFERENCE_DATE)
    ax.set_xlim(boundaries[0], boundaries[-1])
    event_line(ax, label=True)
    color_ax = fig.add_axes([.9, .24, .018, .59])
    fig.colorbar(image, cax=color_ax, label='Coherence (0–1)')
    fig.text(.28, .04, 'Each zone retains its own fixed common-pixel support.', fontsize=8, color=D['reference'])
    return fig


def zone_series(data, selected):
    fig = plt.figure(figsize=(7, 4.0))
    ax = fig.add_axes([.12, .36, .85, .46])
    values = []
    for i, (name, sample) in enumerate(selected.items()):
        mean = sample.sums/sample.count
        values.append(mean)
        ax.plot([o.reference for o in data.obs], mean,
                color=plt.get_cmap('cividis')(.12+.7*i/max(len(selected)-1, 1)),
                ls=('-', '--', ':', '-.')[i%4], marker=('o','s','^','D')[i%4],
                ms=3.4, lw=1.2, label=name)
    clean_axis(ax)
    frame_values(ax, values, theoretical=True)
    date_axis(ax, data.obs, reference=data.c.REFERENCE_DATE)
    event_line(ax, label=True)
    ax.set_ylabel('Zonal mean coherence (unitless)')
    ax.legend(loc='upper left', bbox_to_anchor=(0, -.37), frameon=False,
              ncol=2, fontsize=7.5, borderaxespad=0)
    return fig
