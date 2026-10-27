"""Complementary diagnostics using the unchanged original statistics."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from .style import DESIGN as D, clean_axis, date_axis, frame_values, event_line


def anomaly_panel(ax, data, compact=False):
    q25, median, q75 = np.quantile(data.anomaly, [.25, .5, .75], axis=1)
    dates = mdates.date2num([o.reference for o in data.obs])
    ax.fill_between(dates, q25, q75, color=D['iqr'], alpha=D['iqr_alpha'], lw=0)
    ax.plot(dates, median, color=D['mean'], lw=1.4)
    for i, o in enumerate(data.obs):
        ax.plot(dates[i], median[i], marker=D['markers'][o.phase], color=D[o.phase],
                ms=4.7, ls='none', zorder=4)
    ax.axhline(0, color=D['reference'], lw=.7)
    ax.axhline(-data.c.RECOVERY_TOLERANCE, color=D['median'], lw=.85, ls='--')
    clean_axis(ax)
    frame_values(ax, np.r_[q25, q75], include=[0, -data.c.RECOVERY_TOLERANCE])
    date_axis(ax, data.obs, reference=data.c.REFERENCE_DATE, compact=compact)
    event_line(ax, label=not compact)
    ax.set_ylabel('Coherence anomaly\n(relative to pre baseline)' if compact else 'Coherence anomaly (relative to pre baseline)')


def distribution_panel(ax, data, compact=False):
    for i, (group, phase) in enumerate(zip(data.groups, ('pre', 'event', 'post')), start=1):
        if np.ptp(group) > 0 and group.size > 1:
            violin = ax.violinplot([group], positions=[i], widths=.7, showextrema=False)
            violin['bodies'][0].set(facecolor=D[phase], edgecolor=D[phase], linewidth=.7, alpha=.13)
        ax.boxplot([group], positions=[i], widths=.13, showfliers=False, whis=(25, 75),
                   patch_artist=True, boxprops={'facecolor':'white','edgecolor':D[phase],'linewidth':.8},
                   medianprops={'color':D[phase], 'linewidth':1.3},
                   whiskerprops={'color':D[phase]}, capprops={'color':D[phase]})
    ax.set_xticks([1, 2, 3], ['Pre-disaster', 'Co-event', 'Post-disaster'])
    ax.set_xlim(.4, 3.6)
    ax.set_ylabel('Coherence (unitless)')
    clean_axis(ax)
    # A distribution plot displays its full finite support, never truncates its tails.
    frame_values(ax, np.concatenate(data.groups), theoretical=True)
    if not compact:
        ax.set_xlabel('One value per sampled pixel per phase', labelpad=7)


def recovery_panel(ax, data, compact=False):
    if not data.eligible.any():
        ax.text(.5, .5, 'No eligible sampled pixels', ha='center', va='center', transform=ax.transAxes)
        ax.set_axis_off()
        return
    q25, median, q75 = np.quantile(data.normalized, [.25, .5, .75], axis=1)
    obs = [data.obs[i] for i in data.post_indices]
    dates = mdates.date2num([o.reference for o in obs])
    if len(obs) > 1:
        ax.fill_between(dates, q25, q75, color=D['iqr'], alpha=D['iqr_alpha'], lw=0)
        ax.plot(dates, median, color=D['post'], lw=1.4, marker='^', ms=4.7)
    else:
        ax.errorbar(dates, median, yerr=np.vstack((median-q25, q75-median)),
                    fmt='^', color=D['post'], ms=5., lw=1., capsize=5)
    threshold = float(np.median(data.threshold))
    ax.axhline(1, color=D['reference'], lw=.8, ls=':')
    ax.axhline(threshold, color=D['median'], lw=.85, ls='--')
    clean_axis(ax)
    frame_values(ax, np.r_[q25, q75], include=[1, threshold])
    date_axis(ax, obs, reference=data.c.REFERENCE_DATE, compact=compact)
    ax.set_ylabel('Fraction of event loss\nregained (unclipped)' if compact else 'Fraction of event loss regained (unclipped)')
    if not compact:
        ax.annotate('Pre-event baseline', xy=(1, 1), xycoords=('axes fraction','data'),
                    xytext=(-4, 4), textcoords='offset points', ha='right', fontsize=8, color=D['reference'])
        ax.text(.02, .94, 'Single post pair; persistence unverified' if len(obs) == 1 else 'Observed post pairs; no temporal fit',
                transform=ax.transAxes, va='top', fontsize=8, color=D['reference'])


def fraction_panel(ax, data):
    values = data.fraction[np.isfinite(data.fraction)]
    ax.hist(values, bins='auto', color=D['post'], alpha=.65, edgecolor='white', linewidth=.45)
    ax.axvline(0, lw=.8, color=D['reference'], ls='--')
    ax.axvline(1, lw=.8, color=D['mean'], ls=':')
    ax.set(xlabel='Coherence-return fraction (unclipped)', ylabel='Sampled pixel count')
    clean_axis(ax)


def diagnostic_figure(data, panel):
    fig = plt.figure(figsize=(7, 3.6))
    ax = fig.add_axes([.12, .25, .85, .57])
    panel(ax, data)
    fig.text(.12, .045, 'Spatial summaries; IQR is not a confidence interval.',
             fontsize=8, color=D['reference'])
    return fig


def overview_figure(data):
    from .timeseries import statistical_series
    fig, axes = plt.subplots(2, 2, figsize=(7, 5.9))
    fig.subplots_adjust(left=.12, right=.97, bottom=.18, top=.92, wspace=.43, hspace=.72)
    handles = statistical_series(axes[0, 0], data, compact=True)
    event_line(axes[0, 0])
    anomaly_panel(axes[0, 1], data, compact=True)
    distribution_panel(axes[1, 0], data, compact=True)
    recovery_panel(axes[1, 1], data, compact=True)
    for letter, ax in zip('abcd', axes.ravel()):
        ax.text(0, 1.08, f'({letter})', transform=ax.transAxes, fontweight='bold', fontsize=9)
        ax.xaxis.label.set_fontsize(9)
        ax.yaxis.label.set_fontsize(9)
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(.12, .065),
               ncol=3, frameon=False, columnspacing=1.3, handlelength=2.0, borderaxespad=0)
    fig.text(.12, .025, 'Expanded coherence scales; theoretical range 0–1. Spatial IQR, not confidence intervals.',
             fontsize=8, color=D['reference'])
    return fig
