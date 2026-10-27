"""Publication typography, lossless figure export, and evidence-based captions."""
from pathlib import Path
import matplotlib as mpl
from matplotlib import font_manager

PUBLICATION = {
    'colors': {'mean': '#25282b', 'median': '#4b6570', 'iqr': '#92979c',
               'interval': '#777d82', 'event': '#666a6d', 'grid': '#e5e7e9'},
    'markers': {'pre': 'o', 'event': 's', 'post': '^'},
    'phase_labels': {'pre': 'Pre-disaster', 'event': 'Co-event', 'post': 'Post-disaster'},
    'marker_size': 4.7, 'mean_width': 1.4, 'median_width': 1.0,
    'iqr_alpha': .22, 'interval_width': .65, 'grid_width': .5,
    'date_tick_days': 6, 'single_date_tick_days': 18,
    'legend_anchor': (0., -.32), 'legend_location': 'upper left',
    'figure_inches': (7., 3.6), 'single_inches': (3.5, 3.8),
    'dpi': 600, 'event_date': '2022-09-25',
    'stem': 'fig_coherence_timeseries_20220925',
}


def apply_publication_style() -> str:
    available = {font.name for font in font_manager.fontManager.ttflist}
    family = next((name for name in ('Arial', 'Helvetica', 'DejaVu Sans') if name in available),
                  'DejaVu Sans')
    mpl.rcParams.update({
        'font.family': family, 'font.size': 8.5, 'axes.labelsize': 9.5,
        'xtick.labelsize': 8, 'ytick.labelsize': 8.5, 'legend.fontsize': 8,
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'savefig.facecolor': 'white', 'axes.edgecolor': '#45484b',
        'axes.linewidth': .65, 'xtick.major.width': .65, 'ytick.major.width': .65,
        'xtick.major.size': 3., 'ytick.major.size': 3.,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'axes.spines.top': False, 'axes.spines.right': False,
    })
    return family


def export_publication_figure(fig, out: Path, stem: str) -> None:
    """Keep the specified physical page size; tighten layout inside that page.

    The full-figure invisible bounds prevent bbox_inches='tight' from changing
    journal column width or clipping the external legend and notes.
    """
    from matplotlib.patches import Rectangle
    bounds = Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                       facecolor='none', edgecolor='none', linewidth=0)
    fig.add_artist(bounds)
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(out / f'{stem}.{ext}', dpi=PUBLICATION['dpi'],
                    bbox_inches='tight', pad_inches=0, bbox_extra_artists=[bounds])


def caption_text(table, obs, c, ylim, font: str) -> str:
    phase_labels = PUBLICATION['phase_labels']
    periods = '; '.join(
        f"{phase_labels[o.phase].lower()}: {o.start:%Y-%m-%d} to {o.end:%Y-%m-%d} "
        f"({(o.end-o.start).total_seconds()/86400:g} days)" for o in obs)
    count = int(table.complete_case_pixels.iloc[0])
    sampled = int(table.sample_pixels.iloc[0])
    paths = '\n'.join(f'- {o.name}: `{o.path.name}`' for o in obs)
    return f'''**Figure 1. Three-pair interferometric coherence around Typhoon Noru, Infanta, Quezon.**

The coherence periods are {periods}. Solid charcoal lines and filled markers show the exact spatial mean over {count:,} common valid pixels. Dashed lines and open markers show the spatial median of the existing deterministic sample ({sampled:,} pixels); the neutral gray band shows that sample's spatial interquartile range (P25–P75). The IQR describes spatial variability, not a confidence interval for the mean or measurement uncertainty. Circles, squares, and triangles identify the pre-disaster, co-event, and post-disaster pairs, respectively.

Each observation summarizes two SAR acquisitions rather than an instantaneous measurement. Points are positioned at the later (secondary/pair-end) date; the thin capped horizontal segments above the IQR show each full acquisition interval. The thin vertical dashed line marks the supplied Typhoon Noru reference date, 2022-09-25 (UTC date reference). No impact-window shading is used. Straight connections between the three observed statistics and the band edges are visual guides only; no intermediate measurements, fitted recovery trend, duration, persistent return, or structural recovery are inferred.

The displayed vertical limits are {ylim[0]:.6f}–{ylim[1]:.6f}; this is an expanded vertical scale within the theoretical coherence range of 0–1. The limits include all three means, medians, and sample quartiles with padding. The primary figure measures 7.0 × 3.6 inches; the optional single-column version measures 3.5 × 3.8 inches. PNG export is 600 DPI. PDF uses embedded TrueType fonts and SVG retains editable text; selected font: {font}.

**Source and processing provenance.** User-described Sentinel-1 InSAR coherence GeoTIFF inputs (sensor identity has not been independently verified from SLC records):

{paths}

Pair dates are supplied by the user and cross-checked against the configured filenames. Existing analysis statistics, common-pixel support, seeded priority sample, reference-grid alignment, masks, coherence scaling, and analytical methods are unchanged. Rendering reads the existing statistics rather than estimating new values. Configured resampling is `{c.RESAMPLING}`; exact-zero exclusion is `{c.ZERO_IS_INVALID}`. No study-area polygons are configured unless specified in config.py; the common-pixel support is not a validated built-up mask. Orbit, polarization, acquisition time of day, processing software versions, and acquisition comparability are unverified and are not inferred here.
'''
