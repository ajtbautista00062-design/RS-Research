"""Generate complete provenance alongside every applicable figure."""
from .style import DESIGN

DESCRIPTIONS = {
    'main': ('Interferometric coherence across acquisition pairs',
             'Filled markers and a solid line show the exact spatial mean over all common valid pixels; open markers and a dashed line show the sampled spatial median. The neutral band is the sampled P25–P75 spatial IQR. Circles, squares, and triangles denote pre, co-event, and post pairs. The separate coverage strip shows each full acquisition interval, keeping the statistical panel uncluttered.'),
    'anomaly': ('Coherence anomalies relative to each pixel’s pre-event reference',
                'The line shows the sampled median of pixel-level anomalies; the band is the sampled spatial IQR. Each pixel’s baseline is the existing pre-event median (the single pre value in the simulation). The horizontal reference at zero and dashed baseline-minus-tolerance level use the original configuration.'),
    'distribution': ('Spatial coherence distributions by disaster phase',
                     'Violin shapes depict the existing per-pixel phase values; embedded boxes show quartiles and the median. Pre and post groups use the pixel temporal median across their respective phase observations; the co-event group uses the event pair. Thus the 8/1/8 archive is not treated as unequal independent sample replicates. The simulation uses its single observation per phase. Density smoothing is for the spatial distribution only, never temporal interpolation.'),
    'recovery': ('Observed post-event coherence response',
                 'The plotted median and P25–P75 summarize the existing normalized fractions (post coherence minus event coherence)/(pre baseline minus event coherence), restricted to the original loss-eligible sampled pixels. Values are not clipped. The dotted line at one denotes the pre-event baseline; the dashed line is the median pixel-specific return threshold. One post observation is shown as a point with a spatial IQR bar, not a recovery trajectory.'),
    'fraction': ('Distribution of end-state coherence-return fraction',
                 'Histogram bins retain the original automatic-bin rule and include only finite values supplied by the existing recovery-fraction calculation. Fractions are unclipped. Reference lines at zero and one denote no net regain and return to the pre-event baseline; these are coherence-response proxies.'),
    'heatmap': ('Zonal coherence across actual observation dates',
                'Cells show the existing exact zonal means on fixed common-pixel support within each real zone. Column positions follow actual configured pair-reference dates; cell boundaries are adjacent-date midpoints, not additional measurements. The color scale uses the theoretical coherence range 0–1. Zones may have different pixel supports; no synthetic geographic zones are introduced.'),
    'zones': ('Observed zonal coherence series',
              'Each line shows the existing exact mean on a zone’s fixed complete-case pixel support. Actual date spacing is retained. Zone means are descriptive spatial summaries, not independent pixel-level temporal replicates.'),
    'overview': ('Complementary coherence diagnostics',
                 'Panels show (a) exact spatial mean, sampled median and spatial IQR; (b) sampled pixel-level coherence anomalies; (c) existing phase distributions; and (d) loss-eligible normalized post-event response. The associated standalone figure provides the full acquisition-interval strip. A single post-event point cannot establish a recovery trend.'),
}


def make_caption(data, kind, number, font, metadata):
    title, explanation = DESCRIPTIONS[kind]
    sample = data.sample
    periods = '; '.join(f'{o.name}: {o.start:%Y-%m-%d} to {o.end:%Y-%m-%d} '
                        f'({(o.end-o.start).total_seconds()/86400:g} days)' for o in data.obs)
    sources = '\n'.join(f'- `{o.path.name}`' for o in data.obs)
    limit = ('Only one pre and one post pair are available: baseline variability, recovery trend, duration and persistent recovery cannot be established.'
             if len(data.obs) == 3 else
             'The observations are sparse, spatial pixels are autocorrelated, and coherence return does not verify structural reconstruction. No new temporal smoothing or fitted curve is introduced.')
    origin = 'SYNTHETIC COMPONENT TEST; NOT AN ACTUAL RESEARCH RESULT.\n\n' if getattr(data.c, 'SYNTHETIC_VISUALIZATION_TEST', False) else ''
    return f'''{origin}**Figure {number}. {title}.** {explanation}

Common valid-pixel population: {sample.count:,}. Existing deterministic statistical sample: {sample.values.shape[1]:,} pixels. Post-response subset: {int(data.eligible.sum()):,} sampled loss-eligible pixels. The IQR is spatial variability, not a confidence interval for the mean or measurement uncertainty. Straight line connections and band edges connect observed summaries only; they do not add observations. Acquisition intervals: {periods}.

Typhoon Noru reference date: {DESIGN['event_date']}. A thin dashed event line is included only when the date lies inside the displayed time domain; no impact-window shading is used. The main figure uses pair-end dates. Anomaly, post-response and zonal views preserve the existing `{data.c.REFERENCE_DATE}` pair-reference convention. Coherence axes use an expanded vertical scale where appropriate; theoretical coherence range is 0–1. Actual displayed axis limits: {metadata['axis_limits']}. Normalized response and anomaly axes are not restricted to [0,1]. {limit}

Source inputs (sensor, orbit, polarization and processing software lineage are not independently verified by plotting code):

{sources}

Existing raster alignment, coherence scaling, masks, sample selection, baseline calculations, classification, fractions and zonal statistics are unchanged. Configured resampling: `{data.c.RESAMPLING}`. Exact-zero exclusion setting: `{getattr(data.c, 'ZERO_IS_INVALID', 'not configured')}`. Acquisition comparability confirmed in config: `{data.c.ACQUISITION_COMPARABILITY_CONFIRMED}`. No geographic zone inference or physical-recovery claim is made. Export: {font}, embedded TrueType PDF, editable SVG text, 600-DPI PNG; physical size {metadata['figure_inches']} inches.
'''
