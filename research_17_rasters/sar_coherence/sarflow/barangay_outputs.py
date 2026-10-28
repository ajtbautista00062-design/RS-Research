"""Reuse existing rasters; exact streamed zonal summaries and explicit aggregate fits."""
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_window
from rasterio.errors import WindowError
from rasterio.windows import Window, transform as window_transform
from .spatial import windows
from .recovery import assess
from .recovery_io import METRICS, COLORS
from .temporal import day_offsets
from .enhancement_reports import Comparison, connected_region_report, render_comparison


class DiskValues:
    """Scratch-only Float32 values for exact order statistics without full RAM stacks."""
    def __init__(self, directory):
        self.directory = Path(directory)
        self.counts, self.sums = {}, {}

    def add(self, key, values):
        a = np.asarray(values, dtype='float32')
        a = a[np.isfinite(a)]
        if not a.size:
            return
        with (self.directory / (key+'.bin')).open('ab') as stream:
            a.tofile(stream)
        self.counts[key] = self.counts.get(key, 0)+a.size
        self.sums[key] = self.sums.get(key, 0.)+float(a.sum(dtype='float64'))

    def stats(self, key):
        n = self.counts.get(key, 0)
        if not n:
            return dict(count=0, mean=np.nan, q25=np.nan, median=np.nan, q75=np.nan)
        a = np.memmap(self.directory / (key+'.bin'), mode='r+', dtype='float32', shape=(n,))
        positions = np.array([.25,.5,.75])*(n-1)
        lower, upper = np.floor(positions).astype(int), np.ceil(positions).astype(int)
        a.partition(np.unique(np.r_[lower, upper]))  # only temporary scratch data modified
        quantiles = a[lower].astype(float)+(a[upper]-a[lower])*(positions-lower)
        a.flush()
        del a
        return dict(count=n, mean=self.sums[key]/n, q25=quantiles[0], median=quantiles[1], q75=quantiles[2])


def export_barangays(paths, records, layout, obs, c):
    if not records:
        return
    reports = []
    with ExitStack() as context:
        sources = {k: context.enter_context(rasterio.open(p)) for k,p in paths.items()}
        owner = sources['Barangay_Owner']
        for zone_id, record in enumerate(records, 1):
            try:
                crop = geometry_window(owner, [record.geometry]).intersection(Window(0,0,owner.width,owner.height))
            except WindowError:
                reports.append(dict(code=record.code, name=record.name, status='outside raster extent'))
                continue
            # A geometry intersection need not contain any owned pixel centers.
            found = any(np.any(owner.read(1, window=Window(w.col_off+crop.col_off, w.row_off+crop.row_off,
                        w.width,w.height)) == zone_id) for w in windows(int(crop.width),int(crop.height),c.WINDOW_SIZE))
            if not found:
                reports.append(dict(code=record.code, name=record.name, status='no owned pixel centers'))
                continue
            raster_dir = layout.rasters / 'Barangay' / record.folder
            graph_dir = layout.graphs / 'Barangay' / record.folder
            table_dir = layout.tables / 'Barangay' / record.folder
            for folder in (raster_dir, graph_dir, table_dir):
                folder.mkdir(parents=True, exist_ok=False)
            export_zone(sources, crop, zone_id, record, raster_dir, graph_dir, table_dir, obs, c)
            reports.append(dict(code=record.code, name=record.name, folder=record.folder, status='exported'))
    pd.DataFrame(reports).to_csv(layout.tables / 'Barangay_Export_Status.csv', index=False)


def export_zone(sources, crop, zone_id, record, raster_dir, graph_dir, table_dir, obs, c):
    owner = sources['Barangay_Owner']
    pixel_area = abs(owner.transform.a*owner.transform.e)*owner.crs.linear_units_factor[1]**2
    summary, inside_count = Comparison(), 0
    region_counts = {}
    with ExitStack() as context, TemporaryDirectory(prefix='.statistics_', dir=table_dir) as scratch:
        values = DiskValues(scratch)
        writers = {}
        for name, src in sources.items():
            if name in ('Barangay_Owner', 'Spatial_Strata', 'Aligned_Coherence_Observations', 'Raw_ITSA_Metrics'):
                continue
            profile = src.profile.copy()
            profile.update(width=int(crop.width), height=int(crop.height), transform=window_transform(crop, src.transform))
            dst = context.enter_context(rasterio.open(raster_dir / (name+'.tif'), 'w', **profile))
            dst.update_tags(**src.tags())
            dst.update_tags(barangay=record.name, barangay_code=record.code,
                            ownership='first sorted code, pixel centers; models not refitted for clipped rasters')
            for b in range(1,src.count+1):
                dst.set_band_description(b, src.descriptions[b-1])
                dst.update_tags(b, **src.tags(b))
            if name.endswith(('Classes','Classification')):
                dst.write_colormap(1, COLORS)
            writers[name] = dst
        for win in windows(int(crop.width), int(crop.height), c.WINDOW_SIZE):
            source_win = Window(crop.col_off+win.col_off, crop.row_off+win.row_off, win.width, win.height)
            member = owner.read(1, window=source_win) == zone_id
            inside_count += int(member.sum())
            for name, dst in writers.items():
                data = sources[name].read(window=source_win)
                data[:, ~member] = dst.nodata
                dst.write(data, window=win)
            original = sources['Original_ITSA_Recovery_Classification'].read(1, window=source_win)
            enhanced = sources['Enhanced_ITSA_Recovery_Classification'].read(1, window=source_win)
            score = sources['Evidence_Confidence'].read(1, window=source_win, masked=True).filled(np.nan)
            builtup = sources['Builtup_Status'].read(1, window=source_win)
            summary.add(original, enhanced, score, member, builtup)
            stack = sources['Aligned_Coherence_Observations'].read(window=source_win, masked=True).filled(np.nan)
            for i in range(len(obs)):
                values.add(f'date_{i}', stack[i][member])
            metrics = sources['Enhanced_ITSA_Metrics'].read(window=source_win, masked=True).filled(np.nan)
            for i,key in enumerate(METRICS):
                values.add(key, metrics[i][member & (enhanced > 0)])
            regions = sources['Spatial_Recovery_Regions'].read(1, window=source_win)
            for gid, count in zip(*np.unique(regions[member & (regions > 0)], return_counts=True)):
                region_counts[int(gid)] = region_counts.get(int(gid),0)+int(count)
        rows = []
        for i,o in enumerate(obs):
            stats = values.stats(f'date_{i}')
            fraction = stats['count']/inside_count if inside_count else 0
            rows.append(dict(barangay=record.name, code=record.code, observation=o.name, phase=o.phase,
                             master=o.master, slave=o.slave, reference=o.reference, **stats,
                             valid_fraction=fraction, eligible=fraction >= c.MIN_PER_DATE_VALID_COVERAGE,
                             support='per-date available pixels; spatial cohort may change'))
        table = pd.DataFrame(rows)
        table.to_csv(table_dir / 'Coherence_Time_Series.csv', index=False)
        table.to_csv(table_dir / 'Valid_Data_Coverage.csv', index=False)
        metric_rows = [dict(metric=k, **values.stats(k)) for k in METRICS]
        pd.DataFrame(metric_rows).to_csv(table_dir / 'ITSA_Indicator_Statistics.csv', index=False)
        # Classifications use pixel/object models; this separate aggregate fit is descriptive.
        series = np.where(table.eligible, table['median'], np.nan).astype(float)
        result = assess(series[:,None], obs, c)
        pd.DataFrame([dict(barangay=record.name, model_scope='separate fit to available-pixel medians; not classification model',
                           **{k:float(v[0]) for k,v in result.metrics.items()}, recovery_class=int(result.classes[0]))]).to_csv(
                               table_dir / 'ITSA_Regression_Results.csv', index=False)
        graph_series(table, result, obs, graph_dir, c, record.name)
    summary.export(table_dir, pixel_area, record.name)
    render_comparison(summary, graph_dir, c)
    pd.DataFrame([dict(region_id=k, pixels_in_barangay=v, area_in_barangay_m2=v*pixel_area)
                  for k,v in region_counts.items()]).to_csv(table_dir / 'Analysis_Object_Statistics.csv', index=False)
    connected_region_report(raster_dir / 'Enhanced_ITSA_Recovery_Classification.tif', table_dir, c, pixel_area)
    (table_dir / 'metadata.json').write_text(json.dumps(dict(name=record.name,code=record.code,
       owned_pixels=inside_count, coverage_denominator='owned pixel centers within reference extent; excludes off-raster polygon area',
       missing_support='below 80% excluded from aggregate fit, never filled',
       warning='Raster clips retain whole-study model values. Regions crossing a barangay are counted by their intersection only.'), indent=2), encoding='utf-8')


def graph_series(table, result, obs, folder, c, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    dates = pd.to_datetime(table.reference)
    days = day_offsets(obs,c.EVENT_ORIGIN)
    origin = pd.Timestamp(c.EVENT_ORIGIN)
    beta = result.fit.beta[:,0]
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
        fig, ax = plt.subplots(figsize=(10,5.5), layout='constrained')
        for phase,marker in (('pre','o'),('event','D'),('post','s')):
            use = (table.phase == phase) & table.eligible
            ax.errorbar(dates[use], table['median'][use], yerr=np.vstack((table['median'][use]-table.q25[use],
                        table.q75[use]-table['median'][use])), fmt=marker, color='black', capsize=2, label=phase.upper()+' median / spatial IQR')
        use = table.eligible
        ax.plot(dates[use], table['mean'][use], '+',color='.5',label='Available-pixel mean')
        ax.plot(dates[~use], table['median'][~use], 'x',color='.65',label='Below 80% coverage; excluded from fit')
        if np.isfinite(beta).all():
            for phase in ('pre','post'):
                t = days[np.array([o.phase==phase for o in obs])]
                y = beta[0]+beta[1]*t+(beta[2]+beta[3]*t if phase=='post' else 0)
                ax.plot(origin+pd.to_timedelta(t,unit='D'),y,color='black',label=phase.upper()+' aggregate ITSA fit')
            post = days[days >= 0]
            ax.plot(origin+pd.to_timedelta(post,unit='D'), beta[0]+beta[1]*post,'--',color='.5',label='PRE counterfactual')
        baseline = float(result.metrics['baseline'][0])
        if np.isfinite(baseline):
            ax.axhline(baseline,color='.4',ls=':',label='Fixed PRE median')
            ax.axhline(baseline-c.RECOVERY_TOLERANCE,color='.5',ls='-.',label='Return tolerance')
        ax.axvline(origin,color='black',ls='--',label='25 September 2022 intervention')
        ax.set(title=title, xlabel='Acquisition-pair reference date (UTC)',ylabel='Interferometric coherence')
        ax.set_xlim(min(o.start for o in obs), max(o.end for o in obs))
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
        ax.legend(fontsize=7,loc='upper center',bbox_to_anchor=(.5,-.18),ncol=3,frameon=False)
        ax.grid(axis='y',color='.9')
        fig.text(.01,-.20,'Separate descriptive fit to regional medians; changing spatial support. Spatial IQR is not a regression confidence interval.\n'
                 'Event-straddling observation is displayed and supports detection; excluded from pure-phase regression.',fontsize=8)
        for ext in c.FIGURE_FORMATS:
            for name in ('Coherence_Time_Series','ITSA_Fitted_Trends','Pre_Post_Recovery_Trends'):
                fig.savefig(folder / (name+'.'+ext),dpi=c.FIGURE_DPI,bbox_inches='tight')
        plt.close(fig)
