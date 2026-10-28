"""Future-run tiled enhancement, transparent transitions and global components."""
from contextlib import ExitStack
import csv
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from .recovery import Assessment, CLASSES
from .recovery_io import METRICS, write
from .spatial import windows
from .spatial_enhancement import refine


class Comparison:
    def __init__(self):
        self.counts = {s: np.zeros((2, 6), dtype='int64') for s in ('all', 'builtup', 'nonbuiltup', 'unknown')}
        self.transitions = np.zeros((6, 6), dtype='int64')
        self.hist = np.zeros(1001, dtype='int64')
        self.score_sum = 0.

    def add(self, original, enhanced, confidence, inside, landcover):
        for name, mask in (('all', inside), ('builtup', inside & (landcover == 1)),
                           ('nonbuiltup', inside & (landcover == 0)), ('unknown', inside & (landcover == -1))):
            for i, data in enumerate((original, enhanced)):
                self.counts[name][i] += np.bincount(data[mask], minlength=6)
        self.transitions += np.bincount((original[inside]*6+enhanced[inside]), minlength=36).reshape(6, 6)
        valid = inside & (enhanced > 0) & np.isfinite(confidence)
        values = confidence[valid]
        self.hist += np.bincount(np.rint(np.clip(values, 0, 1)*1000).astype(int), minlength=1001)
        self.score_sum += float(values.sum(dtype='float64'))

    def export(self, folder, pixel_area, identity='study'):
        rows = []
        for stratum, counts in self.counts.items():
            for mode, values in zip(('original', 'enhanced'), counts):
                for code, n in enumerate(values):
                    rows.append(dict(zone=identity, stratum=stratum, product=mode, code=code,
                                     label=CLASSES[code], pixels=n, area_m2=n*pixel_area,
                                     percent_valid=100*n/values[1:].sum() if code and values[1:].sum() else np.nan,
                                     denominator='valid classes 1-5 in this stratum; NoData excluded'))
        frame = pd.DataFrame(rows)
        frame.to_csv(folder / 'Recovery_Class_Area_Statistics.csv', index=False)
        frame.to_csv(folder / 'Raw_vs_Enhanced_Comparison.csv', index=False)
        pd.DataFrame([dict(original=a, enhanced=b, pixels=self.transitions[a, b],
                           area_m2=self.transitions[a, b]*pixel_area)
                      for a in range(6) for b in range(6)]).to_csv(folder / 'Class_Transition_Matrix.csv', index=False)
        n = self.hist.sum()
        median = np.searchsorted(np.cumsum(self.hist), (n+1)/2)/1000 if n else np.nan
        pd.DataFrame([dict(zone=identity, count=n, mean=self.score_sum/n if n else np.nan,
                           median_rounded_0_001=median, probability_calibrated=False)]).to_csv(
                               folder / 'Evidence_Score_Summary.csv', index=False)


def render_comparison(summary, folder, c):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .recovery_io import COLORS
    values = summary.counts['all'][:,1:]
    with plt.rc_context({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False}):
        fig, ax = plt.subplots(figsize=(9,4.5), layout='constrained')
        x = np.arange(5)
        for i, label in enumerate(('Original','Enhanced')):
            total = values[i].sum()
            ax.bar(x+(i-.5)*.36, 100*values[i]/max(total,1), .36,
                   color=[np.array(COLORS[k][:3])/255 for k in range(1,6)],
                   edgecolor='black', linewidth=.5, hatch='//' if i == 0 else None, label=label)
        ax.set_xticks(x, ['Recovered','Partial','Not recovered','No loss','Indeterminate'])
        ax.set_ylabel('Share of valid classified pixels (%)')
        ax.set_title('Original versus enhanced ITSA recovery; coverage is not accuracy')
        ax.legend(frameon=False)
        for ext in c.FIGURE_FORMATS:
            fig.savefig(folder / ('Recovery_Class_Comparison.'+ext), dpi=c.FIGURE_DPI)
        plt.close(fig)


def read_metrics(dataset, win):
    data = dataset.read(window=win, masked=True).filled(np.nan)
    return dict(zip(METRICS, data))


def enhance_outputs(paths, layout, c):
    summary = Comparison()
    with ExitStack() as context:
        keys = ('Original_ITSA_Recovery_Classification', 'Original_Evidence_Confidence',
                'Raw_Pixel_ITSA_Classes', 'Raw_ITSA_Metrics', 'Raw_ITSA_Reasons',
                'Spatial_Strata', 'Builtup_Status', 'ITSA_Metrics', 'Raw_ITSA_Quality_Flags',
                'ITSA_Quality_Flags', 'ITSA_Exclusion_Reasons')
        readers = {key: context.enter_context(rasterio.open(paths[key])) for key in keys}
        output_keys = ('Enhanced_ITSA_Recovery_Classification', 'Evidence_Confidence',
                       'Classification_Changes', 'Enhanced_ITSA_Metrics',
                       'ITSA_Disturbance_Magnitude', 'ITSA_Recovery_Magnitude', 'ITSA_Recovery_Fraction',
                       'Enhanced_ITSA_Quality_Flags', 'Enhanced_ITSA_Reasons')
        writers = {key: context.enter_context(rasterio.open(paths[key], 'r+')) for key in output_keys}
        ref = readers[keys[0]]
        factor = ref.crs.linear_units_factor[1]
        size = (abs(ref.transform.a)*factor, abs(ref.transform.e)*factor)
        area = size[0]*size[1]
        halo = int(np.ceil(c.NEIGHBOR_RADIUS_M/min(size)))+1
        log = context.enter_context((layout.tables / 'Pixel_Reassignments.csv').open('w', newline='', encoding='utf-8'))
        detail = csv.writer(log)
        detail.writerow(('row', 'column', 'original_class', 'enhanced_class', 'evidence_score', 'reason', 'threshold_adjusted'))
        for win in windows(ref.width, ref.height, c.WINDOW_SIZE):
            left, top = max(0, int(win.col_off)-halo), max(0, int(win.row_off)-halo)
            right, bottom = min(ref.width, int(win.col_off+win.width)+halo), min(ref.height, int(win.row_off+win.height)+halo)
            expanded = Window(left, top, right-left, bottom-top)
            core = np.s_[int(win.row_off)-top:int(win.row_off-top+win.height),
                         int(win.col_off)-left:int(win.col_off-left+win.width)]
            original = readers[keys[0]].read(1, window=expanded)
            confidence = readers[keys[1]].read(1, window=expanded, masked=True).filled(0)
            metrics = read_metrics(readers['Raw_ITSA_Metrics'], expanded)
            raw = Assessment(readers['Raw_Pixel_ITSA_Classes'].read(1, window=expanded),
                             np.zeros(original.shape, 'uint16'), readers['Raw_ITSA_Reasons'].read(1, window=expanded), metrics, None)
            strata = readers['Spatial_Strata'].read(1, window=expanded)
            enhanced, score, reason = refine(original, confidence, raw, strata, c, size)
            if not c.EXPORT_SPATIALLY_ENHANCED_CLASSIFICATION:
                enhanced, score = original.copy(), confidence.copy()
                reason = np.where(original > 0, 0, 255).astype('uint8')
            before, after, score, reason = original[core], enhanced[core], score[core], reason[core]
            changed = before != after
            # Displayed continuous indicators follow each changed pixel's own model.
            primary_metrics = read_metrics(readers['ITSA_Metrics'], win)
            final_metrics = {k: np.where(changed, metrics[k][core], v) for k, v in primary_metrics.items()}
            write(writers['Enhanced_ITSA_Recovery_Classification'], after, win)
            write(writers['Evidence_Confidence'], np.where(after > 0, score, np.nan), win)
            write(writers['Classification_Changes'], np.stack((np.where(before > 0, before, 255),
                  np.where(after > 0, after, 255), reason)), win)
            write(writers['Enhanced_ITSA_Metrics'], np.stack([final_metrics[k] for k in METRICS]), win)
            for destination, original_key, raw_key in (
                ('Enhanced_ITSA_Quality_Flags','ITSA_Quality_Flags','Raw_ITSA_Quality_Flags'),
                ('Enhanced_ITSA_Reasons','ITSA_Exclusion_Reasons','Raw_ITSA_Reasons')):
                value = np.where(changed, readers[raw_key].read(1, window=win),
                                 readers[original_key].read(1, window=win))
                write(writers[destination], value, win)
            for name, key in (('ITSA_Disturbance_Magnitude', 'modeled_level_loss'),
                              ('ITSA_Recovery_Magnitude', 'modeled_progress'), ('ITSA_Recovery_Fraction', 'recovery_fraction')):
                write(writers[name], final_metrics[key], win)
                writers[name].update_tags(units='dimensionless' if key == 'recovery_fraction' else 'coherence',
                                         model_scope='original analytical model, replaced by own pixel fit where spatial change accepted')
            land = readers['Builtup_Status'].read(1, window=win)
            summary.add(before, after, score, strata[core] > 0, land)
            rr, cc = np.where(changed)
            detail.writerows((int(win.row_off+r), int(win.col_off+s), int(before[r,s]), int(after[r,s]),
                              float(score[r,s]), 'own ITSA + connected spatial evidence', False) for r,s in zip(rr,cc))
    summary.export(layout.tables, area)
    connected_region_report(paths['Enhanced_ITSA_Recovery_Classification'], layout.tables, c, area)
    return summary


def connected_region_report(path, folder, c, pixel_area):
    """Exact four-connected class components across ALL tile/object boundaries.

    Row runs keep only active components in memory. IDs are minimum linear cell
    index+1. Mapping-area flags never delete a small legitimate recovery region.
    """
    with rasterio.open(path) as src, (folder / 'Spatial_Region_Statistics.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(('component_id', 'class', 'pixels', 'area_m2', 'below_mmu', 'near_target_20pct', 'exceeds_object_max'))
        previous, active = [], {}
        totals = np.zeros((6, 4), dtype='int64')
        def emit(gid, value):
            code, count = value
            area = count*pixel_area
            flags = (area < c.MIN_MAPPING_AREA_M2, abs(area-c.TARGET_REGION_AREA_M2) <= .2*c.TARGET_REGION_AREA_M2,
                     area > c.MAX_REGION_AREA_M2)
            writer.writerow((gid, code, count, area, *flags))
            totals[code] += (1, *flags)
        for block in range(0, src.height, 128):
            data = src.read(1, window=Window(0, block, src.width, min(128, src.height-block)))
            for offset, row in enumerate(data):
                y = block+offset
                edges = np.r_[0, np.flatnonzero(row[1:] != row[:-1])+1, src.width]
                current, parent = [], {}
                def root(i):
                    while i in parent:
                        i = parent[i]
                    return i
                start_prev = 0
                for left, right in zip(edges[:-1], edges[1:]):
                    code = int(row[left])
                    if code == 0:
                        continue
                    while start_prev < len(previous) and previous[start_prev][1] <= left:
                        start_prev += 1
                    links = set()
                    j = start_prev
                    while j < len(previous) and previous[j][0] < right:
                        pl, pr, pc, gid = previous[j]
                        if pc == code and pr > left:
                            links.add(root(gid))
                        j += 1
                    gid = min(links) if links else y*src.width+int(left)+1
                    active.setdefault(gid, [code, 0])
                    for other in links-{gid}:
                        active[gid][1] += active.pop(other)[1]
                        parent[other] = gid
                    active[gid][1] += int(right-left)
                    current.append((left, right, code, gid))
                current = [(l,r,k,root(g)) for l,r,k,g in current]
                alive = {g for _,_,_,g in current}
                for gid in list(active):
                    if gid not in alive:
                        emit(gid, active.pop(gid))
                previous = current
        for gid, value in active.items():
            emit(gid, value)
    pd.DataFrame([dict(code=k, connected_regions=v[0], below_mmu=v[1], near_target=v[2], above_object_max=v[3])
                  for k,v in enumerate(totals)]).to_csv(folder / 'Connected_Region_Counts.csv', index=False)
