"""Streaming diagnostics: exact counts, bounded descriptive histograms, sensitivity."""
from collections import Counter
import csv
import numpy as np
from .recovery import assess, REASONS


def write_records(writer, ids, records, pixel_area, obs, stage):
    for gid, (result, series, counts, fractions, spread) in records.items():
        row = dict(stage=stage, region_id=gid, area_m2=int((ids == gid).sum())*pixel_area,
                   recovery_class=int(result.classes[0]), quality_flags=int(result.quality[0]),
                   reason=int(result.reasons[0]), **{k: float(v[0]) for k, v in result.metrics.items()})
        for i, o in enumerate(obs):
            row.update({f"{o.name}_median": float(series[i]), f"{o.name}_valid_count": int(counts[i]),
                        f"{o.name}_valid_fraction": float(fractions[i]), f"{o.name}_mad_sigma": float(spread[i])})
        if writer[0] is None:
            writer[0] = csv.DictWriter(writer[1], fieldnames=list(row))
            writer[0].writeheader()
        writer[0].writerow(row)


class Diagnostics:
    def __init__(self):
        self.counts = {}
        self.reasons = Counter()
        self.quality = Counter()
        self.sensitivity = {}
        self.transitions = {}
        self.histograms = {}
        self.gates = Counter()
        self.evidence = Counter()

    def add(self, products, result, mask):
        if result.gates:
            for name, gate in result.gates.items():
                self.gates[(name, True)] += int(gate[mask].sum())
                self.gates[(name, False)] += int((~gate[mask]).sum())
        decided = np.isin(result.classes, [1, 2, 3, 4]) & mask
        self.evidence["provisional"] += int((decided & (result.quality != 0)).sum())
        self.evidence["no_quality_warning"] += int((decided & (result.quality == 0)).sum())
        self.evidence["indeterminate"] += int(((result.classes == 5) & mask).sum())
        for name, classes in products.items():
            self.counts.setdefault(name, np.zeros(6, dtype="int64"))
            self.counts[name] += np.bincount(classes[mask], minlength=6)
        for code, count in zip(*np.unique(result.reasons[mask], return_counts=True)):
            self.reasons[int(code)] += int(count)
        for bit in (1, 2, 4, 8, 16, 32):
            self.quality[bit] += int(((result.quality[mask] & bit) > 0).sum())
        for name in ("modeled_level_loss", "post_slope", "baseline_deficit", "modeled_progress"):
            a = result.metrics[name][mask]
            # Explicit under/overflow bins; no clipping or sampled distributions.
            inner = np.linspace(-.05, .05, 201) if name == "post_slope" else np.linspace(-1, 1, 201)
            edges = np.r_[-np.inf, inner, np.inf]
            hist = np.histogram(a[np.isfinite(a)], edges)[0]
            if name not in self.histograms:
                self.histograms[name] = [edges, hist]
            else:
                self.histograms[name][1] += hist
        if "pixel" in products and "object" in products:
            pair = products["pixel"][mask]*6+products["object"][mask]
            self.transitions.setdefault("pixel_to_object", np.zeros(36, dtype="int64"))
            self.transitions["pixel_to_object"] += np.bincount(pair, minlength=36)
        if "legacy" in products:
            pair = products["legacy"][mask]*6+products["pixel"][mask]
            self.transitions.setdefault("legacy_to_pixel", np.zeros(36, dtype="int64"))
            self.transitions["legacy_to_pixel"] += np.bincount(pair, minlength=36)

    def sweep(self, series, obs, c, areas, main_classes):
        if not series.size:
            return
        scenarios = [(loss, c.RECOVERY_TOLERANCE) for loss in c.ITSA_LOSS_SENSITIVITY]
        scenarios += [(c.EVENT_LOSS_MINIMUM, tol) for tol in c.SENSITIVITY_TOLERANCES]
        for loss, tol in sorted(set(scenarios)):
            a = assess(series, obs, c, loss_floor=loss, tolerance=tol)
            counts = np.bincount(a.classes, weights=areas, minlength=6)
            key = (loss, tol)
            self.sensitivity.setdefault(key, [np.zeros(6), 0., 0.])
            self.sensitivity[key][0] += counts
            self.sensitivity[key][1] += float(areas[a.classes != main_classes].sum())
            self.sensitivity[key][2] += float(areas.sum())

    def export(self, out, pixel_area):
        import pandas as pd
        pd.DataFrame([dict(gate=k[0], accepted=k[1], pixels=v) for k, v in self.gates.items()]).to_csv(
            out / "classification_gates.csv", index=False)
        pd.DataFrame([dict(evidence=k, pixels=v) for k, v in self.evidence.items()]).to_csv(
            out / "operational_evidence_counts.csv", index=False)
        pd.DataFrame([dict(reason=k, label=REASONS[k], pixels=v) for k, v in self.reasons.items()]).to_csv(
            out / "exclusion_reasons.csv", index=False)
        pd.DataFrame([dict(quality_bit=k, pixels=v) for k, v in self.quality.items()]).to_csv(
            out / "quality_flags.csv", index=False)
        rows = []
        for (loss, tol), (counts, changed, total) in self.sensitivity.items():
            for code, n in enumerate(counts):
                rows.append(dict(loss_floor=loss, return_tolerance=tol, code=code, pixels=n,
                                 area_m2=n*pixel_area, changed_share=changed/total if total else None,
                                 denominator="PRE-segmented object members only; fixed geometry and cohorts"))
        pd.DataFrame(rows).to_csv(out / "threshold_sensitivity.csv", index=False)
        for name, counts in self.transitions.items():
            pd.DataFrame(counts.reshape(6, 6), index=range(6), columns=range(6)).to_csv(out / f"{name}.csv")
        for name, (edges, hist) in self.histograms.items():
            pd.DataFrame({"lower": edges[:-1], "upper": edges[1:], "pixels": hist}).to_csv(
                out / f"distribution_{name}.csv", index=False)
