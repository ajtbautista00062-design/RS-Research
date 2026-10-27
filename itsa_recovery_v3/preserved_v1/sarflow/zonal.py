"""Streaming descriptive band statistics and optional approved class counts."""
from dataclasses import dataclass, field
import numpy as np
from .metrics import BAND_NAMES, CLASS_LABELS


@dataclass
class ZoneStats:
    name: str
    band_names: tuple = BAND_NAMES
    count: np.ndarray = field(default_factory=lambda: np.zeros(len(BAND_NAMES), dtype="int64"))
    total: np.ndarray = field(default_factory=lambda: np.zeros(len(BAND_NAMES), dtype="float64"))
    squares: np.ndarray = field(default_factory=lambda: np.zeros(len(BAND_NAMES), dtype="float64"))
    low: np.ndarray = field(default_factory=lambda: np.full(len(BAND_NAMES), np.inf))
    high: np.ndarray = field(default_factory=lambda: np.full(len(BAND_NAMES), -np.inf))
    class_counts: np.ndarray = field(default_factory=lambda: np.zeros(6, dtype="int64"))

    def __post_init__(self):
        n = len(self.band_names)
        self.count = np.zeros(n, dtype="int64")
        self.total = np.zeros(n)
        self.squares = np.zeros(n)
        self.low = np.full(n, np.inf)
        self.high = np.full(n, -np.inf)

    def add(self, bands, mask, classes=None):
        for j, band in enumerate(bands):
            values = band[mask & np.isfinite(band)].astype("float64")
            if values.size:
                self.count[j] += values.size
                self.total[j] += values.sum()
                self.squares[j] += np.square(values).sum()
                self.low[j] = min(self.low[j], values.min())
                self.high[j] = max(self.high[j], values.max())
        if classes is not None:
            self.class_counts += np.bincount(classes[mask], minlength=6)[:6]

    def rows(self):
        result = []
        for j, name in enumerate(self.band_names):
            n = int(self.count[j])
            mean = self.total[j] / n if n else np.nan
            std = np.sqrt(max(0., self.squares[j]/n-mean**2)) if n else np.nan
            result.append(dict(zone=self.name, band=j+1, measurement=name, valid_pixel_count=n,
                               mean=mean, std=std, minimum=self.low[j] if n else np.nan,
                               maximum=self.high[j] if n else np.nan))
        return result

    def classes_row(self):
        assessed = int(self.class_counts[1:5].sum())
        eligible = int(self.class_counts[1:4].sum())
        inside_grid = int(self.class_counts.sum())
        row = dict(zone=self.name, nodata_count=int(self.class_counts[0]),
                   indeterminate_count=int(self.class_counts[5]),
                   assessed_count=assessed, recovery_eligible_count=eligible,
                   recovered_share_eligible=(self.class_counts[1]/eligible if eligible else np.nan),
                   assessed_share_inside_grid=(assessed/inside_grid if inside_grid else np.nan))
        for i, label in CLASS_LABELS.items():
            row[f"class_{i}_count"] = int(self.class_counts[i])
            if i <= 4:
                row[f"class_{i}_share_assessed"] = self.class_counts[i]/assessed if assessed else np.nan
        return row
