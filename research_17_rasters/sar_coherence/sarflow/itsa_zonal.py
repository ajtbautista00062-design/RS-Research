"""Windowed class shares and continuous summaries, with explicit denominators."""
import numpy as np


class ZonalSummary:
    def __init__(self, name):
        self.name = name
        self.counts = np.zeros(6, dtype="int64")
        self.n = 0
        self.total = self.squares = 0.
        self.minimum, self.maximum = np.inf, -np.inf

    def add(self, classes, fraction, mask):
        self.counts += np.bincount(classes[mask], minlength=6)
        selected = mask & np.isin(classes, [1, 2, 3]) & np.isfinite(fraction)
        values = fraction[selected].astype(float)
        if values.size:
            self.n += values.size
            self.total += float(values.sum())
            self.squares += float((values**2).sum())
            self.minimum = min(self.minimum, float(values.min()))
            self.maximum = max(self.maximum, float(values.max()))

    def row(self):
        valid = int(self.counts[1:].sum())
        eligible = int(self.counts[1:4].sum())
        result = dict(zone=self.name, study_pixels=int(self.counts.sum()),
                      valid_pixels=valid, eligible_pixels=eligible, fraction_valid_pixels=self.n,
                      recovered_among_loss_eligible=self.counts[1]/eligible if eligible else np.nan,
                      mean_fraction=self.total/self.n if self.n else np.nan,
                      std_fraction=np.sqrt(max(0, self.squares/self.n-(self.total/self.n)**2)) if self.n else np.nan,
                      minimum_fraction=self.minimum if self.n else np.nan,
                      maximum_fraction=self.maximum if self.n else np.nan)
        for code, n in enumerate(self.counts):
            result[f"class_{code}_count"] = int(n)
            result[f"class_{code}_share_valid"] = n/valid if code and valid else np.nan
        return result
