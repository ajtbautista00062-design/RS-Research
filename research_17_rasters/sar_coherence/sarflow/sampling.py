"""Bounded uniform priority sampling; exact means/counts, sampled quantiles."""
from dataclasses import dataclass, field
import numpy as np
import pandas as pd


@dataclass
class Sample:
    limit: int
    seed: int
    values: np.ndarray = field(default_factory=lambda: np.empty((17, 0), dtype="float32"))
    keys: np.ndarray = field(default_factory=lambda: np.empty(0))
    coordinates: np.ndarray = field(default_factory=lambda: np.empty((2, 0), dtype="int64"))
    count: int = 0
    sums: np.ndarray = field(default_factory=lambda: np.zeros(17, dtype="float64"))

    def __post_init__(self):
        self.rng = np.random.default_rng(self.seed)

    def add(self, values: np.ndarray, rows: np.ndarray, cols: np.ndarray) -> None:
        """Inputs must be complete-case values; each pixel visited exactly once."""
        self.count += values.shape[1]
        self.sums += values.sum(axis=1, dtype="float64")
        keys = self.rng.random(values.shape[1])
        combined = np.r_[self.keys, keys]
        indices = np.argpartition(combined, min(self.limit, len(combined))-1)[:self.limit] if len(combined) else []
        self.values = np.concatenate((self.values, values), axis=1)[:, indices]
        self.coordinates = np.concatenate((self.coordinates, np.vstack((rows, cols))), axis=1)[:, indices]
        self.keys = combined[indices]

    def table(self, names: list[str]) -> pd.DataFrame:
        if not self.count:
            return pd.DataFrame()
        q = np.quantile(self.values, [0.25, 0.5, 0.75], axis=1)
        return pd.DataFrame({"observation": names, "complete_case_pixels": self.count,
                             "sample_pixels": self.values.shape[1], "mean_exact": self.sums/self.count,
                             "q25_sample": q[0], "median_sample": q[1], "q75_sample": q[2]})
