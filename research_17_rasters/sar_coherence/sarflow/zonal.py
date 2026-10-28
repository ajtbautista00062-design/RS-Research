"""Windowed class counts and continuous recovery-fraction zonal summaries."""
from dataclasses import dataclass, field
import numpy as np

CLASS_LABELS = {1: "recovered", 2: "partially_recovered", 3: "not_recovered",
                4: "no_detectable_event_loss"}


@dataclass
class ZoneStats:
    name: str
    counts: np.ndarray = field(default_factory=lambda: np.zeros(5, dtype="int64"))
    fractions: list[np.ndarray] = field(default_factory=list)

    def add(self, fractions: np.ndarray, classes: np.ndarray, mask: np.ndarray) -> None:
        values = classes[mask]
        valid_codes = values[(values == -9999) | ((values >= 1) & (values <= 4))]
        if valid_codes.size:
            # Index 0 records NoData; indices 1..4 are the class codes.
            self.counts += np.bincount(np.where(valid_codes == -9999, 0, valid_codes), minlength=5)[:5]
        eligible = mask & np.isin(classes, [1, 2, 3]) & np.isfinite(fractions)
        if eligible.any():
            self.fractions.append(fractions[eligible].astype("float32", copy=True))

    def summary(self) -> dict:
        total = int(self.counts.sum())
        eligible = int(self.counts[1:4].sum())
        row = {"zone": self.name, "zone_pixel_count": total,
               "nodata_pixel_count": int(self.counts[0])}
        for code, label in CLASS_LABELS.items():
            count = int(self.counts[code])
            row[f"class_{code}_{label}_count"] = count
            row[f"class_{code}_{label}_share_of_zone"] = count / total if total else np.nan
        row["recovery_eligible_pixel_count"] = eligible
        row["recovered_fraction_among_eligible"] = self.counts[1] / eligible if eligible else np.nan
        if not self.fractions:
            row.update({key: np.nan for key in ("recovery_fraction_mean", "recovery_fraction_median",
                "recovery_fraction_std", "recovery_fraction_q25", "recovery_fraction_q75",
                "recovery_fraction_min", "recovery_fraction_max")})
            return row
        values = np.concatenate(self.fractions)
        row.update(recovery_fraction_mean=float(np.mean(values)),
                   recovery_fraction_median=float(np.median(values)),
                   recovery_fraction_std=float(np.std(values, ddof=0)),
                   recovery_fraction_q25=float(np.quantile(values, .25)),
                   recovery_fraction_q75=float(np.quantile(values, .75)),
                   recovery_fraction_min=float(np.min(values)),
                   recovery_fraction_max=float(np.max(values)))
        return row
