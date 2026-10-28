"""Pair-aware configuration validation."""
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
import pandas as pd
import numpy as np
import re


@dataclass(frozen=True)
class Observation:
    name: str
    path: Path
    phase: str
    master: pd.Timestamp
    slave: pd.Timestamp
    start: pd.Timestamp
    end: pd.Timestamp
    reference: pd.Timestamp


def local_path(value: str | Path, root: Path) -> Path:
    """Reject access outside the designated workspace, including resolved symlinks."""
    path = Path(value).expanduser().resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"Path outside PROJECT_ROOT: {path}")
    return path


def timestamp(value: str) -> pd.Timestamp:
    result = pd.Timestamp(value)
    if pd.isna(result) or result.tzinfo is not None:
        raise ValueError("Use explicit, timezone-naive ISO dates/timestamps, all in UTC.")
    return result


def configured_pair(c, name, path):
    """Explicit master/slave roles take precedence; blank pair dates use filename support."""
    dates = c.PAIR_DATES[name]
    if all(dates):
        return tuple(map(timestamp, dates))
    if any(dates):
        raise ValueError("Supply both master/slave dates, or leave both blank for filename parsing.")
    matches = re.findall(r"(?<!\d)(\d{8})_(\d{8})(?!\d)", path.stem)
    if len(matches) != 1:
        raise ValueError("Blank pair dates require one unambiguous YYYYMMDD_YYYYMMDD filename pair.")
    # Filename date order is recorded as supplied, not proof of SLC master/slave roles.
    return tuple(pd.to_datetime(d, format="%Y%m%d", errors="raise") for d in matches[0])


def validate(c: ModuleType) -> list[Observation]:
    expected = [*(f"Pre{i}" for i in range(1, 9)), "During",
                *(f"Post{i}" for i in range(1, 9))]
    if list(c.PATHS) != expected or set(c.PAIR_DATES) != set(expected):
        raise ValueError("Exactly Pre1..Pre8, During, Post1..Post8 are required in order.")
    errors: list[str] = []
    observations: list[Observation] = []
    event_start, event_end = timestamp(c.EVENT_START), timestamp(c.EVENT_END)
    if not event_start <= timestamp(c.EVENT_ORIGIN) <= event_end:
        errors.append("EVENT_ORIGIN must lie inside EVENT_START..EVENT_END.")
    if c.REFERENCE_DATE not in ("end", "midpoint"):
        errors.append("REFERENCE_DATE must be end or midpoint.")
    for name in expected:
        try:
            if not c.PATHS[name]:
                raise ValueError("raster path is empty")
            path = local_path(c.PATHS[name], c.PROJECT_ROOT)
            if not path.is_file():
                raise ValueError(f"raster does not exist: {path}")
            master, slave = configured_pair(c, name, path)
            start, end = sorted((master, slave))
            if start == end:
                raise ValueError("pair acquisitions must differ")
            phase = "pre" if name.startswith("Pre") else "post" if name.startswith("Post") else "event"
            if phase == "pre" and end >= event_start:
                raise ValueError("pre pair must finish before event window")
            if phase == "post" and start <= event_end:
                raise ValueError("post pair must start after event window")
            if phase == "event" and (end < event_start or start > event_end):
                raise ValueError("During pair must overlap the event window")
            ref = end if c.REFERENCE_DATE == "end" else start + (end - start) / 2
            observations.append(Observation(name, path, phase, master, slave, start, end, ref))
        except (ValueError, TypeError, KeyError) as exc:
            errors.append(f"{name}: {exc}")
    if len(observations) == 17:
        if len({o.path for o in observations}) != 17:
            errors.append("Duplicate input paths are not allowed.")
        if len({(o.start, o.end) for o in observations}) != 17:
            errors.append("Duplicate acquisition pairs are not allowed.")
        if getattr(c, "RECOVERY_MODEL", "legacy_observed") == "itsa":
            observations.sort(key=lambda o: o.reference)
        if not all(a.reference < b.reference and a.end < b.end
                   for a, b in zip(observations, observations[1:])):
            errors.append("Observations must have distinct increasing reference AND pair-end dates.")
        spans = [(o.end - o.start).total_seconds() / 86400 for o in observations]
        if max(spans) / min(spans) > c.MAX_PAIR_BASELINE_RATIO and not c.ALLOW_BASELINE_MISMATCH:
            errors.append("Pair temporal baselines differ too much: harmonize pairs or justify override.")
    if not c.ACQUISITION_COMPARABILITY_CONFIRMED or not c.COMPARABILITY_NOTES.strip():
        errors.append("Confirm acquisition comparability and supply COMPARABILITY_NOTES.")
    if not c.Output_Folder:
        errors.append("Output_Folder is empty.")
    else:
        local_path(c.Output_Folder, c.PROJECT_ROOT)
    if c.Study_Area:
        study = local_path(c.Study_Area, c.PROJECT_ROOT)
        if not study.is_file() or study.suffix.lower() != ".kml":
            errors.append("Study_Area must be an existing polygon .kml.")
    if not (1 <= c.MIN_PRE <= 8 and 2 <= c.CONFIRMATIONS <= c.MIN_POST <= 8):
        errors.append("Require 1<=MIN_PRE<=8 and 2<=CONFIRMATIONS<=MIN_POST<=8.")
    if (getattr(c, "RECOVERY_MODEL", "legacy_observed") == "legacy_observed"
            and not (0 < c.RECOVERY_TOLERANCE < c.LOSS_ABSOLUTE <= 1)):
        errors.append("Require 0 < recovery tolerance < minimum detectable event loss <= 1.")
    if not 0 < c.RECOVERY_TOLERANCE < 1:
        errors.append("Recovery tolerance must lie strictly between zero and one.")
    if not (2 <= c.RECOVERY_CONFIRMATIONS <= c.ENDSTATE_POST_PAIRS <= 8):
        errors.append("Require 2<=RECOVERY_CONFIRMATIONS<=ENDSTATE_POST_PAIRS<=8.")
    if not (0 <= c.BASELINE_MIN <= 1 and 0 <= c.BASELINE_MAD_MAX <= 1
            and c.NOISE_FLOOR > 0 and c.LOSS_MAD_MULTIPLIER >= 0):
        errors.append("Invalid baseline or noise settings.")
    if c.WINDOW_SIZE < 16 or c.SAMPLE_PIXELS < 1 or c.FIGURE_DPI < 300:
        errors.append("WINDOW_SIZE>=16, SAMPLE_PIXELS>=1, FIGURE_DPI>=300 required.")
    if c.MAX_CONFIRM_GAP_DAYS <= 0 or c.MIN_CONFIRM_SPAN_DAYS < 0:
        errors.append("Invalid confirmation gap/span settings.")
    if errors:
        raise ValueError("Configuration errors:\n - " + "\n - ".join(errors))
    return observations


def day_offsets(obs: list[Observation], origin: str) -> np.ndarray:
    return np.array([(o.reference - timestamp(origin)).total_seconds() / 86400 for o in obs])
