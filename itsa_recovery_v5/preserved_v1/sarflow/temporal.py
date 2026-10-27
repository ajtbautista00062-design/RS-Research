"""Filename date discovery, phase assignment, ordering and duplicate protection."""
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
import hashlib
import re
import pandas as pd
import numpy as np

PAIR_PATTERN = re.compile(r"(?<!\d)(\d{8})_(\d{8})(?!\d)")


@dataclass(frozen=True)
class Observation:
    name: str
    path: Path
    phase: str
    master: pd.Timestamp  # filename first date; master/slave roles not independently verified
    slave: pd.Timestamp
    start: pd.Timestamp
    end: pd.Timestamp
    reference: pd.Timestamp


def local_path(value: str | Path, root: Path) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"Path outside Simulation workspace: {path}")
    return path


def study_path(value, root):
    path = Path(value).expanduser().resolve()
    designated = Path(r"C:\Users\Arlan\OneDrive\Desktop\Jhonjhon\Remote Sensing Research\Codes\Boundary\INFANTA_BOUNDARY.kml").resolve()
    if path != designated:
        return local_path(path, root)
    return path  # expressly authorized read-only input; outputs remain under Simulation


def timestamp(value: str) -> pd.Timestamp:
    result = pd.Timestamp(value)
    if pd.isna(result) or result.tzinfo is not None:
        raise ValueError("Use timezone-naive ISO timestamps, interpreted as UTC.")
    return result


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover(c: ModuleType) -> tuple[list[Observation], list[dict]]:
    directory = local_path(c.INPUT_DIRECTORY, c.PROJECT_ROOT)
    if not directory.is_dir():
        raise ValueError(f"Input directory is missing: {directory}")
    start_event, end_event = timestamp(c.EVENT_START), timestamp(c.EVENT_END)
    origin = timestamp(c.EVENT_ORIGIN)
    if not start_event <= origin <= end_event or c.REFERENCE_DATE not in ("end", "midpoint"):
        raise ValueError("Invalid event window/origin or REFERENCE_DATE.")
    observations, duplicates = [], []
    by_pair, by_hash, by_path = {}, {}, set()
    for candidate in sorted(directory.iterdir()):
        if candidate.suffix.lower() not in (".tif", ".tiff") or not candidate.is_file():
            continue
        path = local_path(candidate, c.PROJECT_ROOT)
        if path in by_path:
            duplicates.append(dict(file=str(candidate), reason="duplicate resolved path"))
            continue
        by_path.add(path)
        matches = list(PAIR_PATTERN.finditer(path.name))
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one YYYYMMDD_YYYYMMDD pair: {path.name}")
        first, second = [pd.to_datetime(v, format="%Y%m%d", errors="raise") for v in matches[0].groups()]
        start, end = sorted((first, second))
        if start == end:
            raise ValueError(f"Identical acquisition dates: {path.name}")
        if end < start_event:
            phase = "pre"
        elif start > end_event:
            phase = "post"
        elif start <= origin <= end:
            phase = "event"
        else:
            raise ValueError(f"Ambiguous phase near event window: {path.name}")
        lower = path.name.lower()
        claimed = "pre" if "pre-event" in lower else "post" if "post-event" in lower else None
        if claimed and claimed != phase:
            raise ValueError(f"Filename phase disagrees with acquisition dates: {path.name}")
        digest = file_digest(path)
        key = (start, end)
        if key in by_pair:
            previous, previous_digest = by_pair[key]
            if digest != previous_digest:
                raise ValueError(f"Different rasters for the same date pair: {previous.name}, {path.name}. Select one comparable product.")
            duplicates.append(dict(file=str(path), kept=str(previous), reason="identical duplicate acquisition pair"))
            continue
        if digest in by_hash:
            raise ValueError(f"Identical raster bytes have different pair dates: {by_hash[digest].name}, {path.name}")
        by_pair[key] = (path, digest)
        by_hash[digest] = path
        ref = end if c.REFERENCE_DATE == "end" else start + (end-start)/2
        name = f"{phase.upper()}_{start:%Y%m%d}_{end:%Y%m%d}"
        observations.append(Observation(name, path, phase, first, second, start, end, ref))
    observations.sort(key=lambda o: (o.reference, o.start, o.path.name))
    phases = [o.phase for o in observations]
    if not all(p in phases for p in ("pre", "event", "post")) or phases.count("event") != 1:
        raise ValueError("Require PRE, POST, and exactly one event-straddling pair; never pool event pairs.")
    from .correction import validate as validate_correction
    validate_correction(c)
    local_path(c.Output_Folder, c.PROJECT_ROOT)
    if c.Study_Area and not study_path(c.Study_Area, c.PROJECT_ROOT).is_file():
        raise ValueError("Study_Area must be an existing polygon KML.")
    if c.WINDOW_SIZE < 16 or c.SAMPLE_PIXELS < 1 or not 300 <= c.FIGURE_DPI <= 600:
        raise ValueError("Require WINDOW_SIZE>=16, SAMPLE_PIXELS>=1 and DPI in 300..600.")
    if c.NODATA != -9999.0:
        raise ValueError("This export convention uses NODATA=-9999.0.")
    if c.ENABLE_BASELINE_CLASSIFICATION:
        from .recovery import validate
        validate(c)
        if (phases.count("pre") < c.CLASS_MIN_PRE
                or phases.count("post") < max(c.CLASS_MIN_POST, c.END_POST_COUNT)):
            raise ValueError("Not enough unique PRE/POST pairs for baseline recovery.")
    return observations, duplicates


def validate(c: ModuleType) -> list[Observation]:
    return discover(c)[0]


def day_offsets(obs: list[Observation], origin: str) -> np.ndarray:
    return np.array([(o.reference-timestamp(origin)).total_seconds()/86400 for o in obs])

