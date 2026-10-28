"""Three output roots; one shared run identifier prevents overwriting past results."""
from dataclasses import dataclass
from pathlib import Path
from .temporal import local_path


@dataclass(frozen=True)
class OutputLayout:
    rasters: Path
    graphs: Path
    tables: Path


def create_layout(c, run_id):
    root = local_path(c.Output_Folder, c.PROJECT_ROOT)
    paths = [root / name for name in ("Rasters", "ITSA Graphs", "CSV files")]
    if getattr(c, "OUTPUT_VERSIONED_RUNS", True):
        paths = [p / run_id for p in paths]
    # Preflight all locations before creating any: never replace old products.
    for path in paths:
        if path.exists() and any(path.iterdir()):
            raise FileExistsError(f"Output directory is not empty: {path}; enable versioned runs.")
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)
    return OutputLayout(*paths)
