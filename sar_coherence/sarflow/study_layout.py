"""Requested RRI and ITSA hierarchy, shared immutable versioned run."""
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ProductFolders:
    rasters: Path
    graphs: Path
    csv: Path


def folders(root: Path, method: str, barangay: str | None = None) -> ProductFolders:
    level = Path('Overall') if barangay is None else Path('Per Barangay') / barangay
    if method == 'RRI':
        paths = [root / 'A. RRI' / category / level for category in ('Raster Files', 'Graphs', 'CSV')]
    elif method == 'ITSA':
        paths = [root / 'B. ITSA' / level / category for category in ('Raster Files', 'ITSA Graph', 'CSV')]
    else:
        raise ValueError(method)
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)
    return ProductFolders(*paths)
