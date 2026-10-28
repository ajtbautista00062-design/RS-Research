"""KML polygons, source inventory, and exact reference-grid virtual warping."""
from dataclasses import asdict
from pathlib import Path
from types import ModuleType
import json
import logging
import numpy as np
import rasterio
from rasterio.features import geometry_mask
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.windows import Window
from defusedxml import ElementTree as ET
from pyproj import Transformer
from shapely.geometry import Polygon, MultiPolygon, mapping, box
from shapely.ops import transform, unary_union
from .temporal import Observation, local_path


def read_zones(path: str, crs: rasterio.crs.CRS, root: Path) -> list[tuple[str, dict]]:
    """Parse all polygon Placemarks including holes; ignore altitude. No GDAL KML driver needed."""
    if not path:
        return []
    source = local_path(path, root)
    if not source.is_file():
        raise FileNotFoundError(source)
    tree = ET.parse(source)
    project = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform
    zones: list[tuple[str, dict]] = []
    for i, pm in enumerate(tree.getroot().findall(".//{*}Placemark")):
        polygons = []
        for element in pm.findall(".//{*}Polygon"):
            def ring(boundary: ET.Element) -> list[tuple[float, float]]:
                text = boundary.findtext(".//{*}coordinates", "")
                return [tuple(map(float, token.split(",")[:2])) for token in text.split()]
            outer = element.find("{*}outerBoundaryIs")
            if outer is None:
                raise ValueError("KML Polygon missing outerBoundaryIs.")
            poly = Polygon(ring(outer), [ring(b) for b in element.findall("{*}innerBoundaryIs")])
            if poly.is_empty or not poly.is_valid:
                raise ValueError(f"Invalid KML polygon at Placemark {i+1}; repair in QGIS.")
            polygons.append(poly)
        if polygons:
            merged = unary_union(polygons)
            projected = transform(project, merged)
            name = pm.findtext("{*}name", f"zone_{i+1}")
            zones.append((f"{i+1:03d}_{name}", mapping(projected)))
    if not zones:
        raise ValueError("KML has no polygon Placemarks. Points/lines are not study boundaries.")
    if len(zones) > 1:
        from shapely.geometry import shape
        geoms = [shape(g) for _, g in zones]
        if any(a.intersection(b).area > 0 for i, a in enumerate(geoms) for b in geoms[i+1:]):
            logging.warning("Zonal polygons overlap: shared pixels contribute to multiple zones.")
    return zones


def inventory(srcs: list, obs: list[Observation], c: ModuleType) -> list[dict]:
    rows = []
    ref = srcs[0]
    if not ref.crs or not ref.crs.is_projected:
        raise ValueError("Pre1 must use a projected analysis CRS; prepare a projected reference first.")
    if ref.transform.b != 0 or ref.transform.d != 0 or ref.transform.a <= 0 or ref.transform.e >= 0:
        raise ValueError("Pre1 must have a north-up, non-rotated grid.")
    for src, o in zip(srcs, obs):
        if (getattr(c, "ZERO_IS_INVALID", False) and src.nodata != 0
                and c.RESAMPLING != "nearest"):
            raise ValueError("Undeclared zero holes require nearest resampling or an explicitly masked source; bilinear would mix zeros into neighbors.")
        if not src.crs or c.BAND > src.count or c.BAND < 1:
            raise ValueError(f"Missing CRS or invalid band: {o.name}")
        if np.issubdtype(np.dtype(src.dtypes[c.BAND-1]), np.complexfloating):
            raise ValueError(f"{o.name} contains complex SLC data, not coherence magnitude.")
        rows.append({"name": o.name, "path": str(o.path), "crs": str(src.crs),
                     "resolution": src.res, "width": src.width, "height": src.height,
                     "bounds": tuple(src.bounds), "transform": tuple(src.transform),
                     "nodata": src.nodata, "dtype": src.dtypes[c.BAND-1],
                     "scale": src.scales[c.BAND-1], "offset": src.offsets[c.BAND-1],
                     "grid_matches_reference": src.crs == ref.crs and src.transform == ref.transform
                     and src.shape == ref.shape})
    return rows


def aligned(src, ref, c: ModuleType) -> WarpedVRT:
    return WarpedVRT(src, crs=ref.crs, transform=ref.transform, width=ref.width,
                     height=ref.height, dtype="float32", nodata=np.nan,
                     resampling=Resampling[c.RESAMPLING], warp_mem_limit=128)


def windows(width: int, height: int, size: int):
    for row in range(0, height, size):
        for col in range(0, width, size):
            yield Window(col, row, min(size, width-col), min(size, height-row))


def mask_for(geometries: list[dict], shape: tuple[int, int], affine, touched: bool) -> np.ndarray:
    if not geometries:
        return np.ones(shape, dtype=bool)
    return geometry_mask(geometries, out_shape=shape, transform=affine,
                         invert=True, all_touched=touched)


def read_stack(vrts: list, srcs: list, window: Window, c: ModuleType) -> np.ndarray:
    layers = []
    for vrt, src in zip(vrts, srcs):
        arr = vrt.read(c.BAND, window=window, masked=True).filled(np.nan)
        arr = (arr * src.scales[c.BAND-1] + src.offsets[c.BAND-1])
        arr = arr * c.COHERENCE_SCALE + c.COHERENCE_OFFSET
        if getattr(c, "ZERO_IS_INVALID", False):
            # This policy is separately validated against processing records in ITSA mode.
            arr[arr == 0] = np.nan
        # Fail rather than silently truncate integer-scaled or incorrect input values.
        finite = arr[np.isfinite(arr)]
        if finite.size and (finite.min() < -1e-6 or finite.max() > 1 + 1e-6):
            raise ValueError(f"Coherence outside [0,1] in {src.name}; check band, scaling, and NoData.")
        layers.append(np.clip(arr, 0, 1).astype("float32"))
    return np.stack(layers)


def verify_source_values(src, c: ModuleType) -> None:
    """Validate BEFORE warping, so interpolation cannot conceal out-of-range source cells."""
    for window in windows(src.width, src.height, c.WINDOW_SIZE):
        data = src.read(c.BAND, window=window, masked=True).astype("float32").filled(np.nan)
        data = (data * src.scales[c.BAND-1] + src.offsets[c.BAND-1]) * c.COHERENCE_SCALE + c.COHERENCE_OFFSET
        valid = data[np.isfinite(data)]
        if valid.size and (valid.min() < -1e-6 or valid.max() > 1+1e-6):
            raise ValueError(f"Invalid raw coherence range in {src.name}; check scale/offset/NoData.")


def export_profile(ref, c: ModuleType) -> dict:
    return dict(driver="GTiff", crs=ref.crs, transform=ref.transform,
                width=ref.width, height=ref.height, count=2, dtype="float32",
                nodata=c.NODATA, tiled=True, blockxsize=256, blockysize=256,
                compress="deflate", predictor=3, BIGTIFF="IF_SAFER", interleave="band")
