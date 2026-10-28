"""Verified external land-cover strata and actual named KML barangay polygons."""
from dataclasses import dataclass
import json
import logging
import re
import numpy as np
import rasterio
from defusedxml import ElementTree as ET
from pyproj import Transformer
from shapely.geometry import Polygon, mapping, shape
from shapely.ops import transform, unary_union
from .spatial import mask_for, read_zones
from .temporal import local_path


@dataclass(frozen=True)
class Barangay:
    name: str
    code: str
    folder: str
    geometry: dict


def safe_name(name):
    text = re.sub(r'[<>:"/\\|?*\x00-\x1f]', '_', name).strip(' .')[:100]
    if not text or re.fullmatch(r'(?i:CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\..*)?', text):
        text = '_' + (text or 'unnamed')
    return text


def load_barangays(path, crs, c):
    if not path or not c.ENABLE_BARANGAY_OUTPUTS:
        return [], dict(status="disabled", valid_features=0)
    tree = ET.parse(local_path(path, c.PROJECT_ROOT))
    project = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform
    records, problems = [], []
    for i, pm in enumerate(tree.getroot().findall('.//{*}Placemark')):
        attributes = {e.get('name'): (e.text or '').strip() for e in pm.findall('.//{*}SimpleData')}
        attributes.update({e.get('name'): e.findtext('{*}value', '').strip() for e in pm.findall('.//{*}Data')})
        name, code = attributes.get(c.BARANGAY_NAME_FIELD), attributes.get(c.BARANGAY_CODE_FIELD)
        if not name or not code:
            problems.append(f"feature {i+1}: lacks explicit barangay name/code fields")
            continue
        polygons = []
        for element in pm.findall('.//{*}Polygon'):
            def ring(boundary):
                return [tuple(map(float, t.split(',')[:2])) for t in boundary.findtext('.//{*}coordinates', '').split()]
            outer = element.find('{*}outerBoundaryIs')
            if outer is None:
                continue
            poly = Polygon(ring(outer), [ring(b) for b in element.findall('{*}innerBoundaryIs')])
            if poly.is_empty or not poly.is_valid:
                problems.append(f"{name}: invalid polygon; not repaired automatically")
                polygons = []
                break
            polygons.append(poly)
        if polygons:
            projected = transform(project, unary_union(polygons))
            if not projected.is_valid:
                problems.append(f"{name}: invalid projected geometry")
                continue
            records.append(Barangay(name, code, safe_name(name)+'_'+safe_name(code), mapping(projected)))
    if len({b.code for b in records}) != len(records):
        raise ValueError("Duplicate barangay codes: dissolve/verify the boundary before analysis.")
    # Stable ownership of shared centers; overlaps are reported, never double counted.
    records.sort(key=lambda b: b.code)
    overlaps = []
    for i, a in enumerate(records):
        for b in records[i+1:]:
            area = shape(a.geometry).intersection(shape(b.geometry)).area
            if area > 0:
                overlaps.append(dict(left=a.code, right=b.code, overlap_crs_area=area))
    report = dict(status="available" if records else "no valid barangay layer; supply one",
                  name_field=c.BARANGAY_NAME_FIELD, code_field=c.BARANGAY_CODE_FIELD,
                  valid_features=len(records), skipped=problems, overlaps=overlaps,
                  ownership="pixel center; first sorted ADM4 code wins shared centers")
    if problems or overlaps:
        logging.warning("Barangay boundary audit: %s", report)
    return records, report


def owners(records, grid_shape, affine):
    result = np.zeros(grid_shape, dtype='uint16')
    for i, record in enumerate(records, 1):
        inside = mask_for([record.geometry], grid_shape, affine, False)
        result[inside & (result == 0)] = i
    return result


class BuiltupLayer:
    """Binary external mask; unknown stays -1, known non-built-up=0, built-up=1."""
    def __init__(self, context, ref, c):
        self.dataset, self.geometries = None, None
        if not c.BUILTUP_MASK:
            logging.warning("No verified built-up mask: land-cover uncertainty retained; no urban inference.")
            return
        if not c.BUILTUP_MASK_NOTES.strip():
            raise ValueError("BUILTUP_MASK requires provenance in BUILTUP_MASK_NOTES.")
        path = local_path(c.BUILTUP_MASK, c.PROJECT_ROOT)
        if path.suffix.lower() in ('.tif', '.tiff'):
            self.dataset = context.enter_context(rasterio.open(path))
            ds = self.dataset
            if (ds.crs, ds.transform, ds.shape, ds.count) != (ref.crs, ref.transform, ref.shape, 1):
                raise ValueError("Built-up raster must match the reference grid exactly and have one band.")
        elif path.suffix.lower() == '.kml':
            self.geometries = [g for _, g in read_zones(path, ref.crs, c.PROJECT_ROOT)]
        elif path.suffix.lower() in ('.json', '.geojson'):
            data = json.loads(path.read_text(encoding='utf-8'))
            features = data.get('features', [])
            project = Transformer.from_crs(c.BUILTUP_GEOJSON_CRS, ref.crs, always_xy=True).transform
            self.geometries = []
            for item in features:
                geom = shape(item['geometry'])
                if geom.geom_type not in ('Polygon', 'MultiPolygon') or not geom.is_valid:
                    raise ValueError("Built-up GeoJSON requires valid polygons.")
                self.geometries.append(mapping(transform(project, geom)))
            if not self.geometries:
                raise ValueError("Empty built-up polygon layer.")
        else:
            raise ValueError("Built-up mask must be a binary TIFF or polygon KML/GeoJSON.")

    def read(self, window, affine):
        size = (int(window.height), int(window.width))
        if self.dataset:
            values = self.dataset.read(1, window=window, masked=True)
            if not np.isin(values.compressed(), (0, 1)).all():
                raise ValueError("Built-up raster must contain only valid 0/1 codes.")
            return values.astype('int16').filled(-1)
        if self.geometries is not None:
            return mask_for(self.geometries, size, affine, False).astype('int16')
        return np.full(size, -1, dtype='int16')


def combine_strata(barriers, builtup, enabled=True):
    # No invented land-cover identity. Unknown, built-up, non-built-up cannot mix.
    if not enabled:
        return barriers.astype('uint64')
    return np.where(barriers > 0, barriers.astype('uint64')*3 + builtup + 1, 0).astype('uint64')
