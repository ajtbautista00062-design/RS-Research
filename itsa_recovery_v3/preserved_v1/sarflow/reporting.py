"""Simulation disclosures and a categorical overview of the GIS deliverable."""
from pathlib import Path
from types import ModuleType
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

from .metrics import CLASS_LABELS
from .figures import save

COLORS = ['#009E73', '#E6AB02', '#D55E00', '#0072B2']


def disclosures(out: Path) -> None:
    rows = [
        ('pre_event_variability', 'One pre pair cannot estimate temporal variability or baseline stability'),
        ('recovery_trend', 'One post pair cannot define a post-event slope'),
        ('recovery_duration', 'No temporal recovery trajectory or confirmation bracket'),
        ('persistent_recovery', 'No second post pair to confirm persistence'),
        ('structural_recovery', 'Coherence alone cannot validate reconstruction'),
        ('acquisition_comparability', 'SLC orbit, polarization, bursts, and processing records not supplied'),
    ]
    pd.DataFrame([dict(analysis=name, status='unavailable', reason=reason)
                  for name, reason in rows]).to_csv(out / 'unavailable_analyses.csv', index=False)


def gis_overview(output: Path, c: ModuleType) -> None:
    """Read categorical overview with nearest-neighbor; source raster is unchanged."""
    with rasterio.open(output) as src:
        scale = min(1., 1100 / max(src.shape))
        image = src.read(2, out_shape=(max(1, int(src.height*scale)),
                                      max(1, int(src.width*scale))),
                         resampling=Resampling.nearest, masked=True)
        bounds = src.bounds
        crs = str(src.crs)
    fig, ax = plt.subplots(figsize=(10, 8))
    cmap = ListedColormap(COLORS)
    cmap.set_bad('#eeeeee')
    ax.imshow(image, cmap=cmap, norm=BoundaryNorm([.5, 1.5, 2.5, 3.5, 4.5], 4),
              extent=(bounds.left, bounds.right, bounds.bottom, bounds.top))
    ax.set(title='Three-raster simulation: single-post coherence response',
           xlabel=f'X coordinate ({crs})', ylabel=f'Y coordinate ({crs})')
    handles = [Patch(color=color, label=f'{i}: {CLASS_LABELS[i]}')
               for i, color in enumerate(COLORS, start=1)]
    handles.append(Patch(color='#eeeeee', label='NoData: invalid, missing or low pre-coherence'))
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, -.12), fontsize=8)
    ax.text(.5, -.35, 'Functionality test only. No trend, duration, persistent\nor structural recovery inference.',
            transform=ax.transAxes, ha='center', va='top', fontsize=9)
    save(fig, output.parent, '08_single_post_class_map', c)
    qgis = ET.Element('qgis', version='3.34', styleCategories='Symbology')
    pipe = ET.SubElement(qgis, 'pipe')
    renderer = ET.SubElement(pipe, 'rasterrenderer', type='paletted', band='2', opacity='1', alphaBand='-1')
    palette = ET.SubElement(renderer, 'colorPalette')
    for code, color in enumerate(COLORS, start=1):
        ET.SubElement(palette, 'paletteEntry', value=str(code), color=color,
                      alpha='255', label=CLASS_LABELS[code])
    ET.ElementTree(qgis).write(output.with_suffix('.qml'), encoding='utf-8', xml_declaration=True)
