"""Recovery exports and independent reread checks, on the original reference grid."""
import json
import numpy as np
import rasterio
from .recovery import LABELS, COLORS, settings, baseline_recovery
from .spatial import windows

def paths(out, dates):
    return (out / f"Infanta_Baseline_Recovery_Classes_{dates}.partial.tif",
            out / f"Infanta_Baseline_Recovery_Fraction_{dates}.partial.tif")

def open_outputs(resources, profile, out, dates, c):
    class_path, fraction_path = paths(out, dates)
    classes = resources.enter_context(rasterio.open(
        class_path, "w", **(profile | dict(count=1, dtype="uint8", nodata=0, predictor=2))))
    fraction = resources.enter_context(rasterio.open(
        fraction_path, "w", **(profile | dict(count=1))))
    classes.set_band_description(1, "Operational baseline coherence recovery classes")
    fraction.set_band_description(1, "Baseline coherence recovery fraction (unclipped)")
    classes.write_colormap(1, COLORS)
    tags = dict(method=json.dumps(settings(c)), zero_is_invalid=str(c.ZERO_IS_INVALID),
                boundary=c.Study_Area,correction="small-hole nearest-neighbor; see correction flags")
    classes.update_tags(**tags, code_0="NoData - missing reference/event or outside mask",
                        **{f"code_{k}": v for k, v in LABELS.items()})
    fraction.update_tags(**tags, units="dimensionless",
                         support="Only assessed detectable-loss pixels (classes 1,2,3)")
    legend = [dict(code=code, label=label, rgba=COLORS[code]) for code, label in LABELS.items()]
    (out / "baseline_recovery_legend.json").write_text(
        json.dumps(dict(classes=legend, nodata=0, settings=settings(c)), indent=2), encoding="utf-8")
    return classes, fraction

def verify(out, dates, source_path, ref, obs, c):
    class_path, fraction_path = paths(out, dates)
    counts = np.zeros(6, dtype="int64")
    with rasterio.open(class_path) as cls, rasterio.open(fraction_path) as frac, rasterio.open(source_path) as source:
        for src, dtype, nodata in ((cls, "uint8", 0), (frac, "float32", c.NODATA)):
            if (src.crs != ref.crs or src.transform != ref.transform or src.shape != ref.shape
                    or src.count != 1 or src.dtypes != (dtype,) or src.nodata != nodata):
                raise RuntimeError("Recovery export grid/dtype/NoData verification failed.")
        for window in windows(ref.width, ref.height, c.WINDOW_SIZE):
            stack = source.read(list(range(1, len(obs)+1)), window=window)
            stack[stack == c.NODATA] = np.nan
            expected_codes, expected_fraction = baseline_recovery(stack, obs, c)
            codes = cls.read(1, window=window)
            fraction = frac.read(1, window=window)
            expected_fraction = np.where(np.isfinite(expected_fraction), expected_fraction, c.NODATA)
            if not np.array_equal(codes, expected_codes) or not np.allclose(fraction, expected_fraction, atol=1e-6):
                raise RuntimeError("Recovery output does not match observed-band classification.")
            counts += np.bincount(codes.ravel(), minlength=6)
    return {str(i): int(n) for i, n in enumerate(counts)}

def finalize(out, dates):
    for path in paths(out, dates):
        target = path.with_name(path.name.replace(".partial.tif", ".tif"))
        path.rename(target)
        write_qgis_style(target, categorical="Classes" in target.name)



def write_qgis_style(path, categorical):
    """Sidecar style auto-loaded by QGIS when it has the same basename as the TIFF."""
    from xml.etree.ElementTree import Element, SubElement, ElementTree
    root=Element("qgis",version="3.34",styleCategories="AllStyleCategories")
    pipe=SubElement(root,"pipe")
    renderer=SubElement(pipe,"rasterrenderer",type="paletted" if categorical else "singlebandpseudocolor",
                        band="1",opacity="1",alphaBand="-1",nodataColor="")
    if categorical:
        palette=SubElement(renderer,"colorPalette")
        for code,label in LABELS.items():
            r,g,b,_=COLORS[code]
            SubElement(palette,"paletteEntry",value=str(code),label=label,
                       color=f"#{r:02x}{g:02x}{b:02x}",alpha="255")
    else:
        renderer.set("classificationMin","-1")
        renderer.set("classificationMax","2")
        shader=SubElement(renderer,"rastershader")
        ramp=SubElement(shader,"colorrampshader",colorRampType="INTERPOLATED",
                        classificationMode="1",minimumValue="-1",maximumValue="2",clip="0")
        for value,color,label in ((-1,"#b2182b","Further decline"),(0,"#ef8a62","No net loss regained"),
                                   (1,"#67a9cf","Baseline level"),(2,"#2166ac","Above baseline")):
            SubElement(ramp,"item",value=str(value),color=color,alpha="255",label=label)
    SubElement(pipe,"rasterresampler",maxOversampling="2",
               zoomedInResampler="nearestNeighbour" if categorical else "bilinear",
               zoomedOutResampler="nearestNeighbour" if categorical else "bilinear")
    ElementTree(root).write(path.with_suffix(".qml"),encoding="utf-8",xml_declaration=True)
