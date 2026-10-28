"""GeoTIFF creation, categorical palette and deterministic class-area accounting."""
import json
import numpy as np
import rasterio
from .spatial import export_profile
from .recovery import CLASSES, QUALITY, REASONS

COLORS = {0: (0, 0, 0, 0), 1: (0, 158, 115, 255), 2: (230, 171, 2, 255),
          3: (213, 94, 0, 255), 4: (0, 114, 178, 255), 5: (153, 153, 153, 255)}
METRICS = ("baseline", "modeled_level_loss", "observed_event_loss", "post_slope",
           "modeled_initial", "modeled_final", "baseline_deficit", "modeled_progress",
           "recovery_fraction", "crossing_days", "adaptive_loss_threshold", "rmse",
           "level_se", "progress_se", "final_prediction_se", "pre_noise", "pre_count", "post_count", "residual_df", "design_condition",
           "beta0", "beta1", "beta2", "beta3")


def unit(name):
    if name == "crossing_days":
        return "days since EVENT_ORIGIN; observed-period modeled crossing only"
    if name in ("post_slope", "beta1", "beta3"):
        return "coherence/day"
    if name in ("pre_count", "post_count", "residual_df"):
        return "count"
    if name == "recovery_fraction":
        return "dimensionless; improvement / modeled first-POST baseline deficit; unclipped"
    if name == "design_condition":
        return "dimensionless; time-scaled regression design condition number"
    return "coherence units"


def create_outputs(context, folder, ref, c):
    outputs = {}
    definitions = {"ITSA_Recovery_Classes": ("uint8", 0, 1),
                   "Raw_Pixel_ITSA_Classes": ("uint8", 0, 1),
                   "Legacy_Observed_Classes": ("uint8", 0, 1),
                   "Object_ITSA_Classes": ("uint8", 0, 1),
                   "Spatially_Regularized_ITSA_Classes": ("uint8", 0, 1),
                   "Spatial_Object_ID": ("uint32", 0, 1),
                   "Regularized_Object_ID": ("uint32", 0, 1),
                   "Pixel_Object_Classification_Change": ("uint8", 255, 1),
                   "Regularization_Classification_Change": ("uint8", 255, 1),
                   "ITSA_Quality_Flags": ("uint16", 65535, 1),
                   "ITSA_Exclusion_Reasons": ("uint8", 255, 1),
                   "ITSA_Metrics": ("float32", c.NODATA, len(METRICS))}
    if not c.EXPORT_REGULARIZED_CLASSIFICATION:
        definitions = {k: v for k, v in definitions.items() if not k.startswith(("Spatially_Regularized", "Regularized_", "Regularization_"))}
    for name, (dtype, nodata, count) in definitions.items():
        profile = export_profile(ref, c)
        profile.update(count=count, dtype=dtype, nodata=nodata, predictor=3 if dtype == "float32" else 2)
        ds = context.enter_context(rasterio.open(folder / f"{name}.tif", "w", **profile))
        outputs[name] = ds
        ds.update_tags(method="ITSA fixed PRE median recovery", time_origin=c.EVENT_ORIGIN,
                       reference_date=c.REFERENCE_DATE, evidence="operational, exploratory pending validation",
                       classification_mode=c.CLASSIFICATION_MODE, source_validity_confirmed=c.SOURCE_VALIDITY_CONFIRMED,
                       source_validity_notes=c.SOURCE_VALIDITY_NOTES, spatial_filling="none",
                       quality_flags=json.dumps(QUALITY), reasons=json.dumps(REASONS))
        if name.endswith("Classes"):
            ds.write_colormap(1, COLORS)
            ds.update_tags(class_definitions=json.dumps(CLASSES))
            if name.startswith("Legacy"):
                ds.update_tags(method="legacy observed PRE/event/POST comparison; not ITSA",
                               class_definitions="1 recovered; 2 partial; 3 not recovered; 4 no loss; 0 unavailable")
        for band in range(1, count+1):
            key = METRICS[band-1] if name == "ITSA_Metrics" else name
            ds.set_band_description(band, key)
            ds.update_tags(band, units=unit(key) if name == "ITSA_Metrics" else "code")
    return outputs


def write(ds, values, window):
    a = np.asarray(values)
    a = np.where(np.isfinite(a), a, ds.nodata).astype(ds.dtypes[0])
    if a.ndim == 2:
        ds.write(a, 1, window=window)
    else:
        ds.write(a, window=window)


def class_rows(name, counts, inside_count, pixel_area):
    return [dict(product=name, code=code, label=CLASSES[code], pixels=int(n),
                 area_m2=float(n*pixel_area), percent_of_study_pixels=100*n/inside_count if inside_count else None,
                 denominator="all pixel centers within study mask; includes structural NoData")
            for code, n in enumerate(counts)]
