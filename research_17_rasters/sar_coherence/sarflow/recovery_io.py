"""GeoTIFF creation, categorical palette and deterministic class-area accounting."""
import json
import numpy as np
import rasterio
from .spatial import export_profile
from .recovery import CLASSES, QUALITY, REASONS

COLORS = {0: (0, 0, 0, 0), 1: (22, 124, 128, 255), 2: (244, 201, 93, 255),
          3: (199, 58, 73, 255), 4: (216, 221, 226, 255), 5: (140, 120, 163, 255)}
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


def create_outputs(context, folder, ref, c, observations=None):
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
    if observations is not None:
        definitions.update({
            "Original_ITSA_Recovery_Classification": ("uint8", 0, 1),
            "Enhanced_ITSA_Recovery_Classification": ("uint8", 0, 1),
            "Spatial_Recovery_Regions": ("uint32", 0, 1),
            "Evidence_Confidence": ("float32", c.NODATA, 1),
            "Original_Evidence_Confidence": ("float32", c.NODATA, 1),
            "Classification_Changes": ("uint8", 255, 3),
            "Valid_Data_Coverage": ("uint8", 255, len(observations)),
            "Aligned_Coherence_Observations": ("float32", c.NODATA, len(observations)),
            "Raw_ITSA_Metrics": ("float32", c.NODATA, len(METRICS)),
            "Raw_ITSA_Reasons": ("uint8", 255, 1),
            "Raw_ITSA_Quality_Flags": ("uint16", 65535, 1),
            "Enhanced_ITSA_Quality_Flags": ("uint16", 65535, 1),
            "Enhanced_ITSA_Reasons": ("uint8", 255, 1),
            "Spatial_Strata": ("uint32", 0, 1),
            "Builtup_Status": ("int16", -9999, 1),
            "Barangay_Owner": ("uint16", 0, 1),
            "ITSA_Disturbance_Magnitude": ("float32", c.NODATA, 1),
            "ITSA_Recovery_Magnitude": ("float32", c.NODATA, 1),
            "ITSA_Recovery_Fraction": ("float32", c.NODATA, 1),
            "Enhanced_ITSA_Metrics": ("float32", c.NODATA, len(METRICS))})
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
        if name.endswith(("Classes", "Classification")):
            ds.write_colormap(1, COLORS)
            ds.update_tags(class_definitions=json.dumps(CLASSES))
            if name.startswith("Legacy"):
                ds.update_tags(method="legacy observed PRE/event/POST comparison; not ITSA",
                               class_definitions="1 recovered; 2 partial; 3 not recovered; 4 no loss; 0 unavailable")
        for band in range(1, count+1):
            key = METRICS[band-1] if name.endswith("ITSA_Metrics") else name
            if name in ("Valid_Data_Coverage", "Aligned_Coherence_Observations"):
                o = observations[band-1]
                key = f"{o.name}_{o.start:%Y%m%d}_{o.end:%Y%m%d}"
                ds.update_tags(band, master=str(o.master), slave=str(o.slave), phase=o.phase,
                               units="0 missing / 1 observed" if name == "Valid_Data_Coverage" else "coherence")
            if name == "Classification_Changes":
                key = ("original_class", "enhanced_class", "change_reason")[band-1]
            ds.set_band_description(band, key)
            if name not in ("Valid_Data_Coverage", "Aligned_Coherence_Observations"):
                ds.update_tags(band, units=unit(key) if name.endswith("ITSA_Metrics") else "code")
        if name in ("Evidence_Confidence", "Original_Evidence_Confidence"):
            ds.update_tags(score="min(temporal completeness, 1/(1+uncertainty/effect_scale)); not probability")
            ds.update_tags(1, units='dimensionless score [0,1]; not probability')
        if name == 'Builtup_Status':
            ds.update_tags(1, units='-1 unknown; 0 externally mapped non-built-up; 1 externally mapped built-up')
        if name in ('Spatial_Recovery_Regions','Spatial_Object_ID','Regularized_Object_ID'):
            ds.update_tags(1, units='analytical object identifier; not a recovery class')
        if name == "Classification_Changes":
            ds.update_tags(reasons="0 unchanged; 1 own pixel ITSA plus connected neighborhood; 255 NoData")
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
